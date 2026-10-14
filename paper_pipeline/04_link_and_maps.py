"""Stage 4: station maps, complexity figures and the FDI-complexity linkage."""

import itertools
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shapefile
from matplotlib.collections import PolyCollection
from matplotlib.lines import Line2D
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, set_link_color_palette

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (LAND_USE_COLORS, LAND_USE_MARKERS, LAND_USE_ORDER, LOCATION_COLORS, LOCATION_ORDER,
                    POLLUTANT_COLORS, POLLUTANT_LABEL, POLLUTANT_ORDER, SEASON_ORDER, SHAPEFILE, TAB_DIR,
                    save_figure, set_style)

CONUS = dict(lon=(-125.5, -66.0), lat=(24.0, 50.0))
METRIC_LABEL = {
    "permutation_entropy": r"Permutation entropy $H$",
    "statistical_complexity": r"Statistical complexity $C_{JS}$",
    "sample_entropy": "Sample entropy",
    "svd_entropy": "SVD entropy",
    "higuchi_fractal_dimension": "Higuchi fractal dimension",
    "correlation_dimension": "Correlation dimension",
    "skewness": "Skewness",
    "kurtosis": "Excess kurtosis",
}
SHORT = {"permutation_entropy": "PE", "statistical_complexity": r"$C_{JS}$", "sample_entropy": "SampEn",
         "svd_entropy": "SVDEn", "higuchi_fractal_dimension": "HFD", "correlation_dimension": r"$D_2$",
         "skewness": "Skew", "kurtosis": "Kurt"}


def state_polygons():
    reader = shapefile.Reader(str(SHAPEFILE))
    polys = []
    for shape in reader.shapes():
        pts = np.asarray(shape.points)
        parts = list(shape.parts) + [len(pts)]
        for a, b in zip(parts[:-1], parts[1:]):
            ring = pts[a:b]
            if ring[:, 0].max() < CONUS["lon"][0] or ring[:, 0].min() > CONUS["lon"][1]:
                continue
            if ring[:, 1].max() < CONUS["lat"][0] or ring[:, 1].min() > CONUS["lat"][1]:
                continue
            polys.append(ring)
    return polys


POLYS = None


def basemap(ax):
    global POLYS
    if POLYS is None:
        POLYS = state_polygons()
    ax.add_collection(PolyCollection(POLYS, facecolor="#eef1f5", edgecolor="#9aa3ad", linewidth=0.35, zorder=0))
    ax.set_xlim(*CONUS["lon"])
    ax.set_ylim(*CONUS["lat"])
    ax.set_aspect(1.25)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def in_conus(frame):
    return frame["Longitude"].between(*CONUS["lon"]) & frame["Latitude"].between(*CONUS["lat"])


def fig_station_map(fdi_st, cx_st):
    set_style()
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 8.6))
    ax = axes[0]
    basemap(ax)
    for lu in LAND_USE_ORDER:
        sub = fdi_st[(fdi_st["Land Use"] == lu) & in_conus(fdi_st)]
        ax.scatter(sub["Longitude"], sub["Latitude"], s=34 if lu != "Forest" else 46,
                   marker=LAND_USE_MARKERS[lu], c=LAND_USE_COLORS[lu], edgecolor="k" if lu != "Forest" else "w",
                   linewidth=0.4, label=f"{lu} ({(fdi_st['Land Use'] == lu).sum()})", zorder=3)
    ax.legend(ncol=4, loc="lower center", bbox_to_anchor=(0.5, -0.17), frameon=False, handletextpad=0.2,
              columnspacing=0.8)
    outside = (~in_conus(fdi_st)).sum()
    ax.set_title(f"(a) Fine Dominance Index network, n = {len(fdi_st)} stations"
                 + (f" ({outside} outside the map extent)" if outside else ""), loc="left")
    ax = axes[1]
    basemap(ax)
    union = cx_st.drop_duplicates("site_id")
    n_par = cx_st.groupby("site_id")["parameter"].nunique()
    union = union.assign(n_par=union["site_id"].map(n_par))
    for loc in LOCATION_ORDER:
        sub = union[(union["Location Setting"] == loc) & in_conus(union)]
        ax.scatter(sub["Longitude"], sub["Latitude"], s=8 + 10 * sub["n_par"], c=LOCATION_COLORS[loc],
                   edgecolor="k", linewidth=0.3, alpha=0.85, label=f"{loc} ({(union['Location Setting'] == loc).sum()})",
                   zorder=3)
    handles = [Line2D([], [], ls="", marker="o", mfc=LOCATION_COLORS[l], mec="k", mew=0.3, ms=6,
                      label=f"{l} ({(union['Location Setting'] == l).sum()})") for l in LOCATION_ORDER]
    handles += [Line2D([], [], ls="", marker="o", mfc="white", mec="k", ms=np.sqrt(8 + 10 * k),
                       label=f"{k} parameter{'s' if k > 1 else ''}") for k in [1, 2, 4]]
    ax.legend(handles=handles, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.17), frameon=False)
    outside = (~in_conus(union)).sum()
    ax.set_title(f"(b) Temporal-complexity network, n = {len(union)} stations"
                 + (f" ({outside} outside the map extent)" if outside else ""), loc="left")
    fig.tight_layout(h_pad=2.5)
    save_figure(fig, "fig_station_map")
    plt.close(fig)
    return dict(fdi_outside=int((~in_conus(fdi_st)).sum()), cx_outside=int((~in_conus(union)).sum()),
                cx_union=int(len(union)))


def cecp_bounds(m, n_points=400):
    """Minimum and maximum statistical-complexity curves for m ordinal states (Martin et al., 2006)."""
    def hc(p):
        p = np.asarray(p, float)
        p = p / p.sum()
        q = 0.5 * (p + 1.0 / m)
        sh = lambda v: -np.sum(v[v > 0] * np.log(v[v > 0]))
        s_p = sh(p)
        h = s_p / np.log(m)
        jsd = sh(q) - 0.5 * s_p - 0.5 * np.log(m)
        q0 = -2.0 / (((m + 1.0) / m) * np.log(m + 1.0) - 2.0 * np.log(2.0 * m) + np.log(m))
        return h, q0 * jsd * h
    cmin = [hc([p] + [(1 - p) / (m - 1)] * (m - 1)) for p in np.linspace(1.0 / m, 1, n_points)]
    cmax = []
    for k in range(m - 1):
        for p in np.linspace(0, 1.0 / (m - k), n_points):
            probs = [0.0] * k + [p] + [(1 - p) / (m - k - 1)] * (m - k - 1)
            cmax.append(hc(probs))
    cmax = np.array(cmax)
    bins = np.linspace(0, 1, 600)
    idx = np.digitize(cmax[:, 0], bins)
    env = np.array([[cmax[idx == i, 0].mean(), cmax[idx == i, 1].max()] for i in np.unique(idx)])
    return np.array(cmin), env[np.argsort(env[:, 0])]


def fig_cecp(cecp):
    set_style()
    cmin, cmax = cecp_bounds(24, n_points=3000)
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.6))
    for ax, mode in zip(axes, ["pollutant", "landuse"]):
        ax.plot(cmin[:, 0], cmin[:, 1], color="0.35", lw=0.9)
        ax.plot(cmax[:, 0], cmax[:, 1], color="0.35", lw=0.9)
        if mode == "pollutant":
            for p in POLLUTANT_ORDER:
                sub = cecp[cecp["parameter"] == p]
                ax.scatter(sub["H"], sub["C"], s=9, alpha=0.6, color=POLLUTANT_COLORS[p], edgecolor="none",
                           label=f"{POLLUTANT_LABEL[p]} ({len(sub)})")
            ax.scatter(cecp["H_surr"].median(), cecp["C_surr"].median(), marker="*", s=90, color="k",
                       label="Shuffled surrogates", zorder=4)
            ax.set_title(r"(a) All parameters, $D=4$, $\tau=1$ day", loc="left")
        else:
            sub = cecp[cecp["parameter"] == "PM10"]
            for lu in LAND_USE_ORDER:
                s2 = sub[sub["Land Use"] == lu]
                if s2.empty:
                    continue
                ax.scatter(s2["H"], s2["C"], s=14, marker=LAND_USE_MARKERS[lu], color=LAND_USE_COLORS[lu],
                           alpha=0.75, edgecolor="none", label=f"{lu} ({len(s2)})")
            ax.set_title(r"(b) PM$_{10}$ stations by land use", loc="left")
        ax.set_xlim(0.925, 1.0)
        ax.set_ylim(0.0, 0.085)
        ax.set_xlabel(r"Normalized permutation entropy $H$")
        ax.set_ylabel(r"Statistical complexity $C_{JS}$")
        ax.legend(frameon=False, fontsize=6.5, loc="upper right")
        ax.grid(alpha=0.25)
        inset = ax.inset_axes([0.07, 0.07, 0.33, 0.33])
        inset.plot(cmin[:, 0], cmin[:, 1], color="0.35", lw=0.7)
        inset.plot(cmax[:, 0], cmax[:, 1], color="0.35", lw=0.7)
        inset.scatter(cecp["H"], cecp["C"], s=1, color="0.2")
        inset.add_patch(plt.Rectangle((0.925, 0), 0.075, 0.085, fill=False, ec="red", lw=0.7))
        inset.set_xlim(0, 1)
        inset.set_ylim(0, 0.5)
        inset.tick_params(labelsize=5.5, length=2)
    fig.tight_layout()
    save_figure(fig, "fig_cecp")
    plt.close(fig)


def fig_seasonal(site_season):
    set_style()
    metrics = ["permutation_entropy", "sample_entropy", "higuchi_fractal_dimension", "kurtosis"]
    fig, axes = plt.subplots(2, 2, figsize=(7.4, 5.6))
    width = 0.19
    for ax, metric in zip(axes.flat, metrics):
        for k, p in enumerate(POLLUTANT_ORDER):
            sub = site_season[site_season["parameter"] == p]
            data = [sub.loc[sub["Season"] == s, metric].dropna() for s in SEASON_ORDER]
            pos = np.arange(4) + (k - 1.5) * width
            bp = ax.boxplot(data, positions=pos, widths=width * 0.9, whis=(5, 95), showfliers=False,
                            patch_artist=True, medianprops=dict(color="white", lw=1.1))
            for b in bp["boxes"]:
                b.set(facecolor=POLLUTANT_COLORS[p], edgecolor="k", lw=0.4, alpha=0.9)
        ax.set_xticks(range(4))
        ax.set_xticklabels(SEASON_ORDER)
        ax.set_ylabel(METRIC_LABEL[metric])
        ax.grid(axis="y", alpha=0.3)
        if metric == "kurtosis":
            ax.set_ylim(-1.5, 14)
    for ax, letter in zip(axes.flat, "abcd"):
        ax.set_title(f"({letter})", loc="left")
    handles = [plt.Rectangle((0, 0), 1, 1, fc=POLLUTANT_COLORS[p], ec="k", lw=0.4, label=POLLUTANT_LABEL[p])
               for p in POLLUTANT_ORDER]
    fig.legend(handles=handles, ncol=4, loc="lower center", frameon=False, bbox_to_anchor=(0.5, -0.03))
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    save_figure(fig, "fig_complexity_seasonal")
    plt.close(fig)


def fig_landuse(site, kw):
    set_style()
    metrics = ["permutation_entropy", "sample_entropy", "kurtosis"]
    fig, axes = plt.subplots(2, 3, figsize=(7.4, 5.4), sharex=False)
    for i, p in enumerate(["PM2.5 FRM/FEM", "PM10"]):
        sub = site[site["parameter"] == p]
        counts = sub.groupby("Land Use")["site_id"].nunique()
        groups = [g for g in LAND_USE_ORDER if counts.get(g, 0) >= 5]
        for j, metric in enumerate(metrics):
            ax = axes[i, j]
            data = [sub.loc[sub["Land Use"] == g, metric].dropna() for g in groups]
            bp = ax.boxplot(data, widths=0.65, whis=(5, 95), showfliers=False, patch_artist=True,
                            medianprops=dict(color="white", lw=1.1))
            for b, g in zip(bp["boxes"], groups):
                b.set(facecolor=LAND_USE_COLORS[g], edgecolor="k", lw=0.4, alpha=0.85)
            ax.set_xticks(range(1, len(groups) + 1))
            ax.set_xticklabels([f"{g} ({counts[g]})" for g in groups], rotation=60, ha="right", fontsize=7)
            row = kw[(kw["parameter"] == p) & (kw["factor"] == "Land Use") & (kw["metric"] == metric)].iloc[0]
            ptxt = "< 0.001" if row["p_holm"] < 0.001 else f"= {row['p_holm']:.3f}"
            ax.set_title(f"({'abcdef'[3 * i + j]}) {POLLUTANT_LABEL[p]}", loc="left")
            ax.text(0.03, 0.97, r"$\varepsilon^2$" + f" = {row['eps2']:.3f}\n" + r"$p_{\mathrm{Holm}}$ " + ptxt,
                    transform=ax.transAxes, ha="left", va="top", fontsize=7,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.5), zorder=5)
            ax.set_ylabel(METRIC_LABEL[metric])
            ax.grid(axis="y", alpha=0.3)
            if metric == "kurtosis":
                ax.set_ylim(-1, 16)
    fig.tight_layout()
    save_figure(fig, "fig_complexity_landuse")
    plt.close(fig)


def fig_maps(site, clusters):
    set_style()
    fig, axes = plt.subplots(3, 1, figsize=(6.6, 11.2))
    for ax, p, letter in zip(axes[:2], ["PM2.5 FRM/FEM", "PM10"], "ab"):
        basemap(ax)
        sub = site[(site["parameter"] == p) & in_conus(site)].sort_values("permutation_entropy")
        lo, hi = site.loc[site["parameter"].isin(["PM2.5 FRM/FEM", "PM10"]), "permutation_entropy"].quantile([0.05, 0.95])
        sc = ax.scatter(sub["Longitude"], sub["Latitude"], c=sub["permutation_entropy"], cmap="viridis",
                        vmin=lo, vmax=hi, s=16, edgecolor="k", linewidth=0.25, zorder=3)
        cb = fig.colorbar(sc, ax=ax, fraction=0.03, pad=0.01)
        cb.set_label("Median permutation entropy $H$ ($D=3$)")
        ax.set_title(f"({letter}) {POLLUTANT_LABEL[p]}, {len(site[site['parameter'] == p])} stations", loc="left")
    ax = axes[2]
    basemap(ax)
    palette = {1: "#1f77b4", 2: "#d62728", 3: "#2ca02c"}
    sub = clusters[in_conus(clusters)]
    for c in sorted(clusters["cluster"].unique()):
        s2 = sub[sub["cluster"] == c]
        ax.scatter(s2["Longitude"], s2["Latitude"], s=14, color=palette[c], edgecolor="k", linewidth=0.25,
                   label=f"Cluster {c} ({(clusters['cluster'] == c).sum()})", zorder=3)
    ax.legend(ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.12), frameon=False)
    ax.set_title(r"(c) Ward clusters of PM$_{2.5}$ FRM/FEM dynamical profiles ($k=3$)", loc="left")
    fig.tight_layout(h_pad=2.0)
    save_figure(fig, "fig_complexity_maps")
    plt.close(fig)


def fig_structure(corr, link, clusters):
    set_style()
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.6), gridspec_kw=dict(width_ratios=[1, 1.15]))
    ax = axes[0]
    order = list(METRIC_LABEL)
    c = corr.loc[order, order].to_numpy()
    mask = np.triu(np.ones_like(c, bool), k=1)
    im = ax.imshow(np.ma.masked_where(mask, c), cmap="RdBu_r", vmin=-1, vmax=1)
    for i, j in itertools.product(range(len(order)), repeat=2):
        if j <= i:
            ax.text(j, i, f"{c[i, j]:.2f}", ha="center", va="center", fontsize=6.2,
                    color="white" if abs(c[i, j]) > 0.6 else "black")
    ax.set_xticks(range(len(order)))
    ax.set_yticks(range(len(order)))
    ax.set_xticklabels([SHORT[m] for m in order], rotation=60, ha="right")
    ax.set_yticklabels([SHORT[m] for m in order])
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02, label=r"Spearman $\rho$")
    ax.set_title("(a) Metric redundancy (all windows)", loc="left")
    ax = axes[1]
    palette = {1: "#1f77b4", 2: "#d62728", 3: "#2ca02c"}
    leaves = dendrogram(link, no_plot=True)["leaves"]
    sequence = list(dict.fromkeys(clusters["cluster"].to_numpy()[leaves]))
    set_link_color_palette([palette[c] for c in sequence])
    dendrogram(link, ax=ax, no_labels=True, color_threshold=link[-(3 - 1), 2], above_threshold_color="0.4")
    set_link_color_palette(None)
    ax.legend(handles=[Line2D([], [], color=palette[c], lw=2,
                              label=f"Cluster {c} ({(clusters['cluster'] == c).sum()})") for c in sorted(palette)],
              frameon=False, loc="upper right", fontsize=7)
    ax.set_ylabel("Ward distance")
    ax.set_title(r"(b) Ward dendrogram, PM$_{2.5}$ FRM/FEM stations", loc="left")
    ax.spines["bottom"].set_visible(False)
    fig.tight_layout()
    save_figure(fig, "fig_complexity_structure")
    plt.close(fig)


def partial_spearman(x, y, covars):
    rx = stats.rankdata(x)
    ry = stats.rankdata(y)
    z = np.column_stack([np.ones(len(x))] + [stats.rankdata(c) for c in covars])
    res_x = rx - z @ np.linalg.lstsq(z, rx, rcond=None)[0]
    res_y = ry - z @ np.linalg.lstsq(z, ry, rcond=None)[0]
    r = np.corrcoef(res_x, res_y)[0, 1]
    n, k = len(x), len(covars)
    t = r * np.sqrt((n - 2 - k) / (1 - r ** 2))
    return r, 2 * stats.t.sf(abs(t), n - 2 - k)


def link(fdi_means, site, cecp):
    rows, data = [], {}
    for p in ["PM2.5 FRM/FEM", "PM10"]:
        sub = site[site["parameter"] == p].set_index("site_id")
        full = cecp[cecp["parameter"] == p].set_index("site_id")[["H", "C"]].add_suffix("_D4")
        merged = fdi_means.set_index("site_id").join(sub[list(METRIC_LABEL)], how="inner").join(full, how="left")
        data[p] = merged
        for metric in list(METRIC_LABEL) + ["H_D4", "C_D4"]:
            ok = merged[["FDI", metric, "Latitude", "Longitude"]].dropna()
            rho, pv = stats.spearmanr(ok["FDI"], ok[metric])
            prho, ppv = partial_spearman(ok["FDI"].to_numpy(), ok[metric].to_numpy(),
                                         [ok["Latitude"].to_numpy(), ok["Longitude"].to_numpy()])
            rows.append(dict(parameter=p, metric=metric, stations=len(ok), rho=rho, p=pv, partial_rho=prho, partial_p=ppv))
    out = pd.DataFrame(rows)
    from importlib import import_module  # noqa: F401  (keeps the module free of heavy imports at top)
    p = out["p"].to_numpy()
    order = np.argsort(p)
    adj = np.empty_like(p)
    running = 0.0
    for rank, idx in enumerate(order):
        running = max(running, (len(p) - rank) * p[idx])
        adj[idx] = min(running, 1.0)
    out["p_holm"] = adj
    out.to_csv(TAB_DIR / "link_fdi_complexity.csv", index=False)
    print(out.round(4).to_string())
    return out, data


def fig_link(data, out):
    set_style()
    panels = [("PM2.5 FRM/FEM", "H_D4"), ("PM10", "kurtosis"), ("PM10", "sample_entropy")]
    labels = {"H_D4": r"Permutation entropy $H$ ($D=4$)", "kurtosis": "Excess kurtosis (median)",
              "sample_entropy": "Sample entropy (median)"}
    fig, axes = plt.subplots(1, 3, figsize=(7.6, 2.9))
    for ax, (p, metric), letter in zip(axes, panels, "abc"):
        d = data[p]
        for lu in LAND_USE_ORDER:
            s = d[d["Land Use"] == lu]
            ax.scatter(s["FDI"], s[metric], s=13, marker=LAND_USE_MARKERS[lu], color=LAND_USE_COLORS[lu],
                       alpha=0.8, edgecolor="none", label=lu)
        row = out[(out["parameter"] == p) & (out["metric"] == metric)].iloc[0]
        ptxt = "< 0.001" if row["p_holm"] < 0.001 else f"= {row['p_holm']:.3f}"
        ax.text(0.97, 0.96, rf"$\rho_s$ = {row['rho']:.2f} ($p_{{\mathrm{{Holm}}}}$ {ptxt})" + "\n"
                + rf"partial $\rho_s$ = {row['partial_rho']:.2f}",
                transform=ax.transAxes, ha="right", va="top", fontsize=6.8,
                bbox=dict(facecolor="white", edgecolor="0.75", alpha=0.9, pad=2), zorder=5)
        ax.set_xlabel("Station mean FDI")
        ax.set_ylabel(labels[metric])
        ax.set_title(f"({letter}) {POLLUTANT_LABEL[p]}", loc="left")
        ax.grid(alpha=0.25)
        if metric == "kurtosis":
            ax.set_yscale("symlog", linthresh=2)
            ax.set_yticks([0, 1, 2, 5, 10, 20])
            ax.set_yticklabels(["0", "1", "2", "5", "10", "20"])
    handles, labs = axes[0].get_legend_handles_labels()
    fig.legend(handles, labs, ncol=7, loc="lower center", frameon=False, bbox_to_anchor=(0.5, -0.06), fontsize=7)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    save_figure(fig, "fig_fdi_complexity_link")
    plt.close(fig)


def paired_d4(cecp):
    rows = []
    for a, b in [("PM2.5 FRM/FEM", "PM10"), ("PM10", "PMc"), ("PM2.5 FRM/FEM", "PMc"),
                 ("PM2.5 FRM/FEM", "PM2.5 non-FRM/FEM")]:
        wa = cecp[cecp["parameter"] == a].set_index("site_id")
        wb = cecp[cecp["parameter"] == b].set_index("site_id")
        common = wa.index.intersection(wb.index)
        for metric in ["H", "C"]:
            x, y = wa.loc[common, metric], wb.loc[common, metric]
            d = x - y
            ranks = stats.rankdata(np.abs(d))
            rows.append(dict(a=a, b=b, metric=metric, sites=len(common), median_a=x.median(), median_b=y.median(),
                             p=stats.wilcoxon(x, y).pvalue,
                             rank_biserial=(ranks[d > 0].sum() - ranks[d < 0].sum()) / ranks.sum()))
    out = pd.DataFrame(rows)
    out.to_csv(TAB_DIR / "cecp_paired_pollutants.csv", index=False)
    kw = []
    for p in POLLUTANT_ORDER:
        sub = cecp[cecp["parameter"] == p]
        for factor, order in [("Land Use", LAND_USE_ORDER), ("Location Setting", LOCATION_ORDER)]:
            n_g = sub.groupby(factor).size()
            groups = [g for g in order if n_g.get(g, 0) >= 5]
            if len(groups) < 2:
                continue
            data = [sub.loc[sub[factor] == g, "H"] for g in groups]
            h, pv = stats.kruskal(*data)
            kw.append(dict(parameter=p, factor=factor, groups=len(groups), sites=sum(map(len, data)), H=h, p=pv,
                           eps2=h / (sum(map(len, data)) - 1),
                           medians="; ".join(f"{g}: {np.median(d):.4f}" for g, d in zip(groups, data))))
    pd.DataFrame(kw).to_csv(TAB_DIR / "cecp_kruskal.csv", index=False)
    print(out.round(4).to_string())
    print(pd.DataFrame(kw).round(4).to_string())


def main():
    fdi_st = pd.read_csv(TAB_DIR / "fdi_network_stations.csv")
    cx_st = pd.read_csv(TAB_DIR / "complexity_network_stations.csv")
    site = pd.read_csv(TAB_DIR / "complexity_site_medians.csv")
    site_season = pd.read_csv(TAB_DIR / "complexity_site_season_medians.csv")
    kw = pd.read_csv(TAB_DIR / "complexity_kruskal_landuse_location.csv")
    clusters = pd.read_csv(TAB_DIR / "complexity_clusters_pm25.csv")
    link_matrix = np.load(TAB_DIR / "ward_linkage_pm25.npy")
    corr = pd.read_csv(TAB_DIR / "complexity_spearman_all.csv", index_col=0)
    cecp = pd.read_csv(TAB_DIR / "cecp_full_period_D4.csv")
    fdi_means = pd.read_csv(TAB_DIR / "fdi_station_means.csv")

    info = fig_station_map(fdi_st, cx_st)
    print(info)
    pd.Series(info).to_csv(TAB_DIR / "map_info.csv")
    fig_cecp(cecp)
    fig_seasonal(site_season)
    fig_landuse(site, kw)
    fig_maps(site, clusters)
    fig_structure(corr, link_matrix, clusters)
    paired_d4(cecp)
    out, data = link(fdi_means, site, cecp)
    fig_link(data, out)


if __name__ == "__main__":
    main()
