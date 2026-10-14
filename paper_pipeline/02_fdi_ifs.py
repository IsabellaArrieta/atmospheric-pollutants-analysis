"""Stage 2: Fine Dominance Index (FDI), cell fingerprints (IFS), robustness and contrasts.

Inputs : table_csv/daily_long.parquet, table_csv/sites_meta.parquet
Outputs: outputs_paper/tables/*.csv, figures fig_fdi_boxplots and fig_ifs_heatmap
"""

import itertools
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Patch, Rectangle
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (DATA_DIR, LAND_USE_COLORS, LAND_USE_ORDER, LOCATION_ORDER, MIN_VALID_DAYS,
                    SEASON_OF_MONTH, SEASON_ORDER, SEED, TAB_DIR, TOTAL_DAYS, save_figure, set_style)

CELL = ["Land Use", "Location Setting", "Season"]
WEIGHTS = {"FDI": 1.5, "r0": -1.0, "skew": -0.5, "kurt": -0.5}
MIN_CELL_DAYS = 300
MIN_CELL_STATIONS = 3
MAX_CI_WIDTH = 0.15  # station-block bootstrap; see robustness for 0.10 and no CI criterion
MIN_UNIT_DAYS = 45  # station-season-year units need at least half a season of ratio days
MIN_GROUP_STATIONS = 3
B_BOOT = 2000


def load_ratio_days(min_days=MIN_VALID_DAYS, keep_inconsistent=False):
    daily = pd.read_parquet(DATA_DIR / "daily_long.parquet")
    meta = pd.read_parquet(DATA_DIR / "sites_meta.parquet").set_index("site_id")
    wide = (daily[daily["parameter"].isin(["PM2.5 FRM/FEM", "PM10"])]
            .pivot_table(index=["site_id", "date"], columns="parameter", values="value")
            .dropna().reset_index()
            .rename(columns={"PM2.5 FRM/FEM": "PM25"}))
    concurrent = wide.groupby("site_id").size()
    eligible = meta.loc[meta.index.intersection(concurrent[concurrent >= min_days].index)]
    eligible = eligible[eligible["is_us"] & eligible["Land Use"].isin(LAND_USE_ORDER)
                        & eligible["Location Setting"].isin(LOCATION_ORDER)]
    wide = wide[wide["site_id"].isin(eligible.index)].copy()
    n_all = len(wide)
    wide = wide[(wide["PM10"] > 0) & (wide["PM25"] > 0)]
    n_pos = len(wide)
    inconsistent = wide["PM25"] > wide["PM10"]
    if keep_inconsistent:
        wide.loc[inconsistent, "PM25"] = wide.loc[inconsistent, "PM10"]
    else:
        wide = wide[~inconsistent]
    x = wide["PM25"] / wide["PM10"]
    wide["r12"] = x
    wide["r21"] = 1.0 / x
    wide["r0"] = (wide["PM10"] - wide["PM25"]) / (wide["PM10"] + wide["PM25"])
    wide["FDI"] = wide["r12"] * (1.0 - wide["r0"])
    wide["Season"] = wide["date"].dt.month.map(SEASON_OF_MONTH)
    wide["season_year"] = wide["date"].dt.year + (wide["date"].dt.month == 12).astype(int)
    wide = wide.join(eligible[["Land Use", "Location Setting", "Latitude", "Longitude", "State Name"]],
                     on="site_id")
    qc = {"stations": int(wide["site_id"].nunique()), "site_days_concurrent": int(n_all),
          "site_days_positive": int(n_pos), "share_pm25_gt_pm10": float(inconsistent.mean()),
          "ratio_days": int(len(wide)), "min_days": int(min_days)}
    return wide, eligible, qc


def block_bootstrap_mean(frame, value="FDI", b=B_BOOT, seed=SEED):
    """95 % CI of the pooled mean, resampling whole stations (respects serial dependence)."""
    per_station = frame.groupby("site_id")[value].agg(["sum", "count"])
    sums, counts = per_station["sum"].to_numpy(), per_station["count"].to_numpy()
    k = len(sums)
    rng = np.random.default_rng(seed)
    if k < 2:
        return np.nan, np.nan
    draws = rng.integers(0, k, size=(b, k))
    boot = sums[draws].sum(axis=1) / counts[draws].sum(axis=1)
    return tuple(np.quantile(boot, [0.025, 0.975]))


def cell_summary(ratio, max_ci_width=MAX_CI_WIDTH):
    rows = []
    for keys, g in ratio.groupby(CELL):
        lo, hi = block_bootstrap_mean(g)
        rows.append(dict(zip(CELL, keys), n=len(g), n_st=g["site_id"].nunique(),
                         FDI_mean=g["FDI"].mean(), FDI_median=g["FDI"].median(), FDI_lo=lo, FDI_hi=hi,
                         r0_mean=g["r0"].mean(), r0_skew=stats.skew(g["r0"]),
                         r0_kurt=stats.kurtosis(g["r0"]), r12_mean=g["r12"].mean(),
                         PM25_mean=g["PM25"].mean(), PM10_mean=g["PM10"].mean()))
    cells = pd.DataFrame(rows)
    cells["ci_width"] = cells["FDI_hi"] - cells["FDI_lo"]
    cells["thin"] = ((cells["n_st"] < MIN_CELL_STATIONS) | (cells["n"] < MIN_CELL_DAYS)
                     | (cells["ci_width"] > max_ci_width) | cells["ci_width"].isna())
    return cells


FEATURES = {"FDI": "FDI_mean", "r0": "r0_mean", "skew": "r0_skew", "kurt": "r0_kurt"}


def zscores(cells, reference_mask):
    ref = cells.loc[reference_mask, list(FEATURES.values())]
    return (cells[list(FEATURES.values())] - ref.mean()) / ref.std(ddof=0)


def ifs(cells, weights=WEIGHTS):
    """Composite index on the network-wide z-scale (all populated cells are the reference)."""
    z = zscores(cells, np.ones(len(cells), bool))
    return sum(weights[k] * z[col] for k, col in FEATURES.items())


def overlap(a, b, k=10):
    return len(set(a[:k]) & set(b[:k])) / k


def ranked_keys(frame, column, ascending):
    ordered = frame.sort_values(column, ascending=ascending)
    return ordered[CELL].astype(str).agg(" | ".join, axis=1).tolist()


def compare(base, alt_col, frame=None):
    frame = base if frame is None else frame
    return dict(tau=stats.kendalltau(frame["IFS"], frame[alt_col])[0],
                top10=overlap(ranked_keys(frame, "IFS", False), ranked_keys(frame, alt_col, False)),
                bottom10=overlap(ranked_keys(frame, "IFS", True), ranked_keys(frame, alt_col, True)))


def compare_tables(base, alt_cells):
    """Rank agreement on cells retained by both specifications, and overlap of each one's own
    top/bottom-10 league table (thin cells excluded under each specification separately)."""
    alt = alt_cells[~alt_cells["thin"]]
    merged = base[CELL + ["IFS"]].merge(alt[CELL + ["IFS"]], on=CELL, suffixes=("", "_alt"))
    return dict(common_cells=len(merged),
                tau=stats.kendalltau(merged["IFS"], merged["IFS_alt"])[0],
                rho=stats.spearmanr(merged["IFS"], merged["IFS_alt"])[0],
                top10=overlap(ranked_keys(base, "IFS", False), ranked_keys(alt, "IFS", False)),
                bottom10=overlap(ranked_keys(base, "IFS", True), ranked_keys(alt, "IFS", True)))


def robustness(cells, ratio):
    z_all = zscores(cells, np.ones(len(cells), bool))
    mask = ~cells["thin"]
    base = cells[mask].copy()
    z = z_all[mask]

    rows = []
    grid = itertools.product([1.0, 1.25, 1.5, 1.75, 2.0], [-0.5, -0.75, -1.0, -1.25, -1.5],
                             [0.0, -0.25, -0.5, -0.75, -1.0], [0.0, -0.25, -0.5, -0.75, -1.0])
    for wf, wr, ws, wk in grid:
        base["alt"] = wf * z["FDI_mean"] + wr * z["r0_mean"] + ws * z["r0_skew"] + wk * z["r0_kurt"]
        rows.append(dict(w_FDI=wf, w_r0=wr, w_skew=ws, w_kurt=wk, **compare(base, "alt")))
    grid = pd.DataFrame(rows)
    grid.to_csv(TAB_DIR / "ifs_weight_grid.csv", index=False)

    named = {"Equal weights": dict(FDI=1, r0=-1, skew=-1, kurt=-1),
             "Location terms only (no shape terms)": dict(FDI=1.5, r0=-1, skew=0, kurt=0),
             "Shape terms with positive sign": dict(FDI=1.5, r0=-1, skew=0.5, kurt=0.5),
             "Mean FDI only": dict(FDI=1, r0=0, skew=0, kurt=0)}
    named_rows = []
    for name, w in named.items():
        base["alt"] = sum(w[k] * z[c] for k, c in FEATURES.items())
        named_rows.append(dict(spec=name, **compare(base, "alt")))

    # First principal component of the standardized features (sign aligned with IFS).
    zf = z.to_numpy() - z.to_numpy().mean(axis=0)
    _, sing, vt = np.linalg.svd(zf, full_matrices=False)
    pc1 = zf @ vt[0]
    sign = np.sign(stats.spearmanr(pc1, base["IFS"])[0])
    base["alt"] = sign * pc1
    named_rows.append(dict(spec="PC1 of standardized features", **compare(base, "alt"),
                           rho=stats.spearmanr(base["alt"], base["IFS"])[0],
                           pc1_var=float(sing[0] ** 2 / (sing ** 2).sum()),
                           loadings=" ".join(f"{k}:{v:+.2f}" for k, v in zip(FEATURES, sign * vt[0]))))

    # Alternative screening, QC and thin-cell choices, recomputed end to end.
    alternatives = {"Completeness 90%": (dict(min_days=int(round(0.90 * TOTAL_DAYS))), MAX_CI_WIDTH),
                    "Retain PM2.5>PM10 days (r12 capped at 1)": (dict(keep_inconsistent=True), MAX_CI_WIDTH),
                    "Thin-cell CI tolerance 0.10": ({}, 0.10),
                    "No CI criterion for thin cells": ({}, np.inf)}
    for label, (kwargs, tol) in alternatives.items():
        alt_ratio, _, qc_alt = load_ratio_days(**kwargs) if kwargs else (ratio, None, {"stations": ratio["site_id"].nunique()})
        alt_cells = cell_summary(alt_ratio, max_ci_width=tol)
        alt_cells["IFS"] = ifs(alt_cells)
        named_rows.append(dict(spec=label, stations=qc_alt["stations"],
                               non_thin_alt=int((~alt_cells["thin"]).sum()), **compare_tables(base, alt_cells)))

    # Year-by-year replication: IFS recomputed separately for each season-year.
    year_rows = []
    for year in [2023, 2024, 2025]:
        yc = cell_summary(ratio[ratio["season_year"] == year])
        yc["IFS"] = ifs(yc)
        year_rows.append(dict(year=year, non_thin_year=int((~yc["thin"]).sum()), **compare_tables(base, yc)))
    named = pd.DataFrame(named_rows)
    named.to_csv(TAB_DIR / "ifs_robustness_named.csv", index=False)
    years = pd.DataFrame(year_rows)
    years.to_csv(TAB_DIR / "ifs_year_replication.csv", index=False)
    summary = grid[["tau", "top10", "bottom10"]].describe(percentiles=[0.05, 0.5]).T
    summary.to_csv(TAB_DIR / "ifs_weight_grid_summary.csv")
    print(summary)
    print(named.to_string())
    print(years)


def cliffs_delta(x, y):
    x, y = np.asarray(x), np.asarray(y)
    ranks = stats.rankdata(np.concatenate([x, y]))
    u = ranks[:len(x)].sum() - len(x) * (len(x) + 1) / 2
    return 2 * u / (len(x) * len(y)) - 1


def holm(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    adj = np.empty_like(p)
    running = 0.0
    for rank, idx in enumerate(order):
        running = max(running, (len(p) - rank) * p[idx])
        adj[idx] = min(running, 1.0)
    return adj


def contrasts(ratio):
    units = (ratio.groupby(["site_id", "Land Use", "Location Setting", "Season", "season_year"])
             .agg(n=("FDI", "size"), FDI=("FDI", "mean"), r0=("r0", "mean"), r12=("r12", "mean"),
                  PM25=("PM25", "mean"), PM10=("PM10", "mean"))
             .reset_index())
    units = units[(units["n"] >= MIN_UNIT_DAYS) & units["season_year"].isin([2023, 2024, 2025])]
    units.to_csv(TAB_DIR / "fdi_station_season_year_units.csv", index=False)
    omni, post = [], []
    for (loc, sea), sub in units.groupby(["Location Setting", "Season"]):
        n_st = sub.groupby("Land Use")["site_id"].nunique()
        groups = [lu for lu in LAND_USE_ORDER if n_st.get(lu, 0) >= MIN_GROUP_STATIONS]
        if len(groups) < 2:
            continue
        for metric in ["FDI", "r0", "PM25", "PM10"]:
            data = {lu: sub.loc[sub["Land Use"] == lu, metric].to_numpy() for lu in groups}
            if len(groups) >= 3:
                h, p = stats.kruskal(*data.values())
            else:
                h, p = np.nan, stats.mannwhitneyu(*data.values(), alternative="two-sided").pvalue
            n_tot = sum(len(v) for v in data.values())
            omni.append(dict(Location=loc, Season=sea, metric=metric, groups=len(groups), units=n_tot,
                             H=h, p=p, eps2=h / (n_tot - 1) if np.isfinite(h) else np.nan))
            pairs = list(itertools.combinations(groups, 2))
            pvals = [stats.mannwhitneyu(data[a], data[b], alternative="two-sided").pvalue for a, b in pairs]
            for (a, b), pv, pa in zip(pairs, pvals, holm(pvals)):
                post.append(dict(Location=loc, Season=sea, metric=metric, a=a, b=b, n_a=len(data[a]),
                                 n_b=len(data[b]), p=pv, p_holm=pa, delta=cliffs_delta(data[a], data[b])))
    omni = pd.DataFrame(omni)
    omni["p_holm"] = holm(omni["p"].to_numpy())
    post = pd.DataFrame(post)
    omni.to_csv(TAB_DIR / "fdi_kruskal_wallis.csv", index=False)
    post.to_csv(TAB_DIR / "fdi_posthoc_cliffs_delta.csv", index=False)
    print(omni[omni.metric == "FDI"].round(4).to_string())
    print(post[post.metric == "FDI"].round(3).to_string())
    return units


def seasonality(ratio):
    """Within-station seasonal contrast of FDI (Friedman test, Kendall's W) by location setting."""
    st = (ratio.groupby(["site_id", "Location Setting", "Season"])["FDI"].mean()
          .unstack("Season").reindex(columns=SEASON_ORDER).dropna())
    rows = []
    for loc, sub in [("All", st)] + [(loc, st.xs(loc, level="Location Setting")) for loc in LOCATION_ORDER]:
        chi2, p = stats.friedmanchisquare(*[sub[s] for s in SEASON_ORDER])
        n = len(sub)
        w_win_sum = stats.wilcoxon(sub["Winter"], sub["Summer"])
        rows.append(dict(Location=loc, stations=n, chi2=chi2, p=p, kendall_W=chi2 / (n * 3),
                         **{f"median_{s}": sub[s].median() for s in SEASON_ORDER},
                         share_winter_max=float((sub.idxmax(axis=1) == "Winter").mean()),
                         winter_vs_summer_p=w_win_sum.pvalue))
    out = pd.DataFrame(rows)
    out.to_csv(TAB_DIR / "fdi_seasonality_friedman.csv", index=False)
    print(out.round(4).to_string())


def season_station_means(ratio):
    """Seasonal FDI means per station (pooled over years), used for the station-level linkage."""
    st = (ratio.groupby(["site_id", "Land Use", "Location Setting", "Latitude", "Longitude"])
          .agg(FDI=("FDI", "mean"), r0=("r0", "mean"), r12=("r12", "mean"), n=("FDI", "size"))
          .reset_index())
    st.to_csv(TAB_DIR / "fdi_station_means.csv", index=False)


def plot_boxplots(ratio, cells):
    set_style()
    fig, axes = plt.subplots(3, 4, figsize=(11, 8.2), sharey=True, sharex=True)
    for i, loc in enumerate(LOCATION_ORDER):
        for j, sea in enumerate(SEASON_ORDER):
            ax = axes[i, j]
            sub = ratio[(ratio["Location Setting"] == loc) & (ratio["Season"] == sea)]
            data, pos, cols = [], [], []
            for k, lu in enumerate(LAND_USE_ORDER):
                vals = sub.loc[sub["Land Use"] == lu, "FDI"].to_numpy()
                if len(vals) == 0:
                    continue
                data.append(vals)
                pos.append(k)
                cols.append(LAND_USE_COLORS[lu])
            if data:
                bp = ax.boxplot(data, positions=pos, widths=0.7, whis=(5, 95), showfliers=False,
                                patch_artist=True, medianprops=dict(color="white", lw=1.4))
                for patch, c in zip(bp["boxes"], cols):
                    patch.set(facecolor=c, alpha=0.85, edgecolor="black", lw=0.6)
                for k, lu in zip(pos, [LAND_USE_ORDER[p] for p in pos]):
                    row = cells[(cells["Land Use"] == lu) & (cells["Location Setting"] == loc)
                                & (cells["Season"] == sea)]
                    if len(row) and bool(row["thin"].iloc[0]):
                        ax.text(k, 0.98, "†", ha="center", va="top", fontsize=9, color="dimgray",
                                transform=ax.get_xaxis_transform())
            ax.set_xticks(range(len(LAND_USE_ORDER)))
            ax.set_xticklabels(LAND_USE_ORDER, rotation=90)
            ax.grid(axis="y", alpha=0.3)
            if i == 0:
                ax.set_title(sea, fontweight="bold")
            if j == 0:
                ax.set_ylabel(f"{loc}\nFDI")
    axes[0, 0].set_ylim(0, 1.0)
    fig.tight_layout()
    save_figure(fig, "fig_fdi_boxplots")
    plt.close(fig)


def plot_ifs_heatmap(cells):
    set_style()
    vmax = np.nanmax(np.abs(cells["IFS"]))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    cmap = plt.get_cmap("coolwarm_r")
    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.8), sharey=True)
    for ax, loc in zip(axes, LOCATION_ORDER):
        sub = cells[cells["Location Setting"] == loc]
        piv = sub.pivot(index="Land Use", columns="Season", values="IFS").reindex(index=LAND_USE_ORDER, columns=SEASON_ORDER)
        thin = sub.pivot(index="Land Use", columns="Season", values="thin").reindex(index=LAND_USE_ORDER, columns=SEASON_ORDER)
        ax.imshow(np.ma.masked_invalid(piv.to_numpy(float)), cmap=cmap, norm=norm, aspect="auto")
        for r, lu in enumerate(LAND_USE_ORDER):
            for c, sea in enumerate(SEASON_ORDER):
                val = piv.loc[lu, sea]
                if np.isnan(val):
                    ax.add_patch(Rectangle((c - 0.5, r - 0.5), 1, 1, color="#ececec", zorder=2))
                    ax.text(c, r, "–", ha="center", va="center", color="gray", fontsize=8, zorder=3)
                    continue
                if bool(thin.loc[lu, sea]):
                    ax.add_patch(Rectangle((c - 0.5, r - 0.5), 1, 1, facecolor="white", alpha=0.65,
                                           hatch="///", edgecolor="gray", lw=0, zorder=2))
                    ax.text(c, r, f"{val:.2f}†", ha="center", va="center", fontsize=7.5, color="dimgray", zorder=3)
                else:
                    ax.text(c, r, f"{val:.2f}", ha="center", va="center", fontsize=8,
                            color="white" if abs(val) > 0.55 * vmax else "black", zorder=3)
        ax.set_xticks(range(4))
        ax.set_xticklabels(SEASON_ORDER)
        ax.set_yticks(range(len(LAND_USE_ORDER)))
        ax.set_yticklabels(LAND_USE_ORDER)
        ax.set_title(loc, fontweight="bold")
        ax.set_xticks(np.arange(-0.5, 4, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, len(LAND_USE_ORDER), 1), minor=True)
        ax.grid(which="minor", color="white", lw=1.2)
        ax.tick_params(which="minor", length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=axes, fraction=0.025, pad=0.02)
    cbar.set_label("IFS (z-scale)")
    fig.legend(handles=[Patch(facecolor="white", hatch="///", edgecolor="gray", label="Thin cell (†)"),
                        Patch(facecolor="#ececec", label="No stations")],
               loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.45, -0.06))
    save_figure(fig, "fig_ifs_heatmap")
    plt.close(fig)


def main():
    ratio, eligible, qc = load_ratio_days()
    print(qc)
    pd.Series(qc).to_csv(TAB_DIR / "fdi_qc.csv")
    eligible.reset_index().to_csv(TAB_DIR / "fdi_network_stations.csv", index=False)
    pd.crosstab(eligible["Land Use"], eligible["Location Setting"], margins=True, margins_name="Total") \
        .reindex(index=LAND_USE_ORDER + ["Total"]).to_csv(TAB_DIR / "fdi_network_counts.csv")

    cells = cell_summary(ratio)
    cells["IFS"] = ifs(cells)
    cells.to_csv(TAB_DIR / "ifs_cells.csv", index=False)
    show = ["Land Use", "Location Setting", "Season", "FDI_mean", "FDI_lo", "FDI_hi", "n", "n_st", "IFS",
            "r0_mean", "r0_skew", "r0_kurt", "thin"]
    print(cells.sort_values("IFS", ascending=False)[show].round(3).to_string())
    robustness(cells, ratio)
    contrasts(ratio)
    seasonality(ratio)
    season_station_means(ratio)
    plot_boxplots(ratio, cells)
    plot_ifs_heatmap(cells)


if __name__ == "__main__":
    main()
