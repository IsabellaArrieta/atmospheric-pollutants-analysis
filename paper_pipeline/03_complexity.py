"""Stage 3: temporal complexity of daily PM series (station x pollutant x season-year windows).

Re-uses the metric implementations of utils_complexity.py on corrected series:
one instrument record (POC) per site, contiguous climatological windows (December is
assigned to the following winter), short gaps (<=3 days) linearly interpolated, and
ties broken at random before ordinal-pattern metrics.
"""

import itertools
import sys
import zlib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import fcluster, linkage
from sklearn.metrics import adjusted_rand_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import (DATA_DIR, LAND_USE_ORDER, LOCATION_ORDER, MIN_VALID_DAYS, POLLUTANT_ORDER,
                    SEASON_OF_MONTH, SEASON_ORDER, SEED, TAB_DIR)

WINDOW_COVERAGE = 0.75
MAX_GAP = 3
N_SURROGATES = 199
METRICS = ["permutation_entropy", "statistical_complexity", "sample_entropy", "svd_entropy",
           "higuchi_fractal_dimension", "correlation_dimension", "skewness", "kurtosis"]
DYNAMIC = METRICS[:6]
CLUSTER_FEATURES = ["permutation_entropy", "sample_entropy", "svd_entropy", "higuchi_fractal_dimension"]


def season_year(dates):
    return dates.dt.year + (dates.dt.month == 12).astype(int)


def window_calendar():
    days = pd.DataFrame({"date": pd.date_range("2023-01-01", "2025-12-31", freq="D")})
    days["Season"] = days["date"].dt.month.map(SEASON_OF_MONTH)
    days["season_year"] = season_year(days["date"])
    days = days[days["season_year"] <= 2025]
    return days.groupby(["season_year", "Season"])["date"].agg(["min", "max", "count"])


def ordinal_H_C(series, order=3):
    """Normalised permutation entropy and MPR statistical complexity for a batch of series (rows)."""
    series = np.atleast_2d(series)
    n = series.shape[1] - order + 1
    emb = np.stack([series[:, i:i + n] for i in range(order)], axis=-1)
    codes = np.argsort(emb, axis=-1, kind="stable")
    keys = (codes * np.array([order ** k for k in range(order)])).sum(axis=-1)
    uniq = np.array(sorted({sum(p[k] * order ** k for k in range(order))
                            for p in itertools.permutations(range(order))}))
    counts = (keys[..., None] == uniq).sum(axis=1).astype(float)
    p = counts / counts.sum(axis=1, keepdims=True)
    m = len(uniq)

    def shannon(q):
        with np.errstate(divide="ignore", invalid="ignore"):
            return -np.nansum(np.where(q > 0, q * np.log(q), 0.0), axis=1)

    s_p = shannon(p)
    h = s_p / np.log(m)
    s_mix = shannon(0.5 * (p + 1.0 / m))
    jsd = s_mix - 0.5 * s_p - 0.5 * np.log(m)
    q0 = -2.0 / (((m + 1.0) / m) * np.log(m + 1.0) - 2.0 * np.log(2.0 * m) + np.log(m))
    return h, q0 * jsd * h


def window_metrics(task):
    import antropy
    import nolds
    from utils_complexity import safe_metric, statistical_complexity_mpr

    key, values = task
    x = np.asarray(values, float)
    rng = np.random.default_rng([SEED, zlib.crc32(repr(key).encode())])
    x = (x - x.mean()) / x.std()
    # Random tie-breaking (Zunino et al., 2017): negligible jitter relative to data resolution.
    x = x + rng.uniform(-1e-6, 1e-6, size=len(x))
    # nolds.corr_dim fits the correlation sum with scikit-learn's RANSAC, which draws from the
    # global NumPy generator; seeding it per window makes D2 reproducible.
    np.random.seed(zlib.crc32(repr(key).encode()))
    out = {
        "n_obs": len(x),
        "permutation_entropy": safe_metric(antropy.perm_entropy, x, order=3, delay=1, normalize=True),
        "statistical_complexity": safe_metric(statistical_complexity_mpr, x, order=3, delay=1),
        "sample_entropy": safe_metric(nolds.sampen, x, emb_dim=2),
        "svd_entropy": safe_metric(antropy.svd_entropy, x, order=3, delay=1, normalize=True),
        "higuchi_fractal_dimension": safe_metric(antropy.higuchi_fd, x, kmax=10),
        "correlation_dimension": safe_metric(nolds.corr_dim, x, emb_dim=2),
        "skewness": float(stats.skew(values)),
        "kurtosis": float(stats.kurtosis(values)),
    }
    surrogate = np.array([rng.permutation(x) for _ in range(N_SURROGATES)])
    h_s, c_s = ordinal_H_C(surrogate)
    h_obs, _ = ordinal_H_C(x)
    out["pe_surrogate_p"] = (1 + np.sum(h_s <= h_obs[0])) / (N_SURROGATES + 1)
    out["pe_surrogate_mean"] = float(h_s.mean())
    out["c_surrogate_mean"] = float(c_s.mean())
    return key, out


def build_windows():
    daily = pd.read_parquet(DATA_DIR / "daily_long.parquet")
    meta = pd.read_parquet(DATA_DIR / "sites_meta.parquet").set_index("site_id")
    counts = daily.groupby(["site_id", "parameter"])["date"].nunique()
    keep = counts[counts >= MIN_VALID_DAYS].reset_index()[["site_id", "parameter"]]
    keep = keep.join(meta, on="site_id")
    keep = keep[keep["is_us"] & keep["Land Use"].isin(LAND_USE_ORDER) & keep["Location Setting"].isin(LOCATION_ORDER)]
    daily = daily.merge(keep[["site_id", "parameter"]], on=["site_id", "parameter"])
    calendar = window_calendar()
    tasks, info = [], []
    for (site, param), g in daily.groupby(["site_id", "parameter"]):
        s = g.set_index("date")["value"].sort_index()
        full = s.reindex(pd.date_range("2023-01-01", "2025-12-31", freq="D"))
        filled = full.interpolate(limit=MAX_GAP, limit_area="inside")
        for (sy, sea), (start, end, n_cal) in calendar.iterrows():
            observed = full.loc[start:end]
            if observed.notna().sum() < WINDOW_COVERAGE * n_cal:
                continue
            seg = filled.loc[start:end].dropna()
            if seg.std() == 0 or len(seg) < 30:
                continue
            key = (site, param, int(sy), sea)
            tasks.append((key, seg.to_numpy()))
            info.append(dict(site_id=site, parameter=param, season_year=int(sy), Season=sea,
                             calendar_days=int(n_cal), observed_days=int(observed.notna().sum()),
                             interpolated=int(seg.size - observed.notna().sum()),
                             mean_conc=float(observed.mean())))
    return tasks, pd.DataFrame(info), keep


def pattern_counts(values, segments, order, rng):
    """Ordinal-pattern counts summed over contiguous segments (no pattern spans a data gap)."""
    weights = order ** np.arange(order)
    lookup = np.full(order ** order, -1)
    for i, perm in enumerate(itertools.permutations(range(order))):
        lookup[np.dot(perm, weights)] = i
    counts = np.zeros(int(np.prod(range(1, order + 1))))
    for start, stop in segments:
        seg = values[start:stop]
        if len(seg) < order:
            continue
        seg = seg + rng.uniform(-1e-6, 1e-6, size=len(seg)) * np.std(values)
        emb = np.stack([seg[i:len(seg) - order + 1 + i] for i in range(order)], axis=1)
        codes = np.argsort(emb, axis=1, kind="stable") @ weights
        counts += np.bincount(lookup[codes], minlength=len(counts))
    return counts


def entropy_complexity(p):
    p = np.asarray(p, float)
    p = p / p.sum()
    m = len(p)
    shannon = lambda q: -np.sum(q[q > 0] * np.log(q[q > 0]))
    s_p = shannon(p)
    h = s_p / np.log(m)
    jsd = shannon(0.5 * (p + 1.0 / m)) - 0.5 * s_p - 0.5 * np.log(m)
    q0 = -2.0 / (((m + 1.0) / m) * np.log(m + 1.0) - 2.0 * np.log(2.0 * m) + np.log(m))
    return h, q0 * jsd * h


def full_period_task(task):
    key, values, segments, order = task
    rng = np.random.default_rng([SEED, zlib.crc32(repr(key).encode())])
    h, c = entropy_complexity(pattern_counts(values, segments, order, rng))
    surrogate = []
    for _ in range(99):
        h_s, c_s = entropy_complexity(pattern_counts(rng.permutation(values), segments, order, rng))
        surrogate.append((h_s, c_s))
    surrogate = np.array(surrogate)
    return key, dict(H=h, C=c, H_surr=surrogate[:, 0].mean(), C_surr=surrogate[:, 1].mean(),
                     p_ordered=(1 + np.sum(surrogate[:, 0] <= h)) / 100, n_used=int(np.isfinite(values).sum()))


def full_period_ordinal(keep, order=4):
    """Station fingerprints in the complexity-entropy plane (D=4) over the full 2023-2025 record."""
    daily = pd.read_parquet(DATA_DIR / "daily_long.parquet")
    daily = daily.merge(keep[["site_id", "parameter"]], on=["site_id", "parameter"])
    tasks = []
    for (site, param), g in daily.groupby(["site_id", "parameter"]):
        s = g.set_index("date")["value"].sort_index().reindex(pd.date_range("2023-01-01", "2025-12-31", freq="D"))
        s = s.interpolate(limit=MAX_GAP, limit_area="inside")
        valid = s.notna().to_numpy()
        edges = np.flatnonzero(np.diff(np.r_[0, valid.astype(int), 0]))
        values = s.to_numpy()[valid]
        lengths = edges[1::2] - edges[::2]
        stops = np.cumsum(lengths)
        segments = list(zip(stops - lengths, stops))
        tasks.append(((site, param), values, segments, order))
    with ProcessPoolExecutor(max_workers=8) as pool:
        results = dict(pool.map(full_period_task, tasks, chunksize=16))
    out = pd.DataFrame([dict(site_id=k[0], parameter=k[1], **v) for k, v in results.items()])
    out = out.merge(keep[["site_id", "parameter", "Land Use", "Location Setting", "Latitude", "Longitude"]],
                    on=["site_id", "parameter"])
    out.to_csv(TAB_DIR / f"cecp_full_period_D{order}.csv", index=False)
    print(out.groupby("parameter")[["H", "C", "H_surr", "C_surr"]].median().round(4))
    print(out.groupby("parameter")["p_ordered"].apply(lambda p: (p < 0.05).mean()).round(3))
    return out


def holm(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    adj = np.empty_like(p)
    running = 0.0
    for rank, idx in enumerate(order):
        running = max(running, (len(p) - rank) * p[idx])
        adj[idx] = min(running, 1.0)
    return adj


def analyses(win, keep):
    site = (win.groupby(["site_id", "parameter"])[METRICS + ["mean_conc"]].median().reset_index()
            .merge(keep[["site_id", "parameter", "Land Use", "Location Setting", "Latitude", "Longitude", "State Name"]],
                   on=["site_id", "parameter"]))
    site.to_csv(TAB_DIR / "complexity_site_medians.csv", index=False)
    site_season = (win.groupby(["site_id", "parameter", "Season"])[METRICS].median().reset_index()
                   .merge(keep[["site_id", "parameter", "Land Use", "Location Setting"]], on=["site_id", "parameter"]))
    site_season.to_csv(TAB_DIR / "complexity_site_season_medians.csv", index=False)

    # Network summary by pollutant.
    desc = win.groupby("parameter")[METRICS].agg(["median", lambda s: s.quantile(0.25), lambda s: s.quantile(0.75)])
    desc.to_csv(TAB_DIR / "complexity_descriptives_windows.csv")
    sites_per = site.groupby("parameter")["site_id"].nunique()
    windows_per = win.groupby("parameter").size()
    sig = win.groupby("parameter")["pe_surrogate_p"].apply(lambda p: (p < 0.05).mean())
    sig_season = win.groupby(["parameter", "Season"])["pe_surrogate_p"].apply(lambda p: (p < 0.05).mean()).unstack()
    pd.concat([sites_per.rename("sites"), windows_per.rename("windows"), sig.rename("share_ordered_p05")], axis=1) \
        .join(sig_season).to_csv(TAB_DIR / "complexity_network_summary.csv")
    print(pd.concat([sites_per, windows_per, sig], axis=1))
    print(sig_season.round(3))

    # Paired pollutant contrasts at collocated sites (site medians).
    pairs = [("PM2.5 FRM/FEM", "PM10"), ("PM10", "PMc"), ("PM2.5 FRM/FEM", "PMc"),
             ("PM2.5 FRM/FEM", "PM2.5 non-FRM/FEM")]
    rows = []
    for a, b in pairs:
        wa = site[site["parameter"] == a].set_index("site_id")
        wb = site[site["parameter"] == b].set_index("site_id")
        common = wa.index.intersection(wb.index)
        for metric in METRICS:
            x, y = wa.loc[common, metric], wb.loc[common, metric]
            ok = x.notna() & y.notna()
            x, y = x[ok], y[ok]
            d = x - y
            res = stats.wilcoxon(x, y)
            ranks = stats.rankdata(np.abs(d))
            rbc = (ranks[d > 0].sum() - ranks[d < 0].sum()) / ranks.sum()
            rows.append(dict(a=a, b=b, metric=metric, sites=int(ok.sum()), median_a=x.median(),
                             median_b=y.median(), median_diff=d.median(), p=res.pvalue, rank_biserial=rbc))
    paired = pd.DataFrame(rows)
    paired["p_holm"] = holm(paired["p"].to_numpy())
    paired.to_csv(TAB_DIR / "complexity_paired_pollutants.csv", index=False)
    print(paired.round(4).to_string())

    # Seasonality within sites (Friedman across the four seasons, Kendall's W).
    rows = []
    for param in POLLUTANT_ORDER:
        sub = site_season[site_season["parameter"] == param]
        for metric in METRICS:
            wide = sub.pivot(index="site_id", columns="Season", values=metric).reindex(columns=SEASON_ORDER).dropna()
            if len(wide) < 5:
                continue
            chi2, p = stats.friedmanchisquare(*[wide[s] for s in SEASON_ORDER])
            rows.append(dict(parameter=param, metric=metric, sites=len(wide), chi2=chi2, p=p,
                             kendall_W=chi2 / (len(wide) * 3),
                             **{f"median_{s}": wide[s].median() for s in SEASON_ORDER}))
    seas = pd.DataFrame(rows)
    seas["p_holm"] = holm(seas["p"].to_numpy())
    seas.to_csv(TAB_DIR / "complexity_seasonality_friedman.csv", index=False)
    print(seas.round(4).to_string())

    # Land-use and location-setting contrasts on site medians (Kruskal-Wallis, epsilon squared).
    rows = []
    for param in POLLUTANT_ORDER:
        sub = site[site["parameter"] == param]
        for factor, order in [("Land Use", LAND_USE_ORDER), ("Location Setting", LOCATION_ORDER)]:
            n_g = sub.groupby(factor)["site_id"].nunique()
            groups = [g for g in order if n_g.get(g, 0) >= 5]
            if len(groups) < 2:
                continue
            for metric in METRICS:
                data = [sub.loc[sub[factor] == g, metric].dropna().to_numpy() for g in groups]
                h, p = stats.kruskal(*data)
                n = sum(len(d) for d in data)
                rows.append(dict(parameter=param, factor=factor, metric=metric, groups=len(groups),
                                 group_list="; ".join(f"{g} ({n_g[g]})" for g in groups), sites=n, H=h, p=p,
                                 eps2=h / (n - 1), **{f"median_{g}": np.median(d) for g, d in zip(groups, data)}))
    kw = pd.DataFrame(rows)
    kw["p_holm"] = holm(kw["p"].to_numpy())
    kw.to_csv(TAB_DIR / "complexity_kruskal_landuse_location.csv", index=False)
    print(kw[["parameter", "factor", "metric", "groups", "sites", "H", "p_holm", "eps2"]].round(4).to_string())

    # Spearman correlation among metrics (window level, all pollutants pooled and by pollutant).
    win[METRICS].corr(method="spearman").to_csv(TAB_DIR / "complexity_spearman_all.csv")
    for param in POLLUTANT_ORDER:
        win.loc[win["parameter"] == param, METRICS].corr(method="spearman") \
            .to_csv(TAB_DIR / f"complexity_spearman_{param.replace('/', '-').replace(' ', '_')}.csv")

    # Year-to-year persistence of site rankings.
    rows = []
    yearly = win.groupby(["site_id", "parameter", "season_year"])[METRICS].median().reset_index()
    for param in POLLUTANT_ORDER:
        sub = yearly[yearly["parameter"] == param]
        for metric in ["permutation_entropy", "statistical_complexity", "sample_entropy", "higuchi_fractal_dimension"]:
            wide = sub.pivot(index="site_id", columns="season_year", values=metric)
            for y1, y2 in [(2023, 2024), (2024, 2025), (2023, 2025)]:
                ok = wide[[y1, y2]].dropna()
                rows.append(dict(parameter=param, metric=metric, years=f"{y1}-{y2}", sites=len(ok),
                                 rho=stats.spearmanr(ok[y1], ok[y2])[0]))
    pd.DataFrame(rows).to_csv(TAB_DIR / "complexity_year_persistence.csv", index=False)
    print(pd.DataFrame(rows).groupby(["parameter", "metric"])["rho"].agg(["min", "max"]).round(3))

    # Ward clustering of PM2.5 FRM/FEM dynamical profiles (site medians, standardized). C_JS is
    # omitted because it is a deterministic function of H for D=3, and D2 because its RANSAC-based
    # estimate is the least stable metric; the four retained metrics are non-redundant.
    pm = site[site["parameter"] == "PM2.5 FRM/FEM"].dropna(subset=CLUSTER_FEATURES).copy()
    z = ((pm[CLUSTER_FEATURES] - pm[CLUSTER_FEATURES].mean()) / pm[CLUSTER_FEATURES].std()).to_numpy()
    link = linkage(z, method="ward")
    np.save(TAB_DIR / "ward_linkage_pm25.npy", link)
    raw = fcluster(link, t=3, criterion="maxclust")
    # Label clusters from west to east by their median longitude.
    order = pd.Series(pm["Longitude"].to_numpy()).groupby(raw).median().sort_values().index
    relabel = {old: new for new, old in enumerate(order, start=1)}
    pm["cluster"] = [relabel[c] for c in raw]
    pm.to_csv(TAB_DIR / "complexity_clusters_pm25.csv", index=False)
    rng = np.random.default_rng(SEED)
    boot_ari = []
    for _ in range(200):
        idx = np.sort(rng.choice(len(pm), size=int(0.8 * len(pm)), replace=False))
        sub = fcluster(linkage(z[idx], method="ward"), t=3, criterion="maxclust")
        boot_ari.append(adjusted_rand_score(pm["cluster"].to_numpy()[idx], sub))
    six = pm.dropna(subset=DYNAMIC)
    z6 = ((six[DYNAMIC] - six[DYNAMIC].mean()) / six[DYNAMIC].std()).to_numpy()
    ari_six = adjusted_rand_score(six["cluster"], fcluster(linkage(z6, method="ward"), t=3, criterion="maxclust"))
    from sklearn.metrics import silhouette_score
    sil = {k: silhouette_score(z, fcluster(link, t=k, criterion="maxclust")) for k in range(2, 7)}
    pd.DataFrame([dict(bootstrap_ARI_median=np.median(boot_ari), bootstrap_ARI_p05=np.quantile(boot_ari, 0.05),
                       ARI_vs_six_metrics=ari_six, **{f"silhouette_k{k}": v for k, v in sil.items()})]) \
        .to_csv(TAB_DIR / "complexity_cluster_stability.csv", index=False)
    print("cluster stability", np.median(boot_ari), np.quantile(boot_ari, 0.05), ari_six, sil)
    prof = pm.groupby("cluster")[DYNAMIC + ["mean_conc", "Latitude", "Longitude"]].median()
    prof["sites"] = pm.groupby("cluster").size()
    prof.to_csv(TAB_DIR / "complexity_cluster_profiles_pm25.csv")
    ari_lu = adjusted_rand_score(pm["Land Use"], pm["cluster"])
    ari_ls = adjusted_rand_score(pm["Location Setting"], pm["cluster"])
    ct_lu = pd.crosstab(pm["Land Use"], pm["cluster"])
    ct_ls = pd.crosstab(pm["Location Setting"], pm["cluster"])
    chi_lu = stats.chi2_contingency(ct_lu)
    chi_ls = stats.chi2_contingency(ct_ls)
    cramer = lambda ct, chi: np.sqrt(chi[0] / (ct.to_numpy().sum() * (min(ct.shape) - 1)))
    lon_kw = stats.kruskal(*[g["Longitude"] for _, g in pm.groupby("cluster")])
    pd.DataFrame([dict(ARI_land_use=ari_lu, ARI_location=ari_ls, cramerV_land_use=cramer(ct_lu, chi_lu),
                       p_land_use=chi_lu[1], cramerV_location=cramer(ct_ls, chi_ls), p_location=chi_ls[1],
                       longitude_KW_H=lon_kw.statistic, longitude_KW_p=lon_kw.pvalue)]) \
        .to_csv(TAB_DIR / "complexity_cluster_association.csv", index=False)
    ct_lu.to_csv(TAB_DIR / "complexity_cluster_by_landuse.csv")
    ct_ls.to_csv(TAB_DIR / "complexity_cluster_by_location.csv")
    print(prof.round(3))
    print(ct_lu, ct_ls, f"ARI LU={ari_lu:.3f} LS={ari_ls:.3f}", lon_kw)


def main():
    tasks, info, keep = build_windows()
    print(f"windows: {len(tasks)}")
    keep.to_csv(TAB_DIR / "complexity_network_stations.csv", index=False)
    with ProcessPoolExecutor(max_workers=8) as pool:
        results = dict(pool.map(window_metrics, tasks, chunksize=64))
    metrics = pd.DataFrame([dict(zip(["site_id", "parameter", "season_year", "Season"], k), **v)
                            for k, v in results.items()])
    win = info.merge(metrics, on=["site_id", "parameter", "season_year", "Season"])
    win = win.merge(keep[["site_id", "parameter", "Land Use", "Location Setting"]], on=["site_id", "parameter"])
    win.to_csv(TAB_DIR / "complexity_windows.csv", index=False)
    analyses(win, keep)
    full_period_ordinal(keep, order=4)


if __name__ == "__main__":
    main()
