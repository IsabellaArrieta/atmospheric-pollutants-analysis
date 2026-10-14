"""Stage 5: write LaTeX table bodies from the CSV outputs (outputs_paper/latex/*.tex).

Generated fragments are pasted into the manuscript so that every number in the tables
comes directly from the pipeline outputs.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA_DIR, LAND_USE_ORDER, LOCATION_ORDER, OUT_DIR, POLLUTANT_ORDER, SEASON_ORDER, TAB_DIR

LTX = OUT_DIR / "latex"
LTX.mkdir(exist_ok=True)
LOC_SHORT = {"Rural": "Rural", "Suburban": "Suburban", "Urban and center city": "Urban and center city"}
POLL_TEX = {"PM2.5 FRM/FEM": r"PM$_{2.5}$ FRM/FEM", "PM2.5 non-FRM/FEM": r"PM$_{2.5}$ non-FRM/FEM",
            "PM10": r"PM$_{10}$", "PMc": r"PM$_{c}$"}


def fmt_p(p):
    if p < 0.001:
        return r"$<$0.001"
    return f"{p:.3f}"


def write(name, text):
    (LTX / name).write_text(text, encoding="utf-8")
    print(f"--- {name}\n{text}")


def inventory():
    sites = pd.read_parquet(DATA_DIR / "sites_meta.parquet")
    sites = sites[sites["is_us"] & sites["is_open"]]
    lu = sites["Land Use"].fillna("Unknown").replace({"": "Unknown"}).value_counts()
    ls = sites["Location Setting"].fillna("Unknown").replace({"": "Unknown"}).value_counts()
    lu_rows = "\n".join(f"{k} & {v} \\\\" for k, v in lu.items())
    ls_rows = "\n".join(f"{k} & {v} \\\\" for k, v in ls.items())
    write("tab_inventory.tex", f"% land use ({lu.sum()} sites)\n{lu_rows}\n% location setting ({ls.sum()} sites)\n{ls_rows}\n")


def networks():
    fdi = pd.read_csv(TAB_DIR / "fdi_network_counts.csv", index_col=0)
    cx = pd.read_csv(TAB_DIR / "complexity_network_stations.csv")
    cx_counts = cx.pivot_table(index="Land Use", columns="parameter", values="site_id", aggfunc="nunique").fillna(0)
    rows = []
    for lu in LAND_USE_ORDER + ["Total"]:
        f = fdi.loc[lu] if lu in fdi.index else pd.Series(0, index=fdi.columns)
        if lu == "Total":
            c = cx.groupby("parameter")["site_id"].nunique()
        else:
            c = cx_counts.loc[lu] if lu in cx_counts.index else pd.Series(0, index=POLLUTANT_ORDER)
        cells = [int(f.get("Rural", 0)), int(f.get("Suburban", 0)), int(f.get("Urban and center city", 0)), int(f.get("Total", 0))]
        cells += [int(c.get(p, 0)) for p in POLLUTANT_ORDER]
        label = r"\textbf{Total}" if lu == "Total" else lu
        body = " & ".join(str(x) for x in cells)
        if lu == "Total":
            rows.append(r"\midrule")
            body = " & ".join(rf"\textbf{{{x}}}" for x in cells)
        rows.append(f"{label} & {body} \\\\")
    union = cx["site_id"].nunique()
    write("tab_networks.tex", "\n".join(rows) + f"\n% union of complexity stations: {union}\n")


def hotcold():
    cells = pd.read_csv(TAB_DIR / "ifs_cells.csv")
    ok = cells[~cells["thin"]]

    def row(r):
        return (f"{r['Land Use']} & {r['Location Setting']} & {r['Season']} & {r['FDI_mean']:.2f} & "
                f"[{r['FDI_lo']:.2f}, {r['FDI_hi']:.2f}] & {int(r['n'])} & {int(r['n_st'])} & {r['IFS']:.2f} & "
                f"{r['r0_mean']:.2f} & {r['r0_skew']:.2f} & {r['r0_kurt']:.2f} \\\\")
    top = "\n".join(row(r) for _, r in ok.sort_values("IFS", ascending=False).head(10).iterrows())
    bot = "\n".join(row(r) for _, r in ok.sort_values("IFS").head(10).iterrows())
    write("tab_hotcold.tex", f"% non-thin cells: {len(ok)} of {len(cells)}\n% TOP\n{top}\n% BOTTOM\n{bot}\n")
    blocks = []
    for loc in LOCATION_ORDER:
        sub = ok[(ok["Location Setting"] == loc) & (ok["IFS"] > 0)].sort_values("IFS", ascending=False).head(5)
        blocks.append(f"% {loc}\n" + "\n".join(row(r) for _, r in sub.iterrows()))
    write("tab_hot_by_loc.tex", "\n".join(blocks) + "\n")
    thin = cells[cells["thin"]].sort_values("IFS")
    write("thin_cells.txt", thin[["Land Use", "Location Setting", "Season", "n", "n_st", "FDI_mean", "FDI_lo",
                                  "FDI_hi", "IFS"]].round(3).to_string())


def fdi_contrasts():
    omni = pd.read_csv(TAB_DIR / "fdi_kruskal_wallis.csv")
    post = pd.read_csv(TAB_DIR / "fdi_posthoc_cliffs_delta.csv")
    units = pd.read_csv(TAB_DIR / "fdi_station_season_year_units.csv")
    rows = []
    for loc in LOCATION_ORDER:
        for sea in SEASON_ORDER:
            o = omni[(omni["Location"] == loc) & (omni["Season"] == sea) & (omni["metric"] == "FDI")]
            if o.empty:
                continue
            o = o.iloc[0]
            pp = post[(post["Location"] == loc) & (post["Season"] == sea) & (post["metric"] == "FDI")].copy()
            pp["abs"] = pp["delta"].abs()
            best = pp.sort_values("abs", ascending=False).iloc[0]
            a, b, d = best["a"], best["b"], best["delta"]
            if d < 0:
                a, b, d = b, a, -d
            sub = units[(units["Location Setting"] == loc) & (units["Season"] == sea)]
            n_st = sub.groupby("Land Use")["site_id"].nunique()
            groups = [g for g in LAND_USE_ORDER if n_st.get(g, 0) >= 3]
            gtxt = ", ".join(f"{g[:3]}.\\,({n_st[g]})" for g in groups)
            rows.append(f"{loc} & {sea} & {gtxt} & {int(o['units'])} & {o['H']:.2f} & {o['eps2']:.3f} & "
                        f"{fmt_p(o['p_holm'])} & {a} $>$ {b} & {d:.2f} & {fmt_p(best['p_holm'])} \\\\")
    write("tab_fdi_contrasts.tex", "\n".join(rows) + "\n")


def robustness():
    grid = pd.read_csv(TAB_DIR / "ifs_weight_grid.csv")
    named = pd.read_csv(TAB_DIR / "ifs_robustness_named.csv")
    years = pd.read_csv(TAB_DIR / "ifs_year_replication.csv")
    lines = [f"Weight grid (625 specifications), median [min] & {grid['tau'].median():.2f} [{grid['tau'].min():.2f}] & "
             f"{grid['top10'].median():.1f} [{grid['top10'].min():.1f}] & {grid['bottom10'].median():.1f} [{grid['bottom10'].min():.1f}] \\\\"]
    for _, r in named.iterrows():
        lines.append(f"{r['spec']} & {r['tau']:.2f} & {r['top10']:.1f} & {r['bottom10']:.1f} \\\\")
    for _, r in years.iterrows():
        lines.append(f"Season-year {int(r['year'])} only & {r['tau']:.2f} & {r['top10']:.1f} & {r['bottom10']:.1f} \\\\")
    write("tab_ifs_robustness.tex", "\n".join(lines) + "\n")


def complexity_summary():
    win = pd.read_csv(TAB_DIR / "complexity_windows.csv")
    cecp = pd.read_csv(TAB_DIR / "cecp_full_period_D4.csv")
    summ = pd.read_csv(TAB_DIR / "complexity_network_summary.csv", index_col=0)
    metrics = [("permutation_entropy", r"Permutation entropy $H$ ($D=3$)", 3),
               ("statistical_complexity", r"Statistical complexity $C_{JS}$ ($D=3$)", 3),
               ("sample_entropy", r"Sample entropy", 2),
               ("svd_entropy", r"SVD entropy", 3),
               ("higuchi_fractal_dimension", r"Higuchi fractal dimension", 2),
               ("correlation_dimension", r"Correlation dimension $D_2$", 2),
               ("skewness", r"Skewness", 2),
               ("kurtosis", r"Excess kurtosis", 2)]
    rows = ["Stations & " + " & ".join(str(int(summ.loc[p, "sites"])) for p in POLLUTANT_ORDER) + r" \\",
            "Season-year windows & " + " & ".join(str(int(summ.loc[p, "windows"])) for p in POLLUTANT_ORDER) + r" \\",
            r"Windows with $H$ below surrogates (\%) & " + " & ".join(f"{100 * summ.loc[p, 'share_ordered_p05']:.1f}" for p in POLLUTANT_ORDER) + r" \\",
            r"\midrule"]
    for m, label, dec in metrics:
        vals = []
        for p in POLLUTANT_ORDER:
            s = win.loc[win["parameter"] == p, m].dropna()
            vals.append(f"{s.median():.{dec}f} [{s.quantile(0.25):.{dec}f}, {s.quantile(0.75):.{dec}f}]")
        rows.append(f"{label} & " + " & ".join(vals) + r" \\")
    rows.append(r"\midrule")
    for m, label in [("H", r"$H$ ($D=4$, full period)"), ("C", r"$C_{JS}$ ($D=4$, full period)")]:
        vals = []
        for p in POLLUTANT_ORDER:
            s = cecp.loc[cecp["parameter"] == p, m]
            vals.append(f"{s.median():.3f} [{s.quantile(0.25):.3f}, {s.quantile(0.75):.3f}]")
        rows.append(f"{label} & " + " & ".join(vals) + r" \\")
    write("tab_complexity_summary.tex", "\n".join(rows) + "\n")


def complexity_tests():
    seas = pd.read_csv(TAB_DIR / "complexity_seasonality_friedman.csv")
    kw = pd.read_csv(TAB_DIR / "complexity_kruskal_landuse_location.csv")
    pair = pd.read_csv(TAB_DIR / "complexity_paired_pollutants.csv")
    metrics = [("permutation_entropy", r"$H$"), ("sample_entropy", "SampEn"), ("svd_entropy", "SVDEn"),
               ("higuchi_fractal_dimension", "HFD"), ("correlation_dimension", r"$D_2$"), ("kurtosis", "Kurt")]

    def star(p):
        return "" if p >= 0.05 else (r"$^{*}$" if p >= 0.01 else (r"$^{**}$" if p >= 0.001 else r"$^{***}$"))
    rows = []
    for m, label in metrics:
        cells = []
        for p in ["PM2.5 FRM/FEM", "PM10", "PMc"]:
            r = seas[(seas["parameter"] == p) & (seas["metric"] == m)].iloc[0]
            cells.append(f"{r['kendall_W']:.3f}{star(r['p_holm'])}")
        for p in ["PM2.5 FRM/FEM", "PM10"]:
            r = kw[(kw["parameter"] == p) & (kw["factor"] == "Land Use") & (kw["metric"] == m)].iloc[0]
            cells.append(f"{r['eps2']:.3f}{star(r['p_holm'])}")
        for a, b in [("PM2.5 FRM/FEM", "PM10"), ("PM10", "PMc")]:
            r = pair[(pair["a"] == a) & (pair["b"] == b) & (pair["metric"] == m)].iloc[0]
            cells.append(f"{r['rank_biserial']:+.2f}{star(r['p_holm'])}")
        rows.append(f"{label} & " + " & ".join(cells) + r" \\")
    write("tab_complexity_tests.tex", "\n".join(rows) + "\n")


def link():
    out = pd.read_csv(TAB_DIR / "link_fdi_complexity.csv")
    metrics = [("H_D4", r"$H$ ($D=4$, full period)"), ("C_D4", r"$C_{JS}$ ($D=4$, full period)"),
               ("sample_entropy", "Sample entropy"), ("svd_entropy", "SVD entropy"),
               ("higuchi_fractal_dimension", "Higuchi fractal dimension"), ("skewness", "Skewness"),
               ("kurtosis", "Excess kurtosis")]
    rows = []
    for m, label in metrics:
        cells = []
        for p in ["PM2.5 FRM/FEM", "PM10"]:
            r = out[(out["parameter"] == p) & (out["metric"] == m)].iloc[0]
            mark = r"$^{\dagger}$" if r["p_holm"] < 0.05 else ""
            cells += [f"{r['rho']:+.2f}{mark}", f"{r['partial_rho']:+.2f}"]
        rows.append(f"{label} & " + " & ".join(cells) + r" \\")
    write("tab_link.tex", "\n".join(rows) + "\n")


def main():
    inventory()
    networks()
    hotcold()
    fdi_contrasts()
    robustness()
    complexity_summary()
    complexity_tests()
    link()


if __name__ == "__main__":
    main()
