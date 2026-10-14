# Atmospheric Pollutants Analysis (EPA 2023–2025)

Analysis of PM2.5 / PM10 / PM Coarse by Land Use × season: cleaning → entropy → 8 complexity metrics → maps/heatmaps.

## Structure

- `preprocessor_v2_copy.ipynb` — EPA download, Rural/Urban/Military classification, ≥75% filter, seasons, EDA.
- `entropy.ipynb` — Sample Entropy per station (`nolds`, `emb_dim=2`).
- `complexity-analysis.ipynb` — 8 metrics per station × pollutant × year × season, pivots, heatmaps, geo-scatters, topological boxplots.
- `utils_complexity.py` — reusable metric functions (`calculate_metrics`, SC Shannon MPR, Higuchi, corr_dim).
- `comparative_topology.py` — Spearman correlation matrix, Kruskal-Wallis by Land Use, Ward dendrogram + clusters (run: `python comparative_topology.py`).
- `table_csv/` — local data (not in repo, see Data).
- `outputs_complexity/` — generated pivots, heatmaps, geo maps, topology plots. **Superseded**: these exploratory outputs were computed on series with duplicated days and merged station IDs; see `md/REVISION_PAPER_2026-10.md`. Use `outputs_paper/` instead.
- `paper_pipeline/` — reproducible pipeline behind the manuscript (FDI/IFS + temporal complexity, 2023–2025).
- `outputs_paper/` — tables (`tables/`), figures (`figures/`) and LaTeX table bodies (`latex/`) produced by `paper_pipeline/`.
- `docs/Clasiffication_Pollution_Sources_EPA/` — manuscript (`cas-sc-template.tex`, `cas-refs.bib`, `figs/`); the submitted version is kept in `original_submission/`.
- `md/` — working docs.

## Data

`table_csv/metadata_table.csv` (~250 MB) is **not** in the repo. Download it:

- Drive (view): https://drive.google.com/file/d/16sAAikAraZOXHQwMhrWZsAhQbnnsu6jY/view?usp=drive_link
- Direct download: https://drive.google.com/uc?id=16sAAikAraZOXHQwMhrWZsAhQbnnsu6jY

Place it at `table_csv/metadata_table.csv`. The link must be shared as "Anyone with the link can view".

## Reproduce

```bash
pip install -r requirements.txt
jupyter notebook
```

Run in order: `preprocessor_v2_copy.ipynb` → `entropy.ipynb` → `complexity-analysis.ipynb`.

## Reproduce the manuscript (`paper_pipeline/`)

Requires Python ≥ 3.12 (`nolds` 0.6.3 does not import on 3.11). The first stage downloads the
public EPA AirData files (~40 MB) into `table_csv/raw/`; no Drive file is needed.

```bash
pip install -r requirements.txt
python paper_pipeline/01_build_dataset.py   # download, one value per site/parameter/day, screening
python paper_pipeline/02_fdi_ifs.py         # FDI, IFS cells, robustness, station-level contrasts
python paper_pipeline/03_complexity.py      # complexity metrics per season-year window, tests, clusters, D=4 plane
python paper_pipeline/04_link_and_maps.py   # maps, complexity figures, FDI-complexity link
python paper_pipeline/05_latex_tables.py    # LaTeX table bodies used in the manuscript
```

Figures are written to `outputs_paper/figures/` (PNG) and `docs/Clasiffication_Pollution_Sources_EPA/figs/` (PDF).
All random procedures use fixed seeds, so two consecutive runs produce identical outputs.
