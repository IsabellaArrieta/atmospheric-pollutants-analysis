# Atmospheric Pollutants Analysis (EPA 2023–2025)

Analysis of PM2.5 / PM10 / PM Coarse by Land Use × season: cleaning → entropy → 8 complexity metrics → maps/heatmaps.

## Structure

- `preprocessor_v2_copy.ipynb` — EPA download, Rural/Urban/Military classification, ≥75% filter, seasons, EDA.
- `entropy.ipynb` — Sample Entropy per station (`nolds`, `emb_dim=2`).
- `complexity-analysis.ipynb` — 8 metrics per station × pollutant × year × season, pivots, heatmaps, geo-scatters, topological boxplots.
- `utils_complexity.py` — reusable metric functions (`calculate_metrics`, SC Shannon MPR, Higuchi, corr_dim).
- `comparative_topology.py` — Spearman correlation matrix, Kruskal-Wallis by Land Use, Ward dendrogram + clusters (run: `python comparative_topology.py`).
- `table_csv/` — local data (not in repo, see Data).
- `outputs_complexity/` — generated pivots, heatmaps, geo maps, topology plots.
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
