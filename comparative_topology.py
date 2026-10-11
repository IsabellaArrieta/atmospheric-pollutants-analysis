"""Análisis comparativo topológico (paso 7 del plan).

Genera, a partir de outputs_complexity/metrics_by_station_pollutant_period.csv:
- Matriz de correlaciones entre métricas (CSV + heatmap).
- Kruskal-Wallis por Land Use para cada métrica (CSV).
- Dendrograma jerárquico (Ward) por estación + asignación de clusters (PNG + CSV).
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, fcluster, linkage
from scipy.stats import kruskal
from sklearn.preprocessing import StandardScaler

OUTPUT_DIR = Path("outputs_complexity")
METRICS = [
    "sample_entropy",
    "permutation_entropy",
    "svd_entropy",
    "statistical_complexity",
    "skewness",
    "kurtosis",
    "higuchi_fractal_dimension",
    "correlation_dimension",
]

sns.set_theme(style="whitegrid")


def main():
    df = pd.read_csv(OUTPUT_DIR / "metrics_by_station_pollutant_period.csv")
    df = df.dropna(subset=METRICS, how="all")

    # 1. Correlaciones entre métricas.
    corr = df[METRICS].corr(method="spearman")
    corr.to_csv(OUTPUT_DIR / "comparative_correlation_matrix.csv")
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(corr, annot=True, fmt=".2f", cmap="coolwarm", vmin=-1, vmax=1, ax=ax)
    ax.set_title("Correlación Spearman entre métricas de complejidad")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "comparative_correlation_heatmap.png", dpi=150)
    plt.close(fig)

    # 2. Kruskal-Wallis por Land Use.
    rows = []
    for metric in METRICS:
        groups = [
            g[metric].dropna().to_numpy()
            for _, g in df.groupby("Land Use")
            if len(g[metric].dropna()) >= 5
        ]
        if len(groups) >= 2:
            stat, pval = kruskal(*groups)
        else:
            stat, pval = float("nan"), float("nan")
        rows.append({"metric": metric, "kruskal_h": float(stat), "p_value": float(pval),
                     "n_groups": len(groups), "n_total": int(df[metric].notna().sum())})
    pd.DataFrame(rows).to_csv(OUTPUT_DIR / "comparative_kruskal_wallis.csv", index=False)

    # 3. Clustering jerárquico por estación (mediana de cada métrica).
    station = df.groupby("Station ID")[METRICS].median(numeric_only=True).dropna()
    land_use = df.groupby("Station ID")["Land Use"].first()
    scaled = StandardScaler().fit_transform(station[METRICS])
    linked = linkage(scaled, method="ward")
    fig, ax = plt.subplots(figsize=(14, 6))
    dendrogram(linked, labels=station.index.tolist(), no_labels=True, ax=ax)
    ax.set_title("Dendrograma Ward de estaciones por complejidad (mediana de 8 métricas)")
    ax.set_xlabel("Estaciones")
    ax.set_ylabel("Distancia Ward")
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / "comparative_dendrogram_by_station.png", dpi=150)
    plt.close(fig)

    clusters = fcluster(linked, t=3, criterion="maxclust")
    pd.DataFrame({"Station ID": station.index, "Land Use": land_use.loc[station.index].to_numpy(),
                  "cluster": clusters}).to_csv(
        OUTPUT_DIR / "comparative_clusters_by_station.csv", index=False)
    print(f"Estaciones en clustering: {len(station)}; archivos en {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
