"""Utilidades de complejidad y topología de series.

Extraído de `complexity-analysis.ipynb` para reutilizar sin duplicar
descarga/clasificación en cada notebook.
"""

from pathlib import Path

import antropy
import nolds
import numpy as np
import pandas as pd
from scipy.stats import kurtosis, skew

OUTPUT_DIR = Path("outputs_complexity")
DATA_DIR = Path("table_csv")
DATA_PATH = DATA_DIR / "metadata_table.csv"

MIN_OBSERVATIONS = 20
MIN_CORRELATION_OBSERVATIONS = 50
POLLUTANT_ORDER = ["PM Coarse", "PM10", "PM2.5 FRM/FEM", "PM2.5 non-FRM/FEM"]
LAND_USE_ORDER = ["Rural", "Urban", "Military"]
SEASON_ORDER = ["Winter", "Spring", "Summer", "Fall"]


def clean_series(values):
    series = pd.to_numeric(pd.Series(values), errors="coerce").dropna().to_numpy(dtype=float)
    if len(series) < MIN_OBSERVATIONS or not np.isfinite(series).all() or np.std(series) == 0:
        return None
    return (series - np.mean(series)) / np.std(series)


def safe_metric(function, series, **kwargs):
    try:
        value = function(series, **kwargs)
        return float(value) if np.isfinite(value) else np.nan
    except (ValueError, FloatingPointError, RuntimeError, ZeroDivisionError):
        return np.nan


def statistical_complexity_mpr(series, order=3, delay=1):
    # Complejidad estadistica de Shannon (Martin-Platino-Rosso, Batista et al. 2014): C = Q_J * H.
    series = np.asarray(series, dtype=float)
    n_vectors = len(series) - (order - 1) * delay
    if n_vectors < 2:
        return np.nan
    embedded = np.array([series[i:i + order * delay:delay] for i in range(n_vectors)])
    ranks = np.argsort(np.argsort(embedded, axis=1, kind="stable"), axis=1, kind="stable")
    _, counts = np.unique(ranks, axis=0, return_counts=True)
    n_states = 1.0
    for i in range(2, order + 1):
        n_states *= i
    full = np.zeros(int(n_states))
    full[:len(counts)] = counts / counts.sum()

    def shannon(q):
        q = q[q > 0]
        return float(-np.sum(q * np.log(q)))

    s_p = shannon(full)
    h = s_p / np.log(n_states)
    s_mix = shannon(0.5 * (full + 1.0 / n_states))
    jsd = s_mix - 0.5 * s_p - 0.5 * np.log(n_states)
    q0 = -2.0 / (((n_states + 1.0) / n_states) * np.log(n_states + 1.0) - 2.0 * np.log(2.0 * n_states) + np.log(n_states))
    return float(max(0.0, q0 * jsd * h))


def calculate_metrics(values):
    series = clean_series(values)
    if series is None:
        return None
    try:
        hjorth_mobility, hjorth_complexity = antropy.hjorth_params(series)
    except (ValueError, FloatingPointError, RuntimeError, ZeroDivisionError):
        hjorth_mobility = hjorth_complexity = np.nan
    result = {
        "n_observations": len(series),
        "sample_entropy": safe_metric(nolds.sampen, series, emb_dim=2),
        "permutation_entropy": safe_metric(antropy.perm_entropy, series, order=3, normalize=True),
        "svd_entropy": safe_metric(antropy.svd_entropy, series, order=3, normalize=True),
        "statistical_complexity": safe_metric(statistical_complexity_mpr, series, order=3, delay=1),
        "skewness": float(skew(series)),
        "kurtosis": float(kurtosis(series)),
        "hjorth_activity": float(np.var(series)),
        "hjorth_mobility": float(hjorth_mobility),
        "hjorth_complexity": float(hjorth_complexity),
        "higuchi_fractal_dimension": safe_metric(antropy.higuchi_fd, series, kmax=min(10, len(series) // 2)),
        "correlation_dimension": np.nan,
    }
    if len(series) >= MIN_CORRELATION_OBSERVATIONS:
        result["correlation_dimension"] = safe_metric(nolds.corr_dim, series, emb_dim=2)
    return result
