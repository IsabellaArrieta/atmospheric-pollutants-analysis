# EJEMPLOS DE CÓDIGO - IMPLEMENTACIÓN DE MÉTRICAS

## Parte 1: Instalación y Imports

```python
# En la primera celda del nuevo notebook

%pip install antropy
%pip install scikit-learn  # Para análisis adicional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import nolds
import antropy
from scipy.stats import skew, kurtosis
from sklearn.preprocessing import StandardScaler
from joblib import Parallel, delayed
import warnings
warnings.filterwarnings('ignore')

# Configuración visual
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)
```

---

## Parte 2: Funciones de Cálculo de Métricas

### 2.1 Complejidad Estadística

```python
def calculate_statistical_complexity(time_series):
    """
    Calcula la Complejidad Estadística usando Fisher Information.
    
    Args:
        time_series: array-like, serie de tiempo
        
    Returns:
        float: Valor de complejidad (0-1 típicamente)
    """
    try:
        # Limpiar NaNs y ceros
        series = pd.Series(time_series).dropna().values
        if len(series) < 10:  # Mínimo de datos
            return np.nan
        
        # Opción 1: Fisher Information Measure (antropy)
        # Requiere serie normalizada
        series_normalized = (series - np.mean(series)) / (np.std(series) + 1e-10)
        fisher_info = antropy.fisher_information_measure(series_normalized, order=3)
        
        return fisher_info
    except Exception as e:
        print(f"Error calculating statistical complexity: {e}")
        return np.nan


def calculate_statistical_complexity_alternative(time_series):
    """
    Alternativa: Complejidad usando Permutation Entropy + Shannon Entropy.
    """
    try:
        series = pd.Series(time_series).dropna().values
        if len(series) < 10:
            return np.nan
        
        # Permutation Entropy
        perm_entropy = antropy.permutation_entropy(series, order=3)
        
        # SVD Entropy
        svd_entropy_val = antropy.svd_entropy(series, order=3)
        
        # Complejidad = media de ambas
        complexity = (perm_entropy + svd_entropy_val) / 2
        
        return complexity
    except Exception as e:
        return np.nan
```

### 2.2 Dimensión Fractal de Higuchi

```python
def calculate_higuchi_fractal_dimension(time_series):
    """
    Calcula la Dimensión Fractal de Higuchi.
    
    Args:
        time_series: array-like, serie de tiempo
        
    Returns:
        float: Dimensión fractal (1.0-2.0)
    """
    try:
        series = pd.Series(time_series).dropna().values
        if len(series) < 10:
            return np.nan
        
        # Usar antropy
        higuchi_dim = antropy.higuchi_fractal_dimension(series, max_k=10)
        
        return higuchi_dim
    except Exception as e:
        print(f"Error calculating Higuchi dimension: {e}")
        return np.nan
```

### 2.3 Dimensión de Correlación (nolds)

```python
def calculate_correlation_dimension(time_series):
    """
    Calcula la Dimensión de Correlación usando nolds.
    
    Args:
        time_series: array-like, serie de tiempo
        
    Returns:
        float: Dimensión de correlación
    """
    try:
        series = pd.Series(time_series).dropna().values
        if len(series) < 50:  # Requiere más datos
            return np.nan
        
        # Normalizar
        series_norm = (series - np.mean(series)) / (np.std(series) + 1e-10)
        
        # Calcular dimensión de correlación
        corr_dim = nolds.corr_dim(series_norm, emb_dim=2, rvals=None)
        
        return corr_dim
    except Exception as e:
        print(f"Error calculating correlation dimension: {e}")
        return np.nan
```

### 2.4 Métricas Topológicas Complementarias

```python
def calculate_topological_metrics(time_series):
    """
    Calcula métricas topológicas de una serie.
    
    Returns:
        dict con: skewness, kurtosis, entropy_permutation, entropy_svd
    """
    try:
        series = pd.Series(time_series).dropna().values
        if len(series) < 10:
            return {}
        
        metrics = {
            'skewness': float(skew(series)),
            'kurtosis': float(kurtosis(series)),
            'entropy_permutation': float(antropy.permutation_entropy(series, order=3)),
            'entropy_svd': float(antropy.svd_entropy(series, order=3)),
            'hjorth_activity': float(antropy.hjorth_activity(series)),
            'hjorth_complexity': float(antropy.hjorth_complexity(series)),
            'hjorth_mobility': float(antropy.hjorth_mobility(series)),
        }
        return metrics
    except Exception as e:
        return {}
```

### 2.5 Función Wrapper Completa

```python
def calculate_all_complexity_metrics(time_series):
    """
    Calcula TODAS las métricas en un solo call.
    
    Returns:
        dict con todos los valores
    """
    series = pd.Series(time_series).dropna().values
    
    if len(series) < 10:
        return None
    
    result = {
        'n_samples': len(series),
        'sample_entropy': nolds.sampen(series, emb_dim=2),
        'statistical_complexity': calculate_statistical_complexity(series),
        'higuchi_fractal_dimension': calculate_higuchi_fractal_dimension(series),
        'correlation_dimension': calculate_correlation_dimension(series),
    }
    
    # Agregar métricas topológicas
    topo_metrics = calculate_topological_metrics(series)
    result.update(topo_metrics)
    
    return result
```

---

## Parte 3: Aplicar Métricas a Estaciones

### 3.1 Calcular por Estación Individual

```python
def calculate_metrics_per_station(metadata_table, pollutant_tables):
    """
    Calcula métricas para CADA estación.
    
    Args:
        metadata_table: metadata_table del preprocessor
        pollutant_tables: dict con tablas por contaminante
        
    Returns:
        DataFrame con columnas: Station_ID, Land_Use, Metrics...
    """
    
    results = []
    stations = metadata_table['Station ID'].unique()
    
    print(f"Calculando métricas para {len(stations)} estaciones...")
    
    for station_id in stations:
        # Obtener datos para esta estación
        station_data = metadata_table[metadata_table['Station ID'] == station_id]
        
        # Info de la estación
        land_use = station_data['Land Use'].iloc[0]
        location_setting = station_data['Location Setting'].iloc[0]
        lat = station_data['Latitude'].iloc[0]
        lon = station_data['Longitude'].iloc[0]
        
        # Concatenar todas las mediciones de contaminantes para esta estación
        all_measurements = station_data['Arithmetic Mean'].values
        
        # Calcular métricas
        metrics = calculate_all_complexity_metrics(all_measurements)
        
        if metrics:
            metrics['Station_ID'] = station_id
            metrics['Land_Use'] = land_use
            metrics['Location_Setting'] = location_setting
            metrics['Latitude'] = lat
            metrics['Longitude'] = lon
            results.append(metrics)
    
    return pd.DataFrame(results)
```

### 3.2 Calcular por Estación y Año

```python
def calculate_metrics_per_station_year(metadata_table):
    """
    Calcula métricas para CADA estación en CADA año.
    """
    
    results = []
    stations = metadata_table['Station ID'].unique()
    years = metadata_table['Year'].unique()
    
    for station_id in stations:
        for year in years:
            # Filtrar datos
            mask = (metadata_table['Station ID'] == station_id) & \
                   (metadata_table['Year'] == year)
            data = metadata_table[mask]
            
            if len(data) < 10:  # Mínimo de observaciones
                continue
            
            # Info
            land_use = data['Land Use'].iloc[0]
            location_setting = data['Location Setting'].iloc[0]
            
            # Mediciones
            measurements = data['Arithmetic Mean'].values
            
            # Métricas
            metrics = calculate_all_complexity_metrics(measurements)
            
            if metrics:
                metrics['Station_ID'] = station_id
                metrics['Year'] = year
                metrics['Land_Use'] = land_use
                metrics['Location_Setting'] = location_setting
                results.append(metrics)
    
    return pd.DataFrame(results)
```

### 3.3 Calcular por Estación y Climatología

```python
def calculate_metrics_per_station_season(metadata_table):
    """
    Calcula métricas para CADA estación en CADA season.
    """
    
    results = []
    stations = metadata_table['Station ID'].unique()
    seasons = metadata_table['season'].unique()
    
    for station_id in stations:
        for season in seasons:
            # Filtrar
            mask = (metadata_table['Station ID'] == station_id) & \
                   (metadata_table['season'] == season)
            data = metadata_table[mask]
            
            if len(data) < 10:
                continue
            
            # Info
            land_use = data['Land Use'].iloc[0]
            
            # Métricas
            measurements = data['Arithmetic Mean'].values
            metrics = calculate_all_complexity_metrics(measurements)
            
            if metrics:
                metrics['Station_ID'] = station_id
                metrics['Season'] = season
                metrics['Land_Use'] = land_use
                results.append(metrics)
    
    return pd.DataFrame(results)
```

---

## Parte 4: Crear Tablas Pivot

### 4.1 Pivot Zona × Contaminante (por Año)

```python
def create_pivot_zone_pollutant(metadata_table, metrics_df, metric_name, year=None):
    """
    Crea tabla pivot: Land Use (filas) × Pollutant (columnas)
    Valores: métrica especificada
    """
    
    # Mergear métricas con metadata
    if year:
        df = metadata_table[metadata_table['Year'] == year].copy()
    else:
        df = metadata_table.copy()
    
    # Agrupar por Land Use y Pollutant, calcular métrica
    pivot_data = df.groupby(['Land Use', 'Pollutant'])['Arithmetic Mean'].agg([
        ('mean', 'mean'),
        ('median', 'median'),
        ('std', 'std'),
        ('count', 'count')
    ]).reset_index()
    
    # Pivotar
    pivot = pivot_data.pivot(index='Land Use', columns='Pollutant', values='mean')
    
    return pivot

# Uso:
for year in [2023, 2024, 2025]:
    pivot = create_pivot_zone_pollutant(metadata_table, None, 'entropy', year=year)
    
    # Plotear
    plt.figure(figsize=(10, 6))
    sns.heatmap(pivot, annot=True, fmt='.2f', cmap='YlOrRd', cbar_kws={'label': 'Entropía'})
    plt.title(f'Entropía por Zona y Contaminante - Año {year}')
    plt.tight_layout()
    plt.show()
```

### 4.2 Pivot Zona × Contaminante (por Climatología)

```python
def create_pivot_zone_pollutant_season(metadata_table, season):
    """
    Idem pero filtrado por season.
    """
    
    df = metadata_table[metadata_table['season'] == season].copy()
    
    pivot_data = df.groupby(['Land Use', 'Pollutant'])['Arithmetic Mean'].agg([
        ('mean', 'mean'),
    ]).reset_index()
    
    pivot = pivot_data.pivot(index='Land Use', columns='Pollutant', values='mean')
    
    return pivot

# Uso:
fig, axes = plt.subplots(2, 2, figsize=(16, 12))
seasons = ['Winter', 'Spring', 'Summer', 'Fall']

for idx, season in enumerate(seasons):
    ax = axes[idx // 2, idx % 2]
    
    pivot = create_pivot_zone_pollutant_season(metadata_table, season)
    
    sns.heatmap(pivot, annot=True, fmt='.2f', cmap='YlOrRd', ax=ax,
                cbar_kws={'label': 'Entropía'})
    ax.set_title(f'Entropía por Zona y Contaminante - {season}')

plt.tight_layout()
plt.show()
```

---

## Parte 5: Mapas Geoespaciales

### 5.1 Scatter Plot por Año

```python
def plot_geospatial_entropy_year(station_metrics_year, year):
    """
    Plotea mapa geoespacial (Lon × Lat) coloreado por entropía.
    """
    
    df = station_metrics_year[station_metrics_year['Year'] == year]
    
    plt.figure(figsize=(12, 8))
    scatter = plt.scatter(
        df['Longitude'], 
        df['Latitude'],
        c=df['sample_entropy'],
        s=200,
        cmap='viridis',
        alpha=0.6,
        edgecolors='black',
        linewidth=0.5
    )
    
    # Añadir etiquetas
    for idx, row in df.iterrows():
        plt.annotate(row['Station_ID'], 
                    (row['Longitude'], row['Latitude']),
                    fontsize=8, alpha=0.7)
    
    plt.colorbar(scatter, label='Sample Entropy')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title(f'Distribución Espacial de Entropía - Año {year}')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

# Uso:
for year in [2023, 2024, 2025]:
    plot_geospatial_entropy_year(station_metrics_year, year)
```

### 5.2 Scatter Plot por Climatología

```python
def plot_geospatial_entropy_season(station_metrics_season, season):
    """
    Mapa geoespacial filtrado por season.
    """
    
    df = station_metrics_season[station_metrics_season['Season'] == season]
    
    plt.figure(figsize=(12, 8))
    scatter = plt.scatter(
        df['Longitude'],
        df['Latitude'],
        c=df['sample_entropy'],
        s=200,
        cmap='coolwarm',
        alpha=0.6,
        edgecolors='black'
    )
    
    plt.colorbar(scatter, label='Sample Entropy')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    plt.title(f'Distribución Espacial de Entropía - {season}')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()

# Uso:
for season in ['Winter', 'Spring', 'Summer', 'Fall']:
    plot_geospatial_entropy_season(station_metrics_season, season)
```

---

## Parte 6: Análisis Comparativo Topológico

```python
def compare_topology_by_landuse(station_metrics_year):
    """
    Compara métricas topológicas entre zonas.
    """
    
    metrics = ['sample_entropy', 'statistical_complexity', 'higuchi_fractal_dimension']
    
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    
    for idx, metric in enumerate(metrics):
        ax = axes[idx]
        
        # Box plot por Land Use
        sns.boxplot(data=station_metrics_year, x='Land_Use', y=metric, ax=ax,
                   palette='Set2')
        ax.set_title(f'{metric} por Zona')
        ax.set_ylabel(metric)
        ax.set_xlabel('Land Use')
    
    plt.tight_layout()
    plt.show()
    
    # Estadísticas descriptivas
    print("\n=== Estadísticas por Zona ===")
    for land_use in ['Rural', 'Urban', 'Military']:
        subset = station_metrics_year[station_metrics_year['Land_Use'] == land_use]
        print(f"\n{land_use}:")
        for metric in metrics:
            print(f"  {metric}: {subset[metric].mean():.4f} ± {subset[metric].std():.4f}")

# Uso:
compare_topology_by_landuse(station_metrics_year)
```

---

## Parte 7: Flujo Completo en Notebook

```python
# === FLUJO COMPLETO ===

# 1. Cargar datos (reutilizando de entropy.ipynb)
# Asume que tienes: metadata_table, station_dataset

# 2. Calcular métricas por estación
print("Calculando métricas por estación...")
station_metrics = calculate_metrics_per_station(metadata_table, final_tables)
print(f"✓ {len(station_metrics)} estaciones procesadas")

# 3. Calcular métricas por estación y año
print("\nCalculando métricas por estación y año...")
station_metrics_year = calculate_metrics_per_station_year(metadata_table)
print(f"✓ {len(station_metrics_year)} registros estación-año")

# 4. Calcular métricas por estación y season
print("\nCalculando métricas por estación y climatología...")
station_metrics_season = calculate_metrics_per_station_season(metadata_table)
print(f"✓ {len(station_metrics_season)} registros estación-season")

# 5. Crear y visualizar pivots por año
print("\nCreando tablas pivot por año...")
for year in [2023, 2024, 2025]:
    pivot = create_pivot_zone_pollutant(metadata_table, None, 'entropy', year=year)
    # Plotear...

# 6. Crear y visualizar pivots por climatología
print("\nCreando tablas pivot por climatología...")
for season in ['Winter', 'Spring', 'Summer', 'Fall']:
    pivot = create_pivot_zone_pollutant_season(metadata_table, season)
    # Plotear...

# 7. Mapas geoespaciales
print("\nCreando mapas geoespaciales...")
for year in [2023, 2024, 2025]:
    plot_geospatial_entropy_year(station_metrics_year, year)

# 8. Análisis topológico
print("\nAnálisis topológico...")
compare_topology_by_landuse(station_metrics_year)

print("\n✓✓✓ Análisis completo!")
```

---

## Notas Importantes

1. **Performance:** El cálculo de todas las métricas puede tomar ~5-10 minutos
   - Considerar usar `Parallel(n_jobs=-1)` de joblib para paralelizar

2. **Validación:** Verificar rangos esperados:
   - Sample Entropy: 0.0 - 2.0
   - Statistical Complexity: 0.0 - 1.0
   - Higuchi Dimension: 1.0 - 2.0
   - Correlation Dimension: 1.0 - 3.0

3. **NaNs:** Algunas series cortas pueden no tener métrica calculable
   - Usar `.dropna()` después de cálculos

4. **Normalización:** Algunas métricas requieren series normalizadas
   - Esto ya está hecho en las funciones

