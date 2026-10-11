# ANÁLISIS DE NOTEBOOKS - Analysis Pollutants Project

## 1. ANÁLISIS DE `preprocessor_v2_copy.ipynb`

### Propósito General
**Preprocesamiento completo de datos EPA (2023-2025)** para caracterizar fuentes de contaminación ambiental basándose en la ubicación de estaciones meteorológicas.

### Estructura y Flujo

#### **Fase 1: Descarga de Datos**
- Descarga 12 datasets ZIP directamente desde EPA (aqs.epa.gov)
- **Datos por año (2023, 2024, 2025):**
  - PM2.5 FRM/FEM (código: 88101)
  - PM2.5 non-FRM/FEM (código: 88502)
  - PM10 (código: 81102)
  - PM Coarse (código: 86101)
- Stations dataset: Metadatos de todas las estaciones de monitoreo
- **Función personalizada** `leer_zip_aqs()` que maneja descargas y extracción de archivos ZIP en memoria

#### **Fase 2: Clasificación de Estaciones**
- Extrae datos de 5 columnas clave: Latitude, Longitude, Land Use, Location Setting
- **Clasificación de Land Use:**
  - Rural (prefijo: R)
  - Urban (prefijo: U)
  - Military (prefijo: M)
- **Asignación de Station ID:** `{prefijo}{número ordenado}` (ej: R01, U02, M03)
- Genera tabla `station_dataset` con metadatos clasificados por tipo de zona

#### **Fase 3: Filtrado por Observaciones Válidas**
- **Umbral:** 75% de días observados en el período 2023-2025 (1096 días totales = 822 obs mínimas)
- Elimina estaciones con datos insuficientes
- Asigna **estaciones meteorológicas** (seasons) basadas en mes:
  - Winter: Dic, Ene, Feb (DJF)
  - Spring: Mar, Abr, May (MAM)
  - Summer: Jun, Jul, Ago (JJA)
  - Fall: Sep, Oct, Nov (SON)

#### **Fase 4: Consolidación de Datos Limpios**
- Para cada contaminante:
  - Concatena 3 años de datos
  - Filtra solo estaciones válidas
  - Crea DataFrames por contaminante: `pm25_frm_fem_table`, `pm10_table`, etc.
- **Estructura final:** Matriz (Fecha Local × Station ID) con valores de concentración

#### **Fase 5: Tablas de Metadatos**
1. **`metadata_table`:** Union de todos los datos procesados
   - Columnas: Date Local, Station ID, Land Use, Location Setting, season, Pollutant, Arithmetic Mean
   - Vinculada con `station_dataset` para obtener info de zona

2. **`station_index_table`:** Tabla de referencia por estación
   - Info de cada estación: Land Use, Location Setting, Latitud, Longitud
   - Estaciones meteorológicas disponibles (seasons)
   - Contaminantes disponibles para cada estación

#### **Fase 6: Análisis Exploratorio Básico**
- **Gráficas generadas:**
  - Distribución de Land Use (countplot)
  - Distribución de Location Setting (countplot)
  - Crosstab Land Use × Location Setting (heatmap)
  - Distribución de Arithmetic Mean por contaminante (KDE)
  - Mean por Location Setting y Land Use (barplots)
  - Entropía por Location Setting (boxplot)

#### **Fase 7: Análisis Robusto para Paper**
- Resúmenes estadísticos con funciones robustas:
  - Mediana, MAD (Median Absolute Deviation)
  - Bootstrap confidence intervals
  - Análisis por: Pollutant × Year × Season
  - Análisis por: Pollutant × Season
  - Análisis por: Pollutant × Location Setting
  - Análisis por: Pollutant × Land Use
  - Resumen de Entropía × Location Setting

### Salidas Principales
```
- station_dataset: Estaciones clasificadas (Rural, Urban, Military)
- metadata_table: Todos los registros con contexto completo
- station_index_table: Índice de estaciones con metadata
- final_tables: Matrices (Fecha × Station ID) por contaminante
- Múltiples resúmenes estadísticos robustos
```

---

## 2. ANÁLISIS DE `entropy.ipynb`

### Propósito General
**Cálculo de entropía de series de tiempo** (usando Sample Entropy de nolds) y **preparación de tabla consolidada** con métricas de complejidad para estaciones.

### Estructura y Flujo

#### **Fase 1: Reutiliza Datos del Preprocessor**
- Descarga los mismos 12 datasets ZIP de EPA
- Clasifica estaciones (Land Use, Location Setting)
- Asigna Station IDs de manera idéntica
- Filtra por umbral de observaciones (822 mínimas)

#### **Fase 2: Prepara Tablas Maestras**
- **Genera las mismas tablas que el preprocessor:**
  - `station_dataset`: Clasificación de estaciones
  - `metadata_table`: Datos consolidados con contexto
  - `station_index_table`: Índice de estaciones
  - `final_tables`: Matrices por contaminante (Fecha × Station ID)

#### **Fase 3: Cálculo de Entropía de Muestra**
- **Función `entropy_calc()`:**
  - Recibe tabla principal + múltiples DataFrames de contaminantes
  - Para cada estación y cada contaminante:
    - Extrae serie de tiempo de concentraciones
    - Calcula **Sample Entropy** (nolds.sampen, embedding dim=2)
    - Promedia entropía entre contaminantes
  - Retorna tabla con columna adicional `Entropy`

- **Métrica:** Sample Entropy
  - Mide la complejidad/regularidad de la serie
  - Valores bajos = serie altamente regular (predecible)
  - Valores altos = serie más aleatoria/compleja

#### **Fase 4: Tabla Consolidada**
- `station_index_table_modified`: 
  - Todas las columnas de `station_index_table`
  - PLUS columna `Entropy` agregada

### Salidas Principales
```
- station_index_table_modified: Índice de estaciones + Entropía
- Preparado para análisis de topología y complejidad
```

---

## 3. ANÁLISIS DE PREPARACIÓN PARA LAS INSTRUCCIONES DEL USUARIO

### ✅ ESTÁ PREPARADO PARA:

1. **Tablas Pivot (Zonas en X, Contaminantes en Y):**
   - ✅ `station_dataset` contiene Land Use (Rural, Urban, Military)
   - ✅ `metadata_table` tiene Pollutant (4 tipos)
   - ✅ Se puede usar `pivot_table()` directamente

2. **Mapas de Calor con Entropía:**
   - ✅ `station_index_table_modified` ya tiene entropía calculada
   - ✅ Se puede pivotar: Land Use × Pollutant, valores = Entropy
   - ✅ Se puede hacer por año (metadata_table tiene "Year")
   - ✅ Se puede hacer por climatología (metadata_table tiene "season")

3. **Análisis por Estaciones con Coordenadas:**
   - ✅ `station_dataset` tiene Latitude y Longitude
   - ✅ Station IDs están ordenados (R01, R02, U01, etc.)
   - ✅ Se puede hacer scatter/heatmap por coordenadas

### ✅ IMPLEMENTADO en `complexity-analysis.ipynb`:

1. **Complejidad Estadística (SC):**
   - ✅ Calculada (`statistical_complexity_mpr`): SC de Shannon Martín-Platino-Rosso (Batista et al. 2014), `C = Q_J · H` sobre patrones ordinales (orden 3)

2. **Fractalidad y Dimensión Fractal de Higuchi:**
   - ✅ Calculada (`antropy.higuchi_fd`)

3. **Dimensión de Correlación (nolds):**
   - ✅ Implementada (`nolds.corr_dim`, `emb_dim=2`, series con n ≥ 50)

4. **Librería `antropy`:**
   - ✅ Importada y en uso: permutation_entropy, svd_entropy, higuchi_fd, hjorth_params

5. **Geometría Topológica de Series:**
   - ✅ Comparativo por zona (celda de topología): 8 métricas × boxplots + `topology_summary_by_land_use.csv`

### ⚠️ MEJORAS RECOMENDADAS EN DATOS EXISTENTES:

1. **Sample Entropy Actual:**
   - Está bien, pero embedding_dim=2 es bajo
   - Considerar: dim=3 o adaptativo según longitud de serie

2. **Tablas Pivot por Climatología:**
   - ✅ Implementadas en `complexity-analysis.ipynb`: pivots y heatmaps por año y por season
   - Salidas: `pivot_{metrica}_{año|season}.csv` + `heatmap_{metrica}_{año|season}.png`

3. **Correlación con Coordenadas:**
   - ✅ Implementada: 60 mapas geoespaciales (`geo_{metrica}_{año}.png` y `geo_{metrica}_{año}_{season}.png`)

---

## 4. RESUMEN COMPARATIVO

| Aspecto | preprocessor_v2 | entropy.ipynb | complexity-analysis.ipynb |
|--------|------------------|---------------|---------------------------|
| Descarga datos | ✅ Completa | ✅ Replica |
| Clasifica estaciones | ✅ Sí (Rural/Urban/Military) | ✅ Sí (igual) |
| Filtra por observaciones | ✅ Sí (75% umbral) | ✅ Sí (igual) |
| Crea Station IDs | ✅ Sí | ✅ Reutiliza |
| Asigna seasons | ✅ Sí | ✅ Reutiliza |
| Tablas por contaminante | ✅ Sí (4 matrices) | ✅ Sí (reutiliza) |
| **Calcula entropía** | ❌ NO | ✅ Sample Entropy (nolds) |
| **Visualiza entropía** | ✅ Sí (boxplot, histograma) | ❌ NO |
| Análisis por Land Use | ✅ Sí (barplot, heatmap) | ⚠️ Tabla sin visualizar |
| Análisis por Season | ⚠️ Parcial | ⚠️ Parcial |
| **Complejidad estadística** | ❌ NO | ❌ NO | ✅ SC Shannon (`complexity-analysis.ipynb`) |
| **Dimensión fractal Higuchi** | ❌ NO | ❌ NO | ✅ Higuchi (`complexity-analysis.ipynb`) |
| **Dimensión de correlación** | ❌ NO | ❌ NO | ✅ `nolds.corr_dim` (`complexity-analysis.ipynb`) |

---

## 5. PLAN DE ACCIÓN — ESTADO (pasos 1-3 completados en `complexity-analysis.ipynb`)

### **PASO 1: Nuevas métricas** ✅
```
→ antropy instalado y en uso
→ Complejidad Estadística Shannon (MPR) implementada
→ Dimensión Fractal de Higuchi implementada
→ Dimensión de Correlación (nolds.corr_dim) implementada
→ Tabla consolidada: metrics_by_station_pollutant_period.csv
```

### **PASO 2: Tablas Pivot y Mapas de Calor** ✅
```
→ Pivots Land Use × Pollutant con Entropía (por año y por season)
→ Pivots con Complejidad Estadística y Dimensión Fractal
→ Heatmaps correspondientes en outputs_complexity/
```

### **PASO 3: Análisis Geoespacial por Estaciones** ✅
```
→ Scatter lat/lon por año y por season-año, color = métrica
→ Cubierto para las 4 métricas principales (60 mapas)
→ Pendiente opcional: agrupar estaciones cercanas por clusters
```

### **PASO 4: Análisis de Topología y Diferencias entre Series** ✅ HECHO (2026-10-11, `comparative_topology.py`)
```
→ Boxplots por zona (8 métricas, ya existían)
→ Matriz de correlación Spearman: comparative_correlation_matrix.csv + heatmap
→ Kruskal-Wallis por Land Use: comparative_kruskal_wallis.csv (las 8 métricas con p < 0.005)
→ Dendrograma Ward por estación (1021 estaciones): comparative_dendrogram_by_station.png
→ Clusters (k=3): comparative_clusters_by_station.csv
```

---

## CONCLUSIÓN

- **`preprocessor_v2_copy.ipynb`** ✅ Carga y prepara datos correctamente
- **`entropy.ipynb`** ✅ Calcula entropía (base reutilizada por el análisis de complejidad)
- **`complexity-analysis.ipynb`** ✅ Complejidad estadística (SC Shannon), fractalidad (Higuchi + correlación), pivots/heatmaps, mapas geoespaciales y análisis topológico
- **`comparative_topology.py`** ✅ HECHO (2026-10-11): correlación entre métricas, Kruskal-Wallis por Land Use, dendrograma + clusters por estación
- **`utils_complexity.py`** ✅ HECHO (2026-10-11): funciones extraídas del notebook (`clean_series`, `safe_metric`, `statistical_complexity_mpr`, `calculate_metrics`)
- **Reutilización:** `utils_complexity.py` elimina la duplicación descarga/clasificación entre notebooks — usarlo en nuevos análisis

