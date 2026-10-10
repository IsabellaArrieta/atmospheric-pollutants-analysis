# PLAN DE IMPLEMENTACIÓN - INSTRUCCIONES DEL USUARIO

## INSTRUCCIONES A IMPLEMENTAR

### 1. Tablas Pivot con Mapas de Calor - Entropía

#### 1.1 Por Año
**Estructura:** Land Use (Rural/Urban/Military) × Contaminantes (PM2.5 FRM/FEM, PM2.5 non-FRM/FEM, PM10, PM Coarse)
- **Valores:** Entropía promedio por combinación zona-contaminante
- **Mapas de calor:** Uno por cada año (2023, 2024, 2025)
- **Variables X (zona):** Rural, Urban, Military
- **Variables Y (contaminante):** PM2.5 FRM/FEM, PM2.5 non-FRM/FEM, PM10, PM Coarse
- **Interior (color/valor):** Entropía

#### 1.2 Por Climatología (Seasons)
**Estructura:** Land Use × Contaminantes, separados por season
- **Mapas de calor:** Uno por season (Winter, Spring, Summer, Fall)
- **Mismo layout que 1.1** pero filtrado por season

#### DATOS NECESARIOS:
```
metadata_table contiene:
- Land Use ✅
- Pollutant ✅
- season (Climatología) ✅
- Arithmetic Mean (valores de contaminante) ✅
- Falta: Entropía desagregada por pollutant/season/land_use
```

**Pseudocódigo:**
```python
# Agrupar por Land Use, Pollutant
grouped = metadata_table.groupby(['Land Use', 'Pollutant'])['Arithmetic Mean'].agg([estadísticas])

# Calcular entropía para cada grupo
entropy_by_group = ...

# Pivotar a matriz
pivot = entropy_by_group.pivot(index='Land Use', columns='Pollutant', values='Entropy')

# Plotear heatmap
sns.heatmap(pivot, annot=True, cmap='YlOrRd')
```

---

### 2. Mapas de Calor por Estaciones - Coordenadas

#### 2.1 Por Año
**Estructura:** Scatter plot con coordenadas geográficas
- **X axis:** Longitude
- **Y axis:** Latitude
- **Color/Tamaño:** Entropía de la estación
- **Un mapa por cada año** (2023, 2024, 2025)
- **Estación ID en etiquetas**

#### 2.2 Por Climatología
**Estructura:** Similar, pero con submapas por season
- **4 mapas por año** (uno por Winter, Spring, Summer, Fall)
- **Total:** 12 mapas geoespaciales

#### DATOS NECESARIOS:
```
station_index_table_modified contiene:
- Station ID ✅
- Latitude ✅
- Longitude ✅
- Entropy ✅

Falta: Entropía desagregada por año
Falta: Entropía desagregada por season
Falta: Entropía desagregada por pollutant
```

**Pseudocódigo:**
```python
# Para cada año
for year in [2023, 2024, 2025]:
    df_year = metadata_table[metadata_table['Year'] == year]
    
    # Calcular entropía por estación para ese año
    entropy_year = df_year.groupby('Station ID').apply(entropy_calc)
    
    # Mergear con coordenadas
    df_plot = entropy_year.merge(station_dataset, on='Station ID')
    
    # Plotear
    plt.scatter(df_plot['Longitude'], df_plot['Latitude'], 
                c=df_plot['Entropy'], s=100, cmap='viridis')
    plt.colorbar(label='Entropía')
```

---

## MÉTRICAS A CALCULAR (Punto 2-4 de las instrucciones)

### 2. Geometría Topológica de Series
**"Mirar geometría topológica para diferencias entre series"**

Marcadores a calcular:
- ✅ **Entropía** (ya existe: Sample Entropy)
- ❌ **Forma** (no existe)
  - Skewness (asimetría)
  - Kurtosis (curtosis/pico)
  - Correlación con distribución normal
  
- ❌ **Complejidad Estadística** (no existe)
  - Definición: Complejidad de Lempel-Ziv o Fisher Information
  - Librerías: antropy.fisher_information_measure(), implementación personalizada
  
- ❌ **Fractalidad - Dimensión Fractal de Higuchi** (no existe)
  - Mide auto-similaridad de la serie
  - Biblioteca: antropy o implementación (antropy.higuchi_fractal_dimension())
  - Rango típico: 1.0-2.0 (1=línea recta, 2=brownian motion)

### 3. Usar nolds para Fractalidad - Dimensión de Correlación
**"Ya lo usamos, usar también para fractalidad (correlation dimension)"**

- ✅ nolds está importado
- ❌ No se usa correlation dimension
- Función: `nolds.corr_dim(serie, emb_dim, rvals=None)`
  - Mide dimensión fractal
  - Diferente a Higuchi, basada en correlaciones espaciales
  - Útil para identificar caos en series

### 4. Complejidad Estadística - Antropy
**"Podemos usar antropy como alternativa"**

- ❌ antropy NO está instalado
- Funciones disponibles en antropy:
  - `permutation_entropy()` - Entropy de permutación
  - `svd_entropy()` - Entropy SVD
  - `fisher_information_measure()` - Complejidad Fisher
  - `higuchi_fractal_dimension()` - Dimensión fractal de Higuchi
  - `detrended_fluctuation_analysis()` - DFA
  - `hjorth_activity()`, `hjorth_complexity()`, `hjorth_mobility()` - Parámetros Hjorth

---

## RESUMEN DE IMPLEMENTACIÓN NECESARIA

### LIBRERÍA Y FUNCIONES A AGREGAR

```python
%pip install antropy

# Importar
import antropy
from scipy.stats import skew, kurtosis
import nolds

# Nuevas funciones a crear:

def calculate_topological_metrics(time_series):
    """
    Calcula métricas topológicas de una serie de tiempo
    Returns: dict con:
      - entropy_sample (nolds.sampen)
      - entropy_permutation (antropy)
      - entropy_svd (antropy)
      - skewness
      - kurtosis
      - fisher_info (antropy)
    """

def calculate_fractal_dimension(time_series):
    """
    Calcula dimensiones fractales
    Returns: dict con:
      - higuchi_dim (antropy)
      - correlation_dim (nolds)
    """

def calculate_statistical_complexity(time_series):
    """
    Calcula complejidad estadística
    Returns: float (valor de complejidad)
    """
```

### TABLAS A CREAR

1. **Tabla de Métricas por Estación:**
   ```
   Columnas: Station ID, Land Use, Pollutant, Year, Season,
             Entropy_Sample, Entropy_Permutation, Entropy_SVD,
             Skewness, Kurtosis, Fisher_Info,
             Higuchi_Dim, Correlation_Dim,
             Statistical_Complexity
   ```

2. **Tablas Pivot - Año:**
   ```
   Para cada año y cada métrica:
   - Rows: Land Use (Rural, Urban, Military)
   - Cols: Pollutant (4 tipos)
   - Values: Métrica promedio
   ```

3. **Tablas Pivot - Climatología:**
   ```
   Para cada season y cada métrica:
   - Rows: Land Use
   - Cols: Pollutant
   - Values: Métrica promedio
   ```

4. **Mapas Geoespaciales:**
   ```
   Para cada año/season:
   - Scatter plot (Lon × Lat)
   - Color: Métrica (entropía, fractalidad, etc.)
   - Etiquetas: Station ID
   ```

---

## ORDEN DE IMPLEMENTACIÓN RECOMENDADO

### **PASO 1: Crear Notebook Extendido "complexity-analysis.ipynb"**
- Importar datos base de entropy.ipynb
- Agregar librería antropy
- Crear función `calculate_topological_metrics()`
- Crear función `calculate_fractal_dimension()`
- Crear función `calculate_statistical_complexity()`

### **PASO 2: Agregar Métricas a Tabla Base**
- Calcular todas las métricas para cada estación
- Crear tabla consolidada: `station_metrics_table`
- Guardar para reutilizar

### **PASO 3: Pivots por Año**
- Para cada métrica (entropía, fractalidad, complejidad):
  - Crear pivot: Land Use × Pollutant
  - Plotear heatmap
  - Repetir para 2023, 2024, 2025

### **PASO 4: Pivots por Climatología**
- Idem Paso 3 pero para Winter, Spring, Summer, Fall
- Para cada métrica

### **PASO 5: Mapas Geoespaciales por Año**
- Para cada año (2023, 2024, 2025):
  - Crear scatter plot: Lon × Lat
  - Colorear por entropía
  - Repetir para cada métrica (min. 4 mapas por año)

### **PASO 6: Mapas Geoespaciales por Climatología**
- Para cada season (Winter, Spring, Summer, Fall):
  - Crear scatter plot por año (12 mapas totales = 3 años × 4 seasons)
  - Colorear por métrica
  - Considerar: panels/subplots

### **PASO 7: Análisis Comparativo**
- Comparar geometría topológica entre zonas
- Dendrogramas de clustering por complejidad
- Matriz de correlaciones entre métricas
- Análisis de diferencias estadísticas (ANOVA/Kruskal-Wallis)

---

## NOTAS IMPORTANTES

### Estado Actual:
- ✅ Datos limpios y preprocesados
- ✅ Sample Entropy calculada
- ✅ Tablas de referencia (Land Use, Location Setting, Seasons)
- ❌ Falta: Complejidad estadística
- ❌ Falta: Dimensión fractal completa
- ❌ Falta: Análisis topológico de diferencias

### Dependencias a Instalar:
```
antropy           # Para complejidad y fractalidad alternativa
nolds             # Ya existe, pero subutilizado
scipy.stats       # Para skewness/kurtosis (ya incluida en pandas)
matplotlib/seaborn # Para visualizaciones
```

### Consideraciones de Performance:
- Cálculo de múltiples métricas por ~800 estaciones puede ser lento
- Considerar: Parallelización con joblib
- Cachear resultados intermedios

### Validación de Resultados:
- Verificar rango de valores esperados:
  - Entropía sample: 0.0 - ~2.0
  - Dimensión fractal: 1.0 - 2.0
  - Complejidad estadística: [0, 1] típicamente
  - Skewness: -∞ a +∞ (usualmente -3 a +3)
  - Kurtosis: 0 (normal) a ∞

