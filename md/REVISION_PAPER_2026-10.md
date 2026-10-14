# Revisión del manuscrito (octubre 2026)

Manuscrito: `docs/Clasiffication_Pollution_Sources_EPA/cas-sc-template.tex` (PDF compilado al lado).
Versión enviada a la revista: `docs/Clasiffication_Pollution_Sources_EPA/original_submission/`.
Pipeline que reproduce todas las tablas y figuras: `paper_pipeline/` (ver README).

## 1. Problemas encontrados en los resultados nuevos (`complexity-analysis.ipynb`, `outputs_complexity/`)

| # | Problema | Evidencia | Consecuencia |
|---|----------|-----------|--------------|
| 1 | Los archivos diarios de la EPA traen varias filas por sitio y día (POC × duración "1 HOUR"/"24-HR BLK AVG" × estándar). No se deduplicaron. | Hasta 10 filas por sitio-día en 88101; en `metrics_by_station_pollutant_period.csv` el 68 % de las series (96 % de PM2.5 FRM/FEM) tiene más de 92 observaciones por estación-año, con un máximo de 868. | Las métricas se calcularon sobre series con días repetidos e intercalados. |
| 2 | Colisión de `Station ID`: la celda de preparación integrada mapea RESIDENTIAL, INDUSTRIAL, DESERT y FOREST a "Rural" y numera por uso de suelo original, así que "Rural1" mezcla cuatro estaciones distintas. | 1104 estaciones únicas en `entropy.ipynb` frente a 1021 en `complexity-analysis.ipynb`, con el mismo número de registros. | Series de estaciones distintas mezcladas. |
| 3 | La categoría "Rural/Urban/Military" no es un uso de suelo: llama "Rural" a sitios residenciales e industriales y "Urban" solo a los comerciales; además se excluyó "Agricultural". | Diccionario `land_use_map` de la celda 4. | Las comparaciones por "zona" no son interpretables. |
| 4 | El resultado principal (PM2.5 mucho más "complejo" que PM10: PE ≈ 0,65 frente a 0,95; C ≈ 0,20 frente a 0,04) es un artefacto de los puntos 1–2. | Con series corregidas, PE ≈ 0,965 para PM2.5 y 0,969 para PM10. | No debe reportarse. |
| 5 | Las ventanas de invierno juntan enero–febrero con diciembre del mismo año (no contiguas). | `season` por mes dentro de `Year`. | Patrones ordinales espurios en el salto feb→dic. |
| 6 | Empates frecuentes (PM10 se reporta en enteros) sesgan la entropía de permutación (Zunino et al., 2017). | — | Sesgo hacia menor entropía. |
| 7 | `nolds.corr_dim` usa RANSAC sin semilla: la dimensión de correlación cambia entre ejecuciones. | Dos corridas consecutivas daban D2 distintos. | Resultados no reproducibles. |
| 8 | Con D = 3, PE y complejidad estadística son redundantes (ρ = −1,00). `hjorth_activity` vale 1 siempre (serie estandarizada). | Matriz de correlación. | No aportan información independiente. |

## 2. Problemas encontrados en el paper enviado (FDI/IFS, 2022–2024)

Se revisó el código original (`github.com/sierraporta/Fine-vs-Coarse-Dominance-from-EPA-PM2.5-PM10-2022-2024`).

1. **El FDI implementado no es el de la Ec. 2.** El código usa `r12 = PM25/(PM10+0.1)`, `r0 = (PM10−PM25+1)/(PM10+PM25+0.1)` y aplica `abs()` a concentraciones negativas. Por eso la Fig. 2 muestra FDI negativos, imposibles con la fórmula publicada.
2. **Filtros no declarados:** PM10 > 250 y PM2.5 > 50 µg/m³ se ponían como faltantes (eliminaba episodios de humo y polvo). La completitud real era "≤5 % de faltantes", no "≥1000 días por parámetro".
3. **Peso del IFS:** el texto dice w_FDI = 1,5; el código usa 1,0. Además el signo positivo de asimetría y curtosis contradice la interpretación ("ráfagas gruesas atenuadas").
4. **Los "hotspots" de verano residencial urbano/suburbano (IFS 4,33 y 3,18) eran artefactos:** curtosis de r0 de 604 y 724, producida por días con concentraciones casi nulas (r0 → 10 con el +1/+0,1).
5. **Celdas delgadas:** el texto dice que se excluyen de las tablas, pero el código calcula la bandera después de ordenar. 4 de las 20 celdas de la Tabla 3 tenían n < 300, y la primera (Commercial Rural Winter) tenía el IC demasiado ancho.
6. **Figura equivocada:** se incluyó `fig1_fdi_heatmap.png` (marginales con título "FDI") en lugar de `fig2_fdi_heatmap.png`, que es la que coincide con el pie y con las tablas.
7. **Inconsistencias de texto:** se mencionaban celdas agrícolas y desérticas suburbanas que no existen, y "desierto urbano siempre negativo" cuando el invierno era +0,13.
8. **Inferencia:** las pruebas sobre ~240 000 días autocorrelacionados (pseudo-replicación) y con Mobile incluido en la Tabla 4, aunque los Métodos dicen que se excluye. El análisis de sensibilidad de pesos comparaba contra la media del FDI, no contra el IFS, y no se hizo Holm en las pruebas ómnibus.
9. **Otros:** 1000/1095 en lugar de 1000/1096 (2024 es bisiesto); "is treated" → "are treated"; DOIs corruptos en el .bib (p. ej. `/SUPPL_FILE/...AVI`); citas fuera de contexto (Nasr2025 para Mann–Whitney, Kala2024 para delta de Cliff, Allajbeu2017 para estacionalidad).
10. **Superposición textual con el artículo publicado en CLEAN (doi:10.1002/clen.70191, mismos autores):** varios párrafos eran casi idénticos y el artículo no se citaba. Se reescribieron y ahora se cita y se diferencia explícitamente.

## 3. Qué cambió en el manuscrito

- **Período unificado 2023–2025** (el de los resultados nuevos) para ambas capas.
- **Pipeline nuevo y determinista** (`paper_pipeline/`):
  - Un valor por sitio-parámetro-día: se usa el POC más completo, se prioriza "24-HR BLK AVG" y se descartan las filas "Excluded" y las que tienen < 75 % de horas.
  - Completitud ≥ 75 % (822/1096 días); se excluyen días con PM2.5 > PM10 (1,6 %).
  - FDI acotado en (0, 1]; bootstrap por bloques de estación; celdas delgadas excluidas de verdad de las tablas.
  - Pruebas sobre unidades estación-estación-año, Holm en todas las familias, ε², δ de Cliff y Friedman/W de Kendall para la estacionalidad.
  - Capa de complejidad corregida: ventanas estacionales contiguas, desempate aleatorio, 199 sustitutos barajados por ventana, plano complejidad–entropía con D = 4 y pruebas pareadas entre fracciones.
  - Clustering con 4 métricas no redundantes, más una prueba de estabilidad.
  - Vínculo FDI ↔ complejidad con correlaciones parciales.
- **9 figuras** (antes 3) y **10 tablas** (antes 5). Todos los números de las tablas provienen de `outputs_paper/latex/`, verificados automáticamente.
- **Bibliografía depurada** (DOIs verificados en doi.org) y **27 referencias nuevas**, todas verificadas.
- Se eliminó `longnamesfirst`, que listaba todos los autores en la primera cita. El resumen tiene unas 250 palabras y los highlights menos de 85 caracteres cada uno.

## 4. Resultados principales (nuevos)

- **Estacionalidad:** el FDI es máximo en invierno en el 83 % de las estaciones (W = 0,42). Los hotspots robustos son celdas de invierno en entornos construidos suburbanos y urbanos; los coldspots son sitios industriales rurales y desierto.
- **Uso de suelo:** con inferencia correcta, las diferencias son pequeñas salvo bosque (más fino) e industrial suburbano en verano. Esto coincide con el artículo de CLEAN.
- **Complejidad:** las series son ruido correlacionado de alta entropía, pero más ordenadas que los sustitutos. La entropía ordena PM2.5 < PM10 < PMc. El verano es más suave y de colas más pesadas (humo 2023).
- **Uso de suelo y complejidad:** explica ε² ≈ 0,13–0,15 de la dinámica de PM10 (desierto, industrial y agrícola son más aleatorios y con ráfagas), pero casi nada de PM2.5, que se organiza por geografía.
- **Vínculo entre capas:** las estaciones con mayor FDI tienen un PM2.5 más estructurado y un PM10 con menos ráfagas (|ρ| ≈ 0,31–0,35, robusto a lat/lon).

## 5. Decisiones pendientes de los autores

1. **Declaración de IA generativa:** Elsevier la exige si se usó IA en la redacción. Se dejó un texto estándar al final del manuscrito; revísenlo o ajústenlo.
2. **Conflicto de intereses:** se añadió "no competing interests" (igual que en CLEAN); confírmenlo.
3. **Disponibilidad de código:** el manuscrito apunta a `github.com/IsabellaArrieta/atmospheric-pollutants-analysis`. Hay que hacer push de `paper_pipeline/`, `outputs_paper/` y `docs/` para que el enlace sea válido.
4. **Título:** se amplió para reflejar la capa de complejidad; si la revista no permite cambiarlo, puede volverse al original.
5. **Carta de respuesta:** si es una revisión, conviene explicar el cambio de período (2022–2024 → 2023–2025) y la corrección del FDI/IFS. Los puntos de la sección 2 sirven de base.
6. **Archivos sin uso en `docs/.../`:** `fig1_fdi_heatmap.png`, `fig2_fdi_box_facets.png`, `fig2_fdi_heatmap.png`, `Stations.pdf` y las figuras de la plantilla (`figs/Fig1-3.pdf`, `grabs.pdf`, `pic1.pdf`) ya no se usan. Sus copias están en `original_submission/` y pueden retirarse del paquete de envío.
