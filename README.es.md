[English](README.md) | **Español**

# Pipeline de meta-análisis de miRNA-seq en Enfermedad Inflamatoria Intestinal (IBD)

Pipeline completo de extremo a extremo para el análisis de datos de secuenciación de
microRNA (miRNA-seq) en **enfermedad inflamatoria intestinal** (colitis ulcerosa y
enfermedad de Crohn) en *Homo sapiens*. Cubre desde el descubrimiento de datasets en
GEO/SRA hasta el reporte final interactivo y un tablero (dashboard) comparativo entre
múltiples estudios.

Proyecto terminal de Ingeniería en Biotecnología.

---

## Tabla de contenido

1. [Visión general](#1-visión-general)
2. [`datasets.yaml` — registro maestro](#2-datasetsyaml--registro-maestro)
3. [Preparación del entorno](#3-preparación-del-entorno)
4. [Etapas del pipeline](#4-etapas-del-pipeline)
   - [00_search — Descubrimiento y curación GEO/SRA](#00_search--descubrimiento-y-curación-geosra)
   - [01_download — Descarga de corridas SRA](#01_download--descarga-de-corridas-sra)
   - [02_qc — Control de calidad y organización](#02_qc--control-de-calidad-y-organización-de-muestras)
   - [03_quant — Cuantificación de miRNA (miRge3)](#03_quant--cuantificación-de-mirna-mirge3)
   - [04_merge — Fusión, filtrado y metadatos](#04_merge--fusión-filtrado-y-metadatos)
   - [05_dea_r — Expresión diferencial (DESeq2)](#05_dea_r--expresión-diferencial-deseq2)
   - [06_targets_enrich — Predicción de blancos (multiMiR)](#06_targets_enrich--predicción-de-blancos-multimir)
   - [07_networks_clin — Enriquecimiento, redes y evaluación clínica](#07_networks_clin--enriquecimiento-redes-y-evaluación-clínica)
   - [08_report — Reporte, dashboard y diagrama](#08_report--reporte-dashboard-y-diagrama)
5. [Librería compartida: kits de preparación](#5-librería-compartida-kits-de-preparación)
6. [Módulo paralelo: motor de búsqueda de literatura](#6-módulo-paralelo-motor-de-búsqueda-de-literatura)
7. [Flujo completo sugerido](#7-flujo-completo-sugerido)

---

## 1. Visión general

El pipeline está organizado en etapas numeradas dentro de `pipelines/`. Cada etapa
produce salidas que alimentan a la siguiente, todas bajo `outputs/<id_dataset>/`.

| Etapa | Carpeta | Descripción | Entorno |
|---|---|---|---|
| 00 | `00_search` | Descubrimiento y curación de datasets en GEO/SRA | `py_env` |
| 01 | `01_download` | Descarga de corridas SRA (prefetch + fasterq-dump) | `py_env` |
| 02 | `02_qc` | QC con FastQC/MultiQC + miRTrace y organización de muestras | `py_env` |
| 03 | `03_quant` | Cuantificación de miRNA con miRge3 | `py_env` |
| 04 | `04_merge` | Fusión de conteos, filtrado por RPM y metadatos DESeq2 | `py_env` |
| 05 | `05_dea_r` | Análisis de expresión diferencial (DESeq2) | `r_env` |
| 06 | `06_targets_enrich` | Predicción de blancos miRNA→gen (multiMiR) | `r_env` |
| 07 | `07_networks_clin` | Enriquecimiento GO/KEGG, redes, HMDD, solapamiento IBD y ROC clínico | `r_env` |
| 08 | `08_report` | Reporte HTML, dashboard Streamlit y diagrama del pipeline | `r_env` / Python |

Además existe un módulo paralelo de minería de literatura
(`pipelines/search_engine_papers/`) independiente del flujo principal.

**Datasets analizados actualmente** (ver `datasets.yaml`): `GSE84779`, `PRJNA471862`,
`PRJNA717018`, `GSE114603`, `GSE272890`.

---

## 2. `datasets.yaml` — registro maestro

`datasets.yaml` (en la raíz) es la **única fuente de verdad** del proyecto: un registro
de cada dataset con su accession, BioProject, PMID, tejido, comparaciones, número de
muestras, rutas de salida y `status` (`analyzed` / `in_progress` / `suspended`).

Sirve para que nada tenga que inferirse de los prefijos de carpeta. Alimenta tanto al
**dashboard** como al **reporte maestro** entre datasets. Al añadir un dataset, se copia
un bloque y se rellenan los campos conocidos (los desconocidos se dejan en `null`).

---

## 3. Preparación del entorno

### Variable de entorno para PubMed (opcional)

```bash
export NCBI_API_KEY="tu_api_key"
```

Mejora los límites de uso al consultar PubMed (usado por el motor de búsqueda de
literatura y por los curadores de la etapa 00).

### Entornos conda (WSL / Linux)

El pipeline corre en **WSL** y usa **dos entornos** (los archivos están en `envs/`):

| Entorno | Etapas | Contenido |
|---|---|---|
| **`py_env`** | 00–04 | miRge3, miRTrace, sra-tools, entrez-direct, FastQC, MultiQC, seqkit, pandas (todo lo Python/CLI/Java). Fijado a **Python 3.8** por miRge3. |
| **`r_env`** | 05–08 | R 4.x + Bioconductor: DESeq2, apeglm, multiMiR, clusterProfiler, enrichplot, org.Hs.eg.db, fgsea, tidyverse, igraph, pROC, etc. |

```bash
conda env create -f envs/py_env.yml   # o: conda env update -n py_env -f envs/py_env.yml
conda env create -f envs/r_env.yml
```

Los scripts resuelven sus rutas (`miRge3_Lib`, `mirtrace.jar`) vía `$CONDA_PREFIX`, así
que basta con tener activo el entorno correcto.

### Dependencias Python para Windows / dashboard

`requirements.txt` lista **solo** los paquetes pip usados fuera de los entornos conda
(p. ej. el dashboard Streamlit y el motor de búsqueda en Windows):

```bash
pip install -r requirements.txt
```

> ⚠️ Se fija `pandas<3` porque `upsetplot 0.9.0` (la última versión) se rompe bajo el
> *Copy-on-Write* de pandas 3.0 — ver `pipelines/08_report/dashboard/upset_compat.py`.

### Verificación de herramientas

Antes de empezar, comprueba que toda la cadena de herramientas esté disponible:

```bash
bash pipelines/00_search/00_tool_check.sh
```

Sale con código 0 si todo lo requerido está presente; en caso contrario sale con 1 y
escribe un log con lo que falta.

---

## 4. Etapas del pipeline

### 00_search — Descubrimiento y curación GEO/SRA

Conjunto de curadores que descubren datasets, extraen metadatos y generan fichas
Markdown en `outputs/docs/`, que luego sincronizan `datasets.yaml`.

| Script | Función |
|---|---|
| `00_tool_check.sh` | Verifica la cadena de herramientas completa de ambos entornos conda. |
| `01_geo_ibd_curator.py` | Descubrimiento automático: busca, filtra y extrae metadatos de datasets GEO/SRA mediante *true pivot* sobre BioProject. Usa `ignore_list.txt` como lista negra y genera listas de SRR para descarga. |
| `02_bioproject_curator.py` | Curación dirigida ("sniper") a partir de uno o varios BioProject (`PRJNA…`), saltándose GEO por completo. |
| `03_geo_id_curator.py` | Curación dirigida a accesiones GEO explícitas (`GSE…`); mismo formato de salida que `01_` pero sin el paso de búsqueda amplia. |
| `04_dedup_docs.py` | Elimina fichas `PRJNA*.md` redundantes cuando ya existe un `GSE*.md` para el mismo BioProject (dry-run por defecto; `--confirm` para borrar). |
| `05_sync_yaml.py` | Sincroniza `datasets.yaml` desde `outputs/docs/*.md`: rellena campos `null` sin sobrescribir valores existentes. Con `--add-new` añade entradas nuevas. Requiere `ruamel.yaml`. |

```bash
# Descubrimiento automático
python pipelines/00_search/01_geo_ibd_curator.py

# Curación dirigida por accesión GEO o BioProject
python pipelines/00_search/03_geo_id_curator.py --id GSE114603
python pipelines/00_search/02_bioproject_curator.py --id PRJNA904924

# Limpiar duplicados y sincronizar el manifiesto
python pipelines/00_search/04_dedup_docs.py --confirm
python pipelines/00_search/05_sync_yaml.py --add-new --confirm
```

> La subcarpeta `complementary_scripts/` contiene utilidades legadas/auxiliares
> (`archive_geo_to_sra.sh`, post-procesado de runinfo, obtención de atributos de
> BioSample, etc.).

---

### 01_download — Descarga de corridas SRA

```bash
bash pipelines/01_download/03_download_runs.sh \
  [runinfo.csv|srr_list.txt] [dir_fastq] [hilos_fasterq] [trabajos_paralelos]
```

- Sin argumentos busca `srr_list.txt` en el directorio actual y descarga a `sra_fastq/`.
- Usa `prefetch` → `fasterq-dump` y comprime a `.fastq.gz`; deja log local `01_download.log`.

La subcarpeta `complementary_scripts/` incluye variantes y utilidades:
corrección de nombres de muestras (`011_correction_name_samples.sh`) y compresión
(`012_compress_fastq.sh`).

---

### 02_qc — Control de calidad y organización de muestras

| Script | Función |
|---|---|
| `04_qc_fastqc.sh` | Ejecuta FastQC sobre todos los `*.fastq.gz` y consolida con MultiQC. |
| `05_create_config_mirtrace.sh` | Genera `mirtrace_config.csv` eligiendo el adaptador 3' según el kit del estudio (ver `kits.tsv`). |
| `06_qc_mirtrace.sh` | Corre miRTrace por lotes usando el protocolo/adaptador del kit. |
| `07_organize_samples.py` | Clasifica los FASTQ en `sra_fastq/muestras_clasificadas/` según el runinfo. |

```bash
bash pipelines/02_qc/04_qc_fastqc.sh
DATASET=GSE272890 bash pipelines/02_qc/05_create_config_mirtrace.sh
bash pipelines/02_qc/06_qc_mirtrace.sh mirtrace_config.csv
python pipelines/02_qc/07_organize_samples.py
```

---

### 03_quant — Cuantificación de miRNA (miRge3)

| Script | Función |
|---|---|
| `08_auto_make_list_subfolders.sh` | Genera la lista de subcarpetas de muestras a procesar. |
| `09_auto_subfolders_mirge3.sh` | Itera de forma segura y secuencial todas las subcarpetas de `muestras_clasificadas/`, corriendo miRge3 una a una. Evita rutas `/mnt/` de WSL (mucho más lentas). |
| `10_single_folder_mirge3.sh` | Procesa una sola carpeta / lista de corridas; controles `THREADS`, `PARALLEL`, `FASTQ_DIR`. |

```bash
# Desde la raíz del dataset
bash pipelines/03_quant/08_auto_make_list_subfolders.sh
bash pipelines/03_quant/09_auto_subfolders_mirge3.sh
```

Salida por muestra con `miR.Counts.csv` y `miR.RPM.csv`.

---

### 04_merge — Fusión, filtrado y metadatos

| Script | Función |
|---|---|
| `11_mrina_cut_off.py` | Busca recursivamente los `miR.Counts.csv` / `miR.RPM.csv`, los unifica en una matriz y filtra por **RPM ≥ 100**. Produce `Master_Filtered_Counts_DESeq2.csv`. |
| `12_deseq2_metadata.py` | Localiza el runinfo y la matriz maestra y genera el `sample_metadata.csv` alineado para DESeq2. |

```bash
python pipelines/04_merge/11_mrina_cut_off.py
python pipelines/04_merge/12_deseq2_metadata.py
```

---

### 05_dea_r — Expresión diferencial (DESeq2)

```bash
Rscript pipelines/05_dea_r/13_run_dea.R
```

- Encuentra automáticamente `Master_Filtered_Counts_DESeq2.csv` y `sample_metadata.csv`
  desde la raíz del dataset.
- Interactivo: muestra las columnas de metadatos disponibles y pregunta cuál analizar.
- Corre DESeq2 con `apeglm` y VST; produce dos modelos (**Default Replace** y
  **Strict No-Replace**), figuras (PCA, UMAP, heatmap, volcanos) y tablas de DEA.

---

### 06_targets_enrich — Predicción de blancos (multiMiR)

```bash
Rscript pipelines/06_targets_enrich/14_multimir_targets.R [dea_csv_o_dir] [padj] [lfc] [s_score_cut] [org]
# s_score_cut: piso opcional de S-score para la consulta de blancos. Default 0 = sin
# gate — consulta TODOS los miRNAs DE-significativos (solo padj/|log2FC|). Pasa p.ej. 2.0 para restringir.
```

- Obtiene interacciones miRNA→gen **validadas y predichas** con multiMiR a partir de los
  miRNA significativos del DEA.
- Las ordena por un score `S_miRNA`, divide los miRNA agrupados y consulta por lotes.
- Escribe las tablas en `multimir_outputs/`.

---

### 07_networks_clin — Enriquecimiento, redes y evaluación clínica

| Script | Función |
|---|---|
| `15_pathway_enrich.R` | Enriquecimiento GO (BP/CC/MF) y KEGG de los genes blanco con clusterProfiler; *smart file finder* para rutas vía CLI. Salida en `pathway_outputs/`. |
| `16_interactive_network.R` | Dashboards tripartitos interactivos miRNA → gen → vía/GO (Cytoscape.js). v3: filtra a interacciones *Functional MTI* fuertes y vista centrada en hubs para evitar la "maraña". |
| `17_treatment_response_roc.R` | Evaluación de desempeño como biomarcador (ROC/AUC) de miRNA candidatos: respondedores vs no respondedores a glucocorticoides, y UC vs controles sanos. Salida en `clinical_outputs/`. |
| `18_hmdd_prep.R` | Prepara listas de miRNA copiables para **HMDD v4.0** y **miRNet 2.0** y abre ambos sitios en el navegador. |
| `19_ibd_target_overlap.R` | Evalúa si los blancos validados están enriquecidos en genes de IBD conocidos vía la API de **OpenTargets** (GraphQL) + prueba exacta de Fisher por miRNA. |

```bash
Rscript pipelines/07_networks_clin/15_pathway_enrich.R
Rscript pipelines/07_networks_clin/16_interactive_network.R
Rscript pipelines/07_networks_clin/17_treatment_response_roc.R
Rscript pipelines/07_networks_clin/18_hmdd_prep.R
Rscript pipelines/07_networks_clin/19_ibd_target_overlap.R
```

---

### 08_report — Reporte, dashboard y diagrama

| Script | Función |
|---|---|
| `20_generate_report.py` | Genera el `Final_Report.html` interactivo (Plotly) de un dataset. Descubre automáticamente las subcarpetas de DEA (Default/Strict) y todas las comparaciones; integra QC, expresión diferencial, enriquecimiento, redes y solapamiento IBD. |
| `21_launch_dashboard.sh` | Lanza el dashboard Streamlit (`dashboard/Dashboard.py`). Busca un Python con Streamlit en los entornos conda del proyecto y sirve en `http://localhost:8501`. |
| `pipeline_diagram.py` | Genera el diagrama completo del pipeline para la tesis (`pipeline_diagram.png` a 300 DPI + `.pdf` vectorial). |

```bash
py -3.11 pipelines/08_report/20_generate_report.py <results_dir>
bash pipelines/08_report/21_launch_dashboard.sh
py -3.11 pipelines/08_report/pipeline_diagram.py
```

#### Dashboard comparativo (`dashboard/`)

Aplicación **Streamlit multi-página** dirigida por `datasets.yaml`: añadir un dataset al
manifiesto lo hace aparecer automáticamente.

- `Dashboard.py` — portada con métricas globales y tabla del registro de datasets.
- `data.py` — capa de acceso a datos; lee el manifiesto y la estructura de salidas en disco.
- `pages/1_Dataset_Explorer.py` — explorador por dataset: QC, expresión diferencial, enriquecimiento y solapamiento IBD.
- `pages/2_Cross_Dataset_Analysis.py` — análisis entre datasets: solapamiento (Venn/UpSet), matriz de direcciones, meta-análisis de Stouffer, concordancia de tamaños de efecto (Spearman), vías compartidas y solapamiento con genes IBD.
- `pages/3_miRNA_Lookup.py` — perfil por miRNA: expresión diferencial en todos los datasets/comparaciones más relevancia de solapamiento con dianas IBD y los genes IBD que regula.
- `upset_compat.py` — *shims* de compatibilidad para `upsetplot 0.9.0` bajo pandas/matplotlib modernos.

---

## 5. Librería compartida: kits de preparación

El **kit de preparación de librería** es una propiedad del estudio (cada corrida usa el
mismo prep), por lo que se asigna por dataset en `pipelines/kits.tsv` en lugar de
autodetectarlo de las lecturas.

- `pipelines/kits.tsv` — mapa `dataset_id → kit`. Kits conocidos: `illumina`,
  `nextflex_v3` (Illumina + UMIs 4,4), `nebnext`. Lo no listado usa `illumina`.
- `pipelines/lib/kit_params.sh` — librería que se hace `source` para resolver el kit
  (`resolve_kit`) y cargar sus parámetros (`load_kit_params`): adaptadores y protocolo
  para miRge3 y miRTrace. La usan los scripts de la etapa 02.

> Nota: NEXTflex v3 comparte el adaptador 3' exacto con Illumina estándar y solo difiere
> por los UMIs 4+4, que el *adapter-sniffing* no puede detectar — de ahí el mapeo
> explícito por dataset.

---

## 6. Módulo paralelo: motor de búsqueda de literatura

`pipelines/search_engine_papers/` es un módulo independiente para minería de literatura
científica (PubMed):

- `PubMed_API_0.1.py` — consultas a la API de PubMed (usa `NCBI_API_KEY`).
- `interactive_cli.py` — asistente interactivo para recopilar parámetros de búsqueda y ejecutar consultas.
- `pt_search.py`, `text_utils.py`, `paper_processing/refineria.py` — búsqueda TF-IDF y refinado/procesamiento de los documentos recuperados.

Consulta su propio `README.md` dentro de la carpeta para detalles.

---

## 7. Flujo completo sugerido

1. `00_search` → descubre/curiosa datasets y sincroniza `datasets.yaml`.
2. `01_download` → descarga `.fastq.gz`.
3. `02_qc` → QC (FastQC/MultiQC + miRTrace) y organiza muestras.
4. `03_quant` → conteos por miRNA con miRge3.
5. `04_merge` → matriz de conteos filtrada + metadatos DESeq2.
6. `05_dea_r` → miRNA diferencialmente expresados.
7. `06_targets_enrich` → blancos génicos con multiMiR.
8. `07_networks_clin` → enriquecimiento, redes, HMDD, solapamiento IBD y ROC clínico.
9. `08_report` → `Final_Report.html` por dataset + dashboard comparativo entre datasets.

Cada script incluye más comentarios y logs que pueden consultarse para
personalizaciones adicionales. ¡Éxito en tu análisis!
