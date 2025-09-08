# Proyecto-Terminal

## Última búsqueda
- Término: ("IBS") AND (Systematic Review[Publication Type])
- Artículos guardados: 200

Generado el 2025-09-02 15:35:31 UTC con versión 0.2.

## Environment Variables

Set the `NCBI_API_KEY` environment variable to use an NCBI API key for PubMed requests:

```bash
export NCBI_API_KEY="your_api_key_here"
```

This increases the request rate limits when running `PubMed_API_0.1.py`.

## Instructivo del pipeline de análisis de miRNA

### 1. Preparación del entorno
1. **Clona o descarga** este repositorio y sitúate en su raíz (`Proyecto-Terminal/`).
2. **Variable de entorno para PubMed** (opcional pero recomendado):
   ```bash
   export NCBI_API_KEY="tu_api_key"
   ``` 
   Mejora los límites de uso al consultar PubMed

3. **Herramientas externas** a tener en el `PATH` (según cada etapa):
   - `esearch`, `elink`, `efetch`, `esummary`, `xtract` (EDirect/Entrez)
   - `prefetch`, `fasterq-dump` (SRA Toolkit)  
   - `jq`, `pigz` (opcional), `awk`, `sort`, `xargs`
   - `fastqc`, `cutadapt`, `multiqc`
   - `mirge3` o `miRDeep2` + `bowtie`, `bowtie-build`
   - Python 3 con `pandas`; R con paquetes `DESeq2`, `clusterProfiler`, `multiMiR`, `edgeR`, etc.

> **Tip:** crea directorios `data_raw/`, `data_ref/`, `results/` y `outputs/` si no existen, ya que varias etapas los usan por defecto.

---

### 2. 00_search – Búsqueda y resolución GEO → SRA
```bash
bash pipelines/00_search/geo_to_sra.sh "<tu_consulta>"
```
- Sin argumentos usa la consulta por defecto para enfermedades inflamatorias y miRNA en humano

.
- Genera:
  - `outputs/run_tables/gds_to_sra_runinfo.csv`
  - `outputs/search/datasets.jsonl`
  - `outputs/search/provenance.txt`, entre otros.

**Ejemplo**
```bash
bash pipelines/00_search/geo_to_sra.sh \
'("ulcerative colitis"[All Fields]) AND (microRNA) AND "Homo sapiens"[Organism]'
```

---

### 3. 01_download – Descarga de corridas SRA
```bash
bash pipelines/01_download/download_runs.sh \
  <runinfo.csv|lista_SRR> <dir_fastq> [hilos_fasterq] [trabajos_en_paralelo]
```
- Parámetros principales

.
- Usa `prefetch` → `fasterq-dump` y comprime en `.fastq.gz`.

**Ejemplo**
```bash
bash pipelines/01_download/download_runs.sh \
  outputs/run_tables/gds_to_sra_runinfo.csv data_raw/sra_fastq 8 4
```

---

### 4. 02_qc – Trimming y control de calidad
```bash
bash pipelines/02_qc/trim_qc.sh \
  <INPUT_DIR> <TRIM_DIR> <QC_DIR> [ADAPTER] [MINLEN] [THREADS]
```
- Parámetros y defaults

.
- FastQC antes/después, recorta adaptador 3' con `cutadapt`, genera reporte `MultiQC`.

**Ejemplo**
```bash
bash pipelines/02_qc/trim_qc.sh \
  data_raw/sra_fastq data_raw/trimmed outputs/qc TGGAATTCTCGGGTGCCAAGG 15 8
```

---

### 5. 03_quant – Cuantificación de miRNA  
#### Opción A: miRge3
```bash
bash pipelines/03_quant/mirge3_run.sh \
  <TRIM_DIR> <OUT_DIR> [SPECIES] [ADAPTER] [THREADS]
```
- Defaults: especie `hsa`, adaptador Illumina, 4 hilos

.
- Salida estandarizada: `results/quant_mirge3/<muestra>/miRNA.counts.tsv`.

**Ejemplo**
```bash
bash pipelines/03_quant/mirge3_run.sh \
  data_raw/trimmed results/quant_mirge3 hsa TGGAATTCTCGGGTGCCAAGG 8
```

#### Opción B: miRDeep2
```bash
bash pipelines/03_quant/mirdeep2_run.sh \
  <TRIM_DIR> <OUT_DIR> <GENOME_FA> <BOWTIE_INDEX_PREFIX> \
  <MAT_FA> <HAIRPIN_FA> [SPECIES] [MINLEN]
```
- Permite especificar genoma, índice Bowtie y archivos miRBase

.
- Salida estandarizada en `results/quant_mirdeep2/<muestra>/miRNA.counts.tsv`.

**Ejemplo**
```bash
bash pipelines/03_quant/mirdeep2_run.sh \
  data_raw/trimmed results/quant_mirdeep2 \
  data_ref/human_ds/genome.fa data_ref/bowtie_index/human/genome \
  data_ref/mirbase/mature_hsa.fa data_ref/mirbase/hairpin_hsa.fa hsa 15
```

---

### 6. 04_merge – Fusión de conteos por muestra
```bash
python pipelines/04_merge/merge_counts.py \
  <quant_dir> <metadata_csv> <out_csv> \
  [--prefer-mirna-col COL] [--prefer-count-col COL] \
  [--counts-filename FILE] [--min-rows N]
```
- Argumentos principales y opciones

.
- Crea una matriz (`miRNA` × `muestra`) en `data/data_processed/miRNA_counts_raw.csv`.

**Ejemplo**
```bash
python pipelines/04_merge/merge_counts.py \
  results/quant_mirge3 data/metadata/sample_metadata.csv \
  data/data_processed/miRNA_counts_raw.csv
```

---

### 7. 05_dea_r – Análisis diferencial (DESeq2, opcional edgeR/voom)
```bash
Rscript pipelines/05_dea_r/run_dea.R \
  <counts_csv> <meta_csv> [control_label] [padj_cut] [lfc_cut] [run_voom]
```
- Defaults y argumentos

.
- Produce tablas DESeq2 (`DEA_DESeq2_*_vs_*`) y figuras (volcano, heatmap) en `results/DEA/`.

**Ejemplo**
```bash
Rscript pipelines/05_dea_r/run_dea.R \
  data/data_processed/miRNA_counts_raw.csv \
  data/metadata/sample_metadata.csv Control 0.05 1.0 TRUE
```

---

### 8. 06_targets_enrich – Targets y análisis de enriquecimiento
#### 8.1 Obtener interacciones miRNA–gen con **multiMiR**
```bash
Rscript pipelines/06_targets_enrich/multimir_targets.R \
  [dea_csv_o_directorio] [padj_cut] [lfc_cut] [top_n] [org]
```
- Uso y parámetros

.
- Genera `results/targets_validated.csv`, `results/targets_predicted.csv`, etc.

**Ejemplo**
```bash
Rscript pipelines/06_targets_enrich/multimir_targets.R \
  results/DEA 0.05 1.0 30 hsa
```

#### 8.2 Enriquecimiento funcional (KEGG/GO/Reactome)
```bash
Rscript pipelines/06_targets_enrich/enrichment.R \
  [targets_file] [organism_kegg] [p_cut] [q_cut] [minGS] [maxGS] [background_file]
```
- Argumentos y defaults

.
- Resultados en `results/enrichment/` con tablas y gráficos.

**Ejemplo**
```bash
Rscript pipelines/06_targets_enrich/enrichment.R \
  results/validated_target_genes_unique.csv hsa 0.05 0.05 10 500
```

---

### 9. 07_networks_clin – Redes y evaluación clínica
#### 9.1 Construir red miRNA–gen
```bash
Rscript pipelines/07_networks_clin/build_network.R \
  [validated_csv] [predicted_csv] [include_predicted] \
  [min_evidence] [min_degree] [mirna_list_file] [string_ppi_file]
```
- Argumentos y defaults

.
- Produce archivos de aristas/nodos y figuras en `results/networks/`.

**Ejemplo**
```bash
Rscript pipelines/07_networks_clin/build_network.R \
  results/targets_validated.csv results/targets_predicted.csv FALSE 1 1
```

#### 9.2 Evaluación clínica de miRNAs
```bash
Rscript pipelines/07_networks_clin/clinical_eval.R \
  [norm_counts_csv] [meta_csv] [mirna_list_file] \
  [outcome_col] [positive_label] [top_n_fallback] [numeric_vars_csv]
```
- Opciones disponibles

.
- Salidas en `results/clinical/` con correlaciones, ROC/AUC y figuras.

**Ejemplo**
```bash
Rscript pipelines/07_networks_clin/clinical_eval.R \
  results/DEA/normalized_counts_deseq2.csv \
  data/metadata/sample_metadata.csv \
  "" "disease" "Case" 30
```

---

### 10. Flujo completo sugerido
1. `00_search` → genera `runinfo.csv`.
2. `01_download` → descarga `.fastq.gz`.
3. `02_qc` → trimming + QC.
4. `03_quant` → conteos por miRNA.
5. `04_merge` → matriz de conteos unificada.
6. `05_dea_r` → miRNAs diferenciales.
7. `06_targets_enrich` → objetivos génicos y enriquecimiento.
8. `07_networks_clin` → red interactiva y evaluación clínica.

Este instructivo cubre la ejecución básica del pipeline; cada script incluye más comentarios y logs que pueden consultarse para personalizaciones adicionales. ¡Éxito en tu análisis!
