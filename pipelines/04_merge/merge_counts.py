#!/usr/bin/env python3
# ===========================================
# Script: 04_merge/merge_counts.py
# Descripción: Fusiona los conteos por muestra (miRge3/miRDeep2) en una matriz única.
#              Mapea SRR→sample_id usando metadata/sample_metadata.csv si existe.
#              Estándar: filas = miRNA, columnas = muestras.
# Fecha: $(date +%Y-%m-%d)
# ===========================================

import sys
import argparse
from pathlib import Path
import csv
import re
from datetime import datetime

try:
    import pandas as pd
except ImportError:
    print("[ERROR] pandas no está instalado. Ejecuta: pip install pandas", file=sys.stderr)
    sys.exit(127)

# --- Variables (por defecto relativas al repo PT) ---
ROOT_DIR = Path(__file__).resolve().parents[2]  # .../PT
LOG_DIR = ROOT_DIR / "outputs" / "logs"
LOGFILE = LOG_DIR / "04_merge.log"

# --- Logging básico (a consola + archivo) ---
def log(*args):
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    msg = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] " + " ".join(map(str, args))
    print(msg)
    with open(LOGFILE, "a", encoding="utf-8") as fh:
        fh.write(msg + "\n")

# --- Utilidades ---
def sanitize_cols(cols):
    out = []
    for c in cols:
        c2 = re.sub(r"[^A-Za-z0-9_]+", "", c.strip().lower())
        out.append(c2)
    return out

def sniff_sep_and_read(path: Path):
    """
    Intenta leer TSV/CSV/TXT. Primero como TSV. Si solo hay una columna y contiene comas, reintenta como CSV.
    Devuelve DataFrame y separador usado.
    """
    try:
        df = pd.read_csv(path, sep="\t", engine="python")
        if df.shape[1] == 1:
            # ¿parece CSV?
            if df.columns.size == 1 and ("," in df.columns[0] or df.iloc[:50, 0].astype(str).str.contains(",").any()):
                df = pd.read_csv(path, sep=",", engine="python")
                return df, ","
        return df, "\t"
    except Exception:
        # como CSV
        df = pd.read_csv(path, sep=",", engine="python")
        return df, ","

def find_counts_file(sample_dir: Path, preferred: str = "miRNA.counts.tsv"):
    """
    Busca el archivo de conteos dentro del directorio de muestra.
    Prioriza el nombre estándar 'miRNA.counts.tsv'.
    Si no existe, intenta patrones comunes (cauto para evitar falsos positivos).
    """
    cand = sample_dir / preferred
    if cand.exists():
        return cand

    patterns = [
        "*miR*count*.tsv", "*mir*count*.tsv", "*miR*reads*.tsv", "*mir*.tsv",
        "*miR*count*.txt", "*reads*per*miRNA*.*", "miRNAs_expressed_*.*",
        "*.csv"
    ]
    for pat in patterns:
        hits = list(sample_dir.glob(pat))
        if hits:
            # Toma el primero determinísticamente (ordenado)
            return sorted(hits)[0]
    return None

def load_run_to_sample_map(meta_csv: Path):
    """
    Carga metadata y devuelve dict run_accession -> sample_id.
    Si no existe o faltan columnas, devuelve {}.
    """
    if not meta_csv.exists():
        log(f"[warn] Metadata no encontrada: {meta_csv} — usaré nombres de carpeta como sample_id.")
        return {}
    try:
        dfm = pd.read_csv(meta_csv)
    except Exception as e:
        log(f"[warn] No pude leer metadata ({meta_csv}): {e} — usaré nombres de carpeta.")
        return {}
    cols = {c.lower(): c for c in dfm.columns}
    if "run_accession" in cols and "sample_id" in cols:
        mapping = dict(zip(dfm[cols["run_accession"]].astype(str), dfm[cols["sample_id"]].astype(str)))
        # También admite una columna 'sample' si existe:
        return mapping
    else:
        log("[warn] Metadata no contiene columnas 'run_accession' y 'sample_id'. Usaré nombres de carpeta.")
        return {}

def infer_sample_name(dir_name: str):
    """
    Normaliza nombres de muestra a partir del nombre del directorio.
    Quita sufijos comunes como _1/_R1.
    """
    base = re.sub(r"(_R?1)$", "", dir_name, flags=re.IGNORECASE)
    return base

def choose_mirna_and_count_cols(df: pd.DataFrame, prefer_mirna=None, prefer_count=None):
    """
    Selecciona columnas de miRNA y conteo. Permite overrides.
    """
    orig_cols = list(df.columns)
    low = sanitize_cols(orig_cols)

    # override si el usuario lo indicó
    if prefer_mirna and prefer_mirna in orig_cols:
        mir_col = prefer_mirna
    else:
        candidates_miR = ["mirna", "mir", "mirnaid", "name", "id", "feature", "gene"]
        mir_col = None
        for i, lc in enumerate(low):
            if lc in candidates_miR or "mir" in lc:
                mir_col = orig_cols[i]
                break
        if mir_col is None:
            mir_col = orig_cols[0]

    if prefer_count and prefer_count in orig_cols:
        cnt_col = prefer_count
    else:
        candidates_cnt = ["count", "counts", "readcount", "reads", "totalreads", "abundance", "expression"]
        cnt_col = None
        for i, lc in enumerate(low):
            if lc in candidates_cnt:
                cnt_col = orig_cols[i]
                break
        if cnt_col is None:
            # intenta segunda columna
            cnt_col = orig_cols[1] if len(orig_cols) > 1 else orig_cols[0]

    return mir_col, cnt_col

def coerce_counts_series(s):
    try:
        return pd.to_numeric(s, errors="coerce").fillna(0).astype(int)
    except Exception:
        return pd.to_numeric(s, errors="coerce").fillna(0).astype(float)

# --- CLI ---
parser = argparse.ArgumentParser(
    description="Merge per-sample miRNA counts into a single matrix (rows=miRNA, cols=samples)."
)
parser.add_argument("quant_dir", nargs="?", default=str(ROOT_DIR / "results" / "quant_mirge3"),
                    help="Directory with per-sample subfolders (miRge3/miRDeep2 outputs).")
parser.add_argument("metadata_csv", nargs="?", default=str(ROOT_DIR / "data" / "metadata" / "sample_metadata.csv"),
                    help="Metadata CSV with columns run_accession and sample_id (optional).")
parser.add_argument("out_csv", nargs="?", default=str(ROOT_DIR / "data" / "data_processed" / "miRNA_counts_raw.csv"),
                    help="Output counts matrix CSV.")
parser.add_argument("--prefer-mirna-col", default=None, help="Preferred miRNA column name if known.")
parser.add_argument("--prefer-count-col", default=None, help="Preferred count column name if known.")
parser.add_argument("--counts-filename", default="miRNA.counts.tsv",
                    help="Exact counts filename expected in each sample directory (default: miRNA.counts.tsv).")
parser.add_argument("--min-rows", type=int, default=1, help="Minimum rows required per sample table to include it.")
args = parser.parse_args()

def main():
    # Inicio
    LOGFILE.write_text("")  # reset log
    log("Iniciando script...")
    log(f"QUANT_DIR = {args.quant_dir}")
    log(f"METADATA  = {args.metadata_csv}")
    log(f"OUT_CSV   = {args.out_csv}")

    quant_dir = Path(args.quant_dir)
    meta_csv = Path(args.metadata_csv)
    out_csv = Path(args.out_csv)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    if not quant_dir.exists():
        log(f"[ERROR] Directorio de cuantificación no existe: {quant_dir}")
        sys.exit(1)

    # subdirectorios = muestras
    samples = [p for p in sorted(quant_dir.iterdir()) if p.is_dir()]
    if not samples:
        log(f"[ERROR] No se encontraron subdirectorios de muestra en {quant_dir}")
        sys.exit(1)
    log(f"Subdirectorios detectados (muestras): {len(samples)}")

    run2sample = load_run_to_sample_map(meta_csv)

    tables = []
    included = 0
    skipped = 0

    for smp_dir in samples:
        smp_name = infer_sample_name(smp_dir.name)
        counts_path = find_counts_file(smp_dir, preferred=args.counts_filename)
        if counts_path is None or not counts_path.exists():
            log(f"[warn] No counts file in {smp_dir} (busqué {args.counts_filename} y patrones comunes).")
            skipped += 1
            continue

        try:
            df_raw, used_sep = sniff_sep_and_read(counts_path)
        except Exception as e:
            log(f"[warn] Error leyendo {counts_path}: {e}")
            skipped += 1
            continue

        if df_raw.shape[0] < args.min_rows or df_raw.shape[1] < 1:
            log(f"[warn] Tabla muy pequeña en {counts_path} — omitida.")
            skipped += 1
            continue

        mir_col, cnt_col = choose_mirna_and_count_cols(df_raw, args.prefer_mirna_col, args.prefer_count_col)
        if mir_col not in df_raw.columns or cnt_col not in df_raw.columns:
            log(f"[warn] Columnas no encontradas en {counts_path} (miR={mir_col}, count={cnt_col}) — omitida.")
            skipped += 1
            continue

        df = df_raw[[mir_col, cnt_col]].copy()
        df.columns = ["miRNA", "count"]
        # agrupa por miRNA por si hay duplicados
        df = df.groupby("miRNA", as_index=False)["count"].sum()

        # nombre final de muestra: mapear SRR→sample_id si procede
        final_name = run2sample.get(smp_name, smp_name)

        s = df.set_index("miRNA")["count"]
        s = coerce_counts_series(s)
        s.name = final_name
        tables.append(s)
        included += 1
        log(f"[ok] {smp_dir.name} → {final_name} | filas={df.shape[0]} | sep='{used_sep}' | archivo={counts_path.name}")

    if not tables:
        log("[ERROR] No se pudo construir ninguna columna de conteos. Revisa entradas.")
        sys.exit(1)

    mat = pd.concat(tables, axis=1).fillna(0)
    # ordena columnas alfabéticamente para estabilidad, pero respeta el orden del metadata si todos están presentes
    if run2sample:
        # intenta ordenar por sample_metadata
        try:
            meta = pd.read_csv(meta_csv)
            if "sample_id" in meta.columns:
                order = [sid for sid in meta["sample_id"].astype(str).tolist() if sid in mat.columns]
                rest = [c for c in mat.columns if c not in order]
                mat = mat[order + rest]
        except Exception:
            pass

    # asegura enteros cuando es posible
    try:
        mat = mat.astype(int)
    except Exception:
        # puede haber floats si hubo NaN previos
        pass

    mat.to_csv(out_csv, quoting=csv.QUOTE_MINIMAL)
    log(f"[done] Matriz escrita: {out_csv} | shape={mat.shape[0]}x{mat.shape[1]}")
    log(f"Samples incluidos: {included} | omitidos: {skipped}")

    # Fin
    log("Script finalizado correctamente.")

if __name__ == "__main__":
    main()
