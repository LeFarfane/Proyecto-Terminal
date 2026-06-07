import os
import re
import glob
import shutil
import pandas as pd

def encontrar_sra_fastq(desde):
    """Busca la carpeta sra_fastq/ desde el directorio indicado."""
    candidato = os.path.join(desde, 'sra_fastq')
    if os.path.isdir(candidato):
        return candidato
    for root, dirs, _ in os.walk(desde):
        if 'sra_fastq' in dirs:
            return os.path.join(root, 'sra_fastq')
    return None

def encontrar_csv(directorio):
    """Busca *_runinfo_merged.csv en el directorio indicado."""
    coincidencias = glob.glob(os.path.join(directorio, '*_runinfo_merged.csv'))
    return coincidencias[0] if coincidencias else None

def clasificar_fastq():
    print("🧬 Organizador de Muestras SRA 🧬\n")

    cwd = os.getcwd()

    ruta_base = encontrar_sra_fastq(cwd)
    if not ruta_base:
        print(f"❌ No se encontró la carpeta sra_fastq/ bajo: {cwd}")
        return
    print(f"📂 Carpeta de lecturas: {ruta_base}")

    # Buscar el CSV en el directorio padre de sra_fastq (donde suele estar)
    directorio_csv = os.path.dirname(ruta_base)
    archivo_excel = encontrar_csv(directorio_csv)
    if not archivo_excel:
        # Intentar también en cwd por si la estructura es distinta
        archivo_excel = encontrar_csv(cwd)
    if not archivo_excel:
        print(f"❌ No se encontró ningún *_runinfo_merged.csv bajo: {directorio_csv}")
        return
    print(f"📊 Metadata: {archivo_excel}")

    carpeta_destino_grupos = os.path.join(ruta_base, "muestras_clasificadas")

    try:
        df_crudo = pd.read_csv(archivo_excel)

        col_run = 'Run' if 'Run' in df_crudo.columns else 'Run_ID'
        if col_run not in df_crudo.columns:
            print("❌ No se encontró la columna identificadora de la muestra ('Run' o 'Run_ID').")
            return

        col_diag = None
        for col in ['disease_status', 'disease_state', 'disease', 'Disease', 'Diagnosis',
                    'diagnosis', 'disease state', 'isolate', 'subject']:
            if col in df_crudo.columns and not df_crudo[col].isna().all():
                col_diag = col
                break

        if not col_diag:
            print("❌ No se encontró una columna válida de diagnóstico con datos en el archivo.")
            print(f"Columnas disponibles: {list(df_crudo.columns)}")
            return

        print(f"🏷️  Columna de clasificación: '{col_diag}'\n")

        df = df_crudo[[col_run, col_diag]].copy()
        df = df.rename(columns={col_run: 'Run_ID', col_diag: 'Diagnosis'})

    except Exception as e:
        print(f"❌ Error al leer {archivo_excel}. Detalle: {e}")
        return

    df = df.dropna(subset=['Run_ID', 'Diagnosis'])
    df['Diagnosis'] = (df['Diagnosis'].astype(str)
                       .str.replace('Attribute "', '', regex=False)
                       .str.replace('"', '', regex=False))

    os.makedirs(carpeta_destino_grupos, exist_ok=True)

    archivos_movidos = 0
    archivos_no_encontrados = 0

    for _, row in df.iterrows():
        run_id = str(row['Run_ID']).strip()
        diagnosis = str(row['Diagnosis']).strip()

        nombre_carpeta = re.sub(r'[^a-z0-9]+', '_', re.sub(r"'", '', diagnosis.lower())).strip('_')
        ruta_grupo = os.path.join(carpeta_destino_grupos, nombre_carpeta)
        os.makedirs(ruta_grupo, exist_ok=True)

        nombre_archivo = f"{run_id}.fastq.gz"
        ruta_origen = os.path.join(ruta_base, nombre_archivo)
        ruta_destino = os.path.join(ruta_grupo, nombre_archivo)

        if os.path.exists(ruta_origen):
            shutil.move(ruta_origen, ruta_destino)
            print(f"✅ Movido: {nombre_archivo} -> muestras_clasificadas/{nombre_carpeta}/")
            archivos_movidos += 1
        else:
            print(f"⚠️  No encontrado: {nombre_archivo}")
            archivos_no_encontrados += 1

    print("\n--- 📊 Resumen ---")
    print(f"Archivos organizados con éxito: {archivos_movidos}")
    print(f"Archivos no encontrados: {archivos_no_encontrados}")
    print("¡Proceso terminado! 🎉")

if __name__ == "__main__":
    clasificar_fastq()
