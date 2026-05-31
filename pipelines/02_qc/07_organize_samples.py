import os
import shutil
import pandas as pd

def clasificar_fastq():
    print("🧬 Organizador de Muestras SRA (Modo Ruta Única) 🧬\n")
    
    # 1. Pedir la ruta base única
    ruta_base = input("📂 Introduce la ruta de la carpeta de trabajo (ej. ~/samples o . para la carpeta actual): ").strip()
    if not ruta_base:
        ruta_base = "."
        
    ruta_base = os.path.expanduser(ruta_base)
    
    # 2. Pedir el nombre del Excel
    nombre_excel = input("📊 Introduce el nombre del archivo Excel (ej. metadatos.csv): ").strip()
    archivo_excel = os.path.join(ruta_base, nombre_excel)
    
    carpeta_destino_grupos = os.path.join(ruta_base, "muestras_clasificadas")
    
    print(f"\nBuscando archivos en: {ruta_base} 🚀")
    
    try:
        # Leemos el archivo completo primero
        df_crudo = pd.read_csv(archivo_excel)
        
        # 1. Identificar la columna de Run (puede llamarse 'Run' o 'Run_ID')
        col_run = 'Run' if 'Run' in df_crudo.columns else 'Run_ID'
        if col_run not in df_crudo.columns:
            print("❌ No se encontró la columna identificadora de la muestra ('Run' o 'Run_ID').")
            return

        # 2. Identificar la columna de diagnóstico correcta (que tiene los datos)
        # Priorizamos 'disease' en minúscula porque en tu CSV 'Disease' en mayúscula está vacía.
        col_diag = None
        for col in ['disease_state', 'disease', 'Disease', 'Diagnosis', 'diagnosis', 'disease state', 'isolate']:
            # Verificamos que la columna exista y que no esté compuesta completamente por valores nulos
            if col in df_crudo.columns and not df_crudo[col].isna().all():
                col_diag = col
                break
                
        if not col_diag:
            print("❌ No se encontró una columna válida de diagnóstico con datos en el archivo.")
            print(f"Columnas disponibles: {list(df_crudo.columns)}")
            return
        
        # Filtramos solo las columnas que necesitamos y las estandarizamos
        df = df_crudo[[col_run, col_diag]].copy()
        df = df.rename(columns={col_run: 'Run_ID', col_diag: 'Diagnosis'})

    except Exception as e:
        print(f"❌ Error al leer el archivo en {archivo_excel}.\nVerifica que el nombre sea correcto. Detalle: {e}")
        return

    # Limpiar datos nulos
    df = df.dropna(subset=['Run_ID', 'Diagnosis'])
    
    # Limpiar el formato 'Attribute "Texto"' para que las carpetas tengan nombres limpios
    df['Diagnosis'] = df['Diagnosis'].astype(str).str.replace('Attribute "', '', regex=False).str.replace('"', '', regex=False)
    
    # Crear la carpeta contenedora
    os.makedirs(carpeta_destino_grupos, exist_ok=True)

    archivos_movidos = 0
    archivos_no_encontrados = 0

    # Iterar sobre las muestras
    for index, row in df.iterrows():
        run_id = str(row['Run_ID']).strip()
        diagnosis = str(row['Diagnosis']).strip()

        # Limpiar nombre del grupo y crear su ruta
        nombre_carpeta = diagnosis.replace(" ", "_").replace("/", "-")
        ruta_grupo = os.path.join(carpeta_destino_grupos, nombre_carpeta)
        os.makedirs(ruta_grupo, exist_ok=True)

        # Mover el archivo
        nombre_archivo = f"{run_id}.fastq.gz"
        ruta_origen = os.path.join(ruta_base, nombre_archivo)
        ruta_destino = os.path.join(ruta_grupo, nombre_archivo)

        if os.path.exists(ruta_origen):
            shutil.move(ruta_origen, ruta_destino)
            print(f"✅ Movido: {nombre_archivo} -> muestras_clasificadas/{nombre_carpeta}/")
            archivos_movidos += 1
        else:
            print(f"⚠️ Advertencia: No se encontró {nombre_archivo} en la carpeta de trabajo.")
            archivos_no_encontrados += 1

    print("\n--- 📊 Resumen ---")
    print(f"Archivos organizados con éxito: {archivos_movidos}")
    print(f"Archivos no encontrados: {archivos_no_encontrados}")
    print("¡Proceso terminado! 🎉")

if __name__ == "__main__":
    clasificar_fastq()