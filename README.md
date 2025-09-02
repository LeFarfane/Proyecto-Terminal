# Proyecto-Terminal

## Última búsqueda
- Término: ("IBS") AND (Systematic Review[Publication Type])
- Artículos guardados: 200

Generado el 2025-09-02 13:58:18 UTC con versión 0.2.

## Environment Variables

Set the `NCBI_API_KEY` environment variable to use an NCBI API key for PubMed requests:

```bash
export NCBI_API_KEY="your_api_key_here"
```

This increases the request rate limits when running `PubMed_API_0.1.py`.

## CLI de búsqueda local

Instala dependencias y ejecuta los comandos desde la raíz del repositorio.

### Ejemplos

```bash
# Búsqueda principal (AND por default), título+abstract, año >= 2020
python pt_search.py search --query "microRNA; IBD" --year-min 2020 --k 30

# Operador OR y filtro por journal
python pt_search.py search --query "microRNA; celiac" --op OR --journal-include "Gut; Scientific Reports"

# Excluir términos
python pt_search.py search --query "microRNA; fibrosis" --exclude "exosome"

# Solo en título y con exportación
python pt_search.py search --query "TGF-beta" --fields ti --export-base fibrosis_tgf

# Facetas por journal
python pt_search.py facets --by journal

# Gráfica por año del corpus completo
python pt_search.py plot-yearly

# Red de coautores sobre una búsqueda
python pt_search.py coauthors --query "microRNA; IBD" --year-min 2020
```
