#!/usr/bin/env python3
"""Interactive assistant for PubMed search parameter collection."""
# Script para recopilar parámetros de búsqueda y ejecutar consultas a PubMed.
import os
import json
import importlib.util
import getpass
import re
import inspect

def prompt_list(message: str, default: str, sep: str) -> list:
    """Obtiene una lista a partir de una cadena separada por un delimitador."""
    while True:
        value = input(f"{message} [{default}]: ").strip()
        if not value:
            value = default
        items = [item.strip() for item in value.split(sep) if item.strip()]
        if items:
            return items
        print("Please enter at least one item.")

def prompt_choice(message: str, choices: list, default: str) -> str:
    """Muestra un menú para elegir una opción numérica o por texto."""
    choices_lower = [c.lower() for c in choices]
    while True:
        print(f"{message} [default={default}]:")
        for idx, choice in enumerate(choices, 1):
            print(f"{idx}) {choice}")
        value = input("Choose: ").strip()
        if not value:
            return default
        if value.isdigit():
            idx = int(value) - 1
            if 0 <= idx < len(choices):
                return choices[idx]
        else:
            if value.lower() in choices_lower:
                return choices[choices_lower.index(value.lower())]
        print(f"Invalid choice. Choose number or one of: {', '.join(choices)}")

def prompt_int(message: str, default: int, min_value: int = None, max_value: int = None) -> int:
    """Pide un entero validando rangos mínimo y máximo."""
    while True:
        value = input(f"{message} [{default}]: ").strip()
        if not value:
            return default
        if not value.isdigit():
            print("Please enter a valid integer.")
            continue
        num = int(value)
        if min_value is not None and num < min_value:
            print(f"Please enter a value ≥ {min_value}.")
            continue
        if max_value is not None and num > max_value:
            print(f"Please enter a value ≤ {max_value}.")
            continue
        return num

def prompt_yes_no(message: str, default: str) -> bool:
    """Interpreta respuestas afirmativas o negativas del usuario."""
    yes = {"y", "yes", "s", "si", "sí"}
    no = {"n", "no"}
    default_lower = default.lower()
    while True:
        value = input(f"{message} [{default}]: ").strip()
        if not value:
            value = default
        value_lower = value.lower()
        if value_lower in yes:
            return True
        if value_lower in no:
            return False
        print("Please answer y/n or sí/no.")


# Lista curada de tipos de publicación
PUB_TYPE_OPTIONS = [
    "Review",
    "Systematic Review",
    "Meta-Analysis",
    "Clinical Trial",
    "Randomized Controlled Trial",
    "Clinical Trial, Phase I",
    "Clinical Trial, Phase II",
    "Clinical Trial, Phase III",
    "Clinical Trial, Phase IV",
    "Pragmatic Clinical Trial",
    "Multicenter Study",
    "Observational Study",
    "Comparative Study",
    "Evaluation Studies",
    "Validation Study",
    "Case Reports",
    "Guideline",
    "Practice Guideline",
    "Consensus Development Conference",
    "Editorial",
    "Letter",
    "Comment",
    "Journal Article",
    "Other...",
]

# Sinónimos y abreviaturas aceptados
SYNONYMS = {
    "rct": "Randomized Controlled Trial",
    "meta": "Meta-Analysis",
    "sr": "Systematic Review",
    "ct": "Clinical Trial",
    "other": "Other...",
    "other...": "Other...",
}


def _normalize(text: str) -> str:
    """Quita espacios extra y normaliza la cadena."""
    return re.sub(r"\s+", " ", text.strip())


def _prompt_other_types(existing: list) -> list:
    """Pregunta al usuario por tipos adicionales cuando elige 'Other'."""
    while True:
        other_input = input("Otros tipos (separados por coma): ").strip()
        if not other_input:
            print("Por favor ingresa al menos un tipo.")
            continue
        tokens = [
            _normalize(SYNONYMS.get(tok.lower(), tok))
            for tok in other_input.split(",")
            if _normalize(tok)
        ]
        result = []
        for tok in tokens:
            if tok not in existing and tok not in result:
                result.append(tok)
        if result:
            return result
        print("Por favor ingresa al menos un tipo válido.")


def prompt_pub_types() -> list:
    """Gestiona la selección múltiple de tipos de publicación."""
    options_lower = {opt.lower(): opt for opt in PUB_TYPE_OPTIONS}
    while True:
        print(
            "Publication types (elige números, rangos o escribe nombres) [Enter = ninguno]:"
        )
        for idx, opt in enumerate(PUB_TYPE_OPTIONS, 1):
            print(f"{idx}) {opt}")
        selection = input("Selección: ").strip()
        if not selection:
            return []
        tokens = [t.strip() for t in selection.split(",") if t.strip()]
        final = []
        asked_other = False
        valid = True
        for token in tokens:
            token_norm = _normalize(token)
            # Manejo de rangos numéricos como 2-5
            if re.fullmatch(r"\d+-\d+", token_norm):
                start, end = map(int, token_norm.split("-"))
                if start < 1 or end > len(PUB_TYPE_OPTIONS) or start > end:
                    valid = False
                    break
                for idx in range(start, end + 1):
                    choice = PUB_TYPE_OPTIONS[idx - 1]
                    if choice == "Other...":
                        if not asked_other:
                            for other in _prompt_other_types(final):
                                final.append(other)
                            asked_other = True
                    elif choice not in final:
                        final.append(choice)
            # Manejo de números individuales
            elif token_norm.isdigit():
                idx = int(token_norm)
                if idx < 1 or idx > len(PUB_TYPE_OPTIONS):
                    valid = False
                    break
                choice = PUB_TYPE_OPTIONS[idx - 1]
                if choice == "Other...":
                    if not asked_other:
                        for other in _prompt_other_types(final):
                            final.append(other)
                        asked_other = True
                elif choice not in final:
                    final.append(choice)
            else:
                # Manejo de nombres y sinónimos
                mapped = SYNONYMS.get(token_norm.lower(), token_norm)
                if mapped.lower() in options_lower:
                    choice = options_lower[mapped.lower()]
                    if choice == "Other...":
                        if not asked_other:
                            for other in _prompt_other_types(final):
                                final.append(other)
                            asked_other = True
                    elif choice not in final:
                        final.append(choice)
                else:
                    valid = False
                    break
        if valid:
            return final
        print("Entrada inválida. Intenta de nuevo.")


def prompt_systematic_sb() -> bool:
    """Pregunta si se debe añadir el filtro "systematic [sb]"."""
    yes = {"y", "yes", "s", "si", "sí"}
    no = {"n", "no", ""}
    while True:
        val = input(
            '¿Añadir filtro "systematic [sb]" a la búsqueda? [y/n] (n): '
        ).strip().lower()
        if val in yes:
            return True
        if val in no:
            return False
        print("Please answer y/n or sí/no.")

def main():
    """Flujo principal para recolectar parámetros y ejecutar la búsqueda."""
    query_terms = prompt_list(
        "Search terms (semicolon-separated)",
        "microRNA; celiac disease; inflammatory bowel disease",
        ";",
    )

    # Operador lógico entre términos
    operator = prompt_choice(
        "Operator between terms",
        ["AND", "OR", "NOT"],
        "AND",
    )

    # Orden de los resultados
    sort = prompt_choice(
        "Result order",
        ["relevance", "pubdate"],
        "relevance",
    )
    # Selección de tipos de publicación y filtro systematic
    pub_types = prompt_pub_types()
    systematic_sb = prompt_systematic_sb()
    pub_types_pubmed = [f"{pt}[Publication Type]" for pt in pub_types]
    print(
        "\nPublication types seleccionados: "
        + (", ".join(pub_types) if pub_types else "ninguno")
    )
    print(f"systematic [sb]: {'sí' if systematic_sb else 'no'}\n")
    # Parámetros de descarga y límite de resultados
    batch_size = prompt_int(
        "efetch batch_size",
        100,
        min_value=50,
        max_value=200,
    )

    max_results = prompt_int(
        "Result limit (max_results)",
        200,
        min_value=1,
    )
    # Lectura del API key desde el entorno o entrada manual
    api_key = os.getenv("NCBI_API_KEY", "").strip()
    api_key_provided = bool(api_key)
    if not api_key:
        if prompt_yes_no("Do you want to paste an API key now? (y/n)", "n"):
            while True:
                key = getpass.getpass("API key: ").strip()
                if not key:
                    api_key = ""
                    break
                if re.fullmatch(r"[A-Za-z0-9]{20,80}", key):
                    api_key = key
                    break
                print("Invalid API key format. Enter 20-80 alphanumeric characters or leave blank to skip.")
            api_key_provided = bool(api_key)

    # Formato y nombre de los archivos de salida
    output_format = prompt_choice(
        "Output format", ["csv", "jsonl", "both"], "both"
    )

    append = prompt_yes_no("Append to existing files", "y")

    output_base = input("Base name for output files [papers]: ").strip() or "papers"

    # Resumen de parámetros elegido por el usuario
    summary = {
        "query_terms": query_terms,
        "operator": operator,
        "sort": sort,
        "pub_types": pub_types,
        "pub_types_pubmed": pub_types_pubmed,
        "batch_size": batch_size,
        "max_results": max_results,
        "api_key_provided": api_key_provided,
        "format": output_format,
        "append": append,
        "output_base": output_base,
        "systematic_sb": systematic_sb,
    }
    print(json.dumps(summary, indent=2))
    # Construcción del query final
    terms_query = " {} ".format(operator).join([f"(\"{t}\")" for t in query_terms])
    if pub_types_pubmed:
        final_query = f"{terms_query} AND ({' OR '.join(pub_types_pubmed)})"
    else:
        final_query = terms_query
    print("\nRunning search with query:")
    print(final_query)

    # Carga dinámica del módulo de consulta a PubMed
    spec = importlib.util.spec_from_file_location(
        "PubMed_API_0_1", "PubMed_API_0.1.py"
    )
    pubmed_api = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pubmed_api)
    run_signature = inspect.signature(pubmed_api.run_query)
    run_kwargs = {
        "terms": query_terms,
        "operator": operator,
        "sort": sort,
        "pub_types": pub_types,
        "batch_size": batch_size,
        "max_results": max_results,
        "api_key": api_key,
        "format": output_format,
        "append": append,
        "output_base": output_base,
    }
    # Se pasan argumentos opcionales si la función los soporta
    if "pub_types_pubmed" in run_signature.parameters:
        run_kwargs["pub_types_pubmed"] = pub_types_pubmed
    if "systematic_sb" in run_signature.parameters:
        run_kwargs["systematic_sb"] = systematic_sb
    pubmed_api.run_query(**run_kwargs)

if __name__ == "__main__":
    main()
