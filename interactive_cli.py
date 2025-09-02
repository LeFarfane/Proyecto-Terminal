#!/usr/bin/env python3
"""Interactive assistant for PubMed search parameter collection."""
import os
import json
import importlib.util
import getpass
import re

def prompt_list(message: str, default: str, sep: str) -> list:
    while True:
        value = input(f"{message} [{default}]: ").strip()
        if not value:
            value = default
        items = [item.strip() for item in value.split(sep) if item.strip()]
        if items:
            return items
        print("Please enter at least one item.")

def prompt_choice(message: str, choices: list, default: str) -> str:
    """Prompt user to choose from a list of options (numeric or text)."""
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

def main():
    terms = prompt_list(
        "Search terms (semicolon-separated)",
        "microRNA; celiac disease; inflammatory bowel disease",
        ";",
    )

    operator = prompt_choice(
        "Operator between terms",
        ["AND", "OR", "NOT"],
        "AND",
    )

    sort = prompt_choice(
        "Result order",
        ["relevance", "pubdate"],
        "relevance",
    )

    pub_types_input = input("Publication types (comma-separated) []: ").strip()
    if pub_types_input:
        pub_types = [pt.strip() for pt in pub_types_input.split(",") if pt.strip()]
    else:
        pub_types = []
    pub_types_pubmed = [f"{pt}[Publication Type]" for pt in pub_types]

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

    output_format = prompt_choice(
        "Output format", ["csv", "jsonl", "both"], "both"
    )

    append = prompt_yes_no("Append to existing files", "y")

    output_base = input("Base name for output files [papers]: ").strip() or "papers"

    summary = {
        "terms": terms,
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
    }
    print(json.dumps(summary, indent=2))

    terms_query = " {} ".format(operator).join([f"(\"{t}\")" for t in terms])
    if pub_types_pubmed:
        final_query = f"{terms_query} AND ({' OR '.join(pub_types_pubmed)})"
    else:
        final_query = terms_query
    print("\nRunning search with query:")
    print(final_query)

    # Invoke the PubMed query script with the collected parameters
    spec = importlib.util.spec_from_file_location(
        "PubMed_API_0_1", "PubMed_API_0.1.py"
    )
    pubmed_api = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pubmed_api)
    pubmed_api.run_query(
        terms=terms,
        operator=operator,
        sort=sort,
        pub_types=pub_types,
        pub_types_pubmed=pub_types_pubmed,
        batch_size=batch_size,
        max_results=max_results,
        api_key=api_key,
        format=output_format,
        append=append,
        output_base=output_base,
    )

if __name__ == "__main__":
    main()
