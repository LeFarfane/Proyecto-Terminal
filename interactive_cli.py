#!/usr/bin/env python3
"""Interactive assistant for PubMed search parameter collection."""
import os
import json
import importlib.util

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
    choices_lower = [c.lower() for c in choices]
    while True:
        value = input(f"{message} [{default}]: ").strip()
        if not value:
            value = default
        if value.lower() in choices_lower:
            return choices[choices_lower.index(value.lower())]
        print(f"Invalid choice. Valid options: {', '.join(choices)}")

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
    pubmed_pub_types = [f"{pt}[Publication Type]" for pt in pub_types]
    # pubmed_pub_types is prepared for potential use

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
            api_key = input("API key: ").strip()
            api_key_provided = bool(api_key)

    output_format = prompt_choice(
        "Output format", ["csv", "jsonl", "both"], "both"
    )

    append = prompt_yes_no("Append to existing files", "y")

    summary = {
        "terms": terms,
        "operator": operator,
        "sort": sort,
        "pub_types": pub_types,
        "batch_size": batch_size,
        "max_results": max_results,
        "api_key_provided": api_key_provided,
        "format": output_format,
        "append": append,
    }
    print(json.dumps(summary, indent=2))

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
        batch_size=batch_size,
        max_results=max_results,
        api_key=api_key,
        format=output_format,
        append=append,
    )

if __name__ == "__main__":
    main()
