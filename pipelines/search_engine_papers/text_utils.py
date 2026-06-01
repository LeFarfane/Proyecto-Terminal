import re
import unicodedata
from typing import List

GREEK_MAP = {
    'β': 'beta',
    'α': 'alpha',
    'κ': 'kappa',
    'γ': 'gamma',
    'δ': 'delta',
}

DOMAIN_KEYWORDS = [
    "nf-kb", "tgf-beta", "il-6", "tnf", "t cell", "epithelial barrier",
    "autophagy", "mucosa", "tight junction", "smad2", "smad3"
]


def normalize_text(text: str) -> str:
    if not isinstance(text, str):
        return ""
    text = unicodedata.normalize('NFKC', text)
    for k, v in GREEK_MAP.items():
        text = text.replace(k, v)
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def tokenize(text: str) -> List[str]:
    text = normalize_text(text).lower()
    return re.findall(r"\b\w+\b", text)


def expand_query_terms(terms: List[str]) -> List[List[str]]:
    expanded = []
    for term in terms:
        norm = normalize_text(term).lower()
        variants = {norm}
        if any(tok in norm for tok in ["microrna", "mirna", "mir-"]):
            variants.update({"microrna", "mirna", "mir"})
        if norm in {"ibd", "inflammatory bowel disease"}:
            variants.update({"ibd", "inflammatory bowel disease"})
        if re.search(r"celiac|coeliac", norm):
            variants.update({"celiac", "coeliac"})
        expanded.append(sorted(variants))
    return expanded
