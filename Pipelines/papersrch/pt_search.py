#!/usr/bin/env python3
import argparse
import json
import os
import re
from collections import defaultdict
from pathlib import Path
from typing import List, Dict, Any

import numpy as np
import pandas as pd
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity
import matplotlib.pyplot as plt
import networkx as nx

try:  # Allow running as module or script
    from .text_utils import (
        normalize_text,
        tokenize,
        expand_query_terms,
        DOMAIN_KEYWORDS,
    )
except ImportError:  # pragma: no cover - fallback for direct execution
    from text_utils import (
        normalize_text,
        tokenize,
        expand_query_terms,
        DOMAIN_KEYWORDS,
    )

BASE_DIR = Path(__file__).resolve().parents[2]


# ----------------------------- Data loading -----------------------------

def load_corpus(
    jsonl_path: str | Path | None = None,
    csv_path: str | Path | None = None,
) -> pd.DataFrame:
    jsonl_path = Path(jsonl_path) if jsonl_path else BASE_DIR / "data" / "papers" / "papers.jsonl"
    csv_path = Path(csv_path) if csv_path else BASE_DIR / "data" / "papers" / "papers.csv"
    records: List[Dict[str, Any]] = []
    if jsonl_path.exists():
        with jsonl_path.open("r", encoding="utf-8") as f:
            for line in f:
                try:
                    records.append(json.loads(line))
                except json.JSONDecodeError:
                    continue
    elif csv_path.exists():
        print("Warning: papers.jsonl not found, using CSV")
        df_csv = pd.read_csv(csv_path)
        records = df_csv.to_dict(orient="records")
    else:
        raise FileNotFoundError("No corpus found")
    if not records:
        raise ValueError("Corpus is empty")
    df = pd.DataFrame(records).drop_duplicates("PMID", keep="last")
    df.reset_index(drop=True, inplace=True)
    df["Year"] = pd.to_numeric(df.get("Year"), errors="coerce").astype("Int64")
    df["Title"] = df["Title"].apply(normalize_text)
    df["Abstract"] = df["Abstract"].apply(normalize_text)
    df["Authors"] = df["Authors"].fillna("").apply(normalize_text)
    df["Authors_list"] = df["Authors"].apply(lambda x: [a.strip() for a in x.split(";") if a.strip()])
    df["Journal"] = df["Journal"].fillna("").apply(normalize_text)
    df["DOI"] = df.get("DOI", "").fillna("").apply(normalize_text)
    df["citation_apa"] = df.get("citation_apa", "").fillna("")
    df["Abstract_len"] = df["Abstract"].apply(len)
    if "PublicationTypes" in df.columns:
        df["PublicationTypes_norm"] = df["PublicationTypes"].apply(
            lambda x: x
            if isinstance(x, list)
            else [t.strip() for t in str(x).split(";") if t.strip()] if isinstance(x, str) and x else []
        )
    else:
        df["PublicationTypes_norm"] = [[] for _ in range(len(df))]
    return df


# ---------------------------- Index building ----------------------------

class SearchEngine:
    def __init__(self, df: pd.DataFrame) -> None:
        self.df = df
        self.vectorizer_ti = TfidfVectorizer(ngram_range=(1, 2))
        self.vectorizer_ab = TfidfVectorizer(ngram_range=(1, 2))
        self.tfidf_ti = self.vectorizer_ti.fit_transform(df["Title"].fillna(""))
        self.tfidf_ab = self.vectorizer_ab.fit_transform(df["Abstract"].fillna(""))
        self.inv_ti: Dict[str, set] = defaultdict(set)
        self.inv_ab: Dict[str, set] = defaultdict(set)
        self.sentences: Dict[str, List[List[str]]] = {}
        self.id_to_pos: Dict[str, int] = {}
        for idx, row in df.iterrows():
            pmid = row["PMID"]
            self.id_to_pos[pmid] = idx
            for tok in tokenize(row["Title"]):
                self.inv_ti[tok].add(pmid)
            for tok in tokenize(row["Abstract"]):
                self.inv_ab[tok].add(pmid)
            self.sentences[pmid] = [tokenize(s) for s in re.split(r"[\.\!?]", row["Abstract"]) if s.strip()]

    # -------- candidate docs --------
    def _candidate_pmids(self, query_groups: List[List[str]], op: str, fields: str) -> set:
        if not query_groups:
            return set(self.df["PMID"])
        def index_for(field: str):
            if field == "ti":
                return self.inv_ti
            if field == "ab":
                return self.inv_ab
            comb = defaultdict(set)
            for k, v in self.inv_ti.items():
                comb[k].update(v)
            for k, v in self.inv_ab.items():
                comb[k].update(v)
            return comb
        index = index_for(fields)
        sets = []
        for group in query_groups:
            group_set = set()
            for term in group:
                group_set.update(index.get(term, set()))
            sets.append(group_set)
        all_pmids = set(self.df["PMID"])
        if op == "AND":
            cand = all_pmids
            for s in sets:
                cand &= s
            return cand
        if op == "OR":
            cand = set()
            for s in sets:
                cand.update(s)
            return cand
        if op == "NOT":
            union = set()
            for s in sets:
                union.update(s)
            return all_pmids - union
        return all_pmids

    # -------- filters --------
    def _apply_filters(self, pmids: set, year_min: int, year_max: int, journal_inc: List[str],
                        journal_exc: List[str], author: str, has_doi: bool,
                        exclude_terms: List[str], fields: str) -> set:
        df = self.df.set_index("PMID")
        res = []
        for pmid in pmids:
            row = df.loc[pmid]
            year = row["Year"] if pd.notna(row["Year"]) else None
            if year is not None and year < year_min:
                continue
            if year_max is not None and year is not None and year > year_max:
                continue
            journal = row["Journal"].lower()
            if journal_inc and not any(j.lower() in journal for j in journal_inc):
                continue
            if journal_exc and any(j.lower() in journal for j in journal_exc):
                continue
            if author and author.lower() not in row["Authors"].lower():
                continue
            if has_doi and not row["DOI"]:
                continue
            text_fields = []
            if fields in ("ti", "tiab"):
                text_fields.append(row["Title"])
            if fields in ("ab", "tiab"):
                text_fields.append(row["Abstract"])
            content = " ".join(text_fields).lower()
            if any(t.lower() in content for t in exclude_terms):
                continue
            res.append(pmid)
        return set(res)

    # -------- boosts --------
    def _proximity(self, pmid: str, terms: List[str]) -> bool:
        sentences = self.sentences.get(pmid, [])
        for sent in sentences:
            if sum(1 for t in terms if t in sent) >= 2:
                return True
        return False

    def _domain_boost(self, abstract: str) -> float:
        text = abstract.lower()
        count = sum(1 for kw in DOMAIN_KEYWORDS if kw in text)
        return min(0.2, 0.1 * count)

    # -------- search --------
    def search(self, query: List[str], op: str = "AND", fields: str = "tiab", year_min: int = 2020,
               year_max: int | None = None, journal_inc: List[str] | None = None,
               journal_exc: List[str] | None = None, author: str | None = None,
               has_doi: bool = False, exclude_terms: List[str] | None = None,
               k: int = 30) -> List[Dict[str, Any]]:
        journal_inc = journal_inc or []
        journal_exc = journal_exc or []
        exclude_terms = exclude_terms or []
        groups = expand_query_terms(query)
        pmids = self._candidate_pmids(groups, op, fields)
        pmids = self._apply_filters(pmids, year_min, year_max, journal_inc, journal_exc,
                                    author or "", has_doi, exclude_terms, fields)
        if not pmids:
            return []
        df = self.df[self.df["PMID"].isin(pmids)].copy()
        positions = [self.id_to_pos[p] for p in df["PMID"]]
        q_text = " ".join([t for g in groups for t in g])
        cos_ti = np.zeros(len(df))
        cos_ab = np.zeros(len(df))
        if fields in ("ti", "tiab"):
            q_ti = self.vectorizer_ti.transform([q_text])
            cos_ti = cosine_similarity(q_ti, self.tfidf_ti[positions]).flatten()
        if fields in ("ab", "tiab"):
            q_ab = self.vectorizer_ab.transform([q_text])
            cos_ab = cosine_similarity(q_ab, self.tfidf_ab[positions]).flatten()
        base_scores = 1.2 * cos_ti + 0.8 * cos_ab
        flat_terms = [t for g in groups for t in g]
        results = []
        df = df.reset_index(drop=True)
        for i, row in df.iterrows():
            pmid = row["PMID"]
            score = base_scores[i]
            expl = []
            ti_tokens = tokenize(row["Title"])
            ab_tokens = tokenize(row["Abstract"])
            matched = sorted({t for t in flat_terms if t in ti_tokens or t in ab_tokens})
            if any(t in ti_tokens for t in flat_terms):
                score += 0.2
                expl.append("title")
            year = row["Year"] if pd.notna(row["Year"]) else None
            if year is not None:
                if year <= 2020:
                    rec = 0.0
                elif year >= 2025:
                    rec = 0.2
                else:
                    rec = ((year - 2020) / 5) * 0.2
                score += rec
                if rec:
                    expl.append(f"recency+{rec:.2f}")
            pub_types = (
                row.get("PublicationTypes_norm")
                or row.get("PublicationTypes")
                or row.get("PublicationType")
                or row.get("Publication Types")
            )
            if isinstance(pub_types, str):
                pub_types = [t.strip() for t in re.split(r";|,", pub_types) if t.strip()]
            if isinstance(pub_types, list) and any(
                re.search(r"review|meta-analysis", pt, re.I) for pt in pub_types
            ):
                score += 0.2
                expl.append("review")
            dom = self._domain_boost(row["Abstract"])
            if dom:
                score += dom
                expl.append(f"domain+{dom:.2f}")
            if self._proximity(pmid, flat_terms):
                score += 0.15
                expl.append("proximity")
            results.append({
                "pmid": pmid,
                "title": row["Title"],
                "abstract": row["Abstract"],
                "journal": row["Journal"],
                "year": row["Year"],
                "doi": row["DOI"],
                "citation_apa": row["citation_apa"],
                "score": float(score),
                "cos_title": float(cos_ti[i]) if len(cos_ti) else 0.0,
                "cos_abstract": float(cos_ab[i]) if len(cos_ab) else 0.0,
                "matched_terms": matched,
                "explanation": expl,
                "abstract_len": row["Abstract_len"],
                "has_doi": 1 if row["DOI"] else 0,
            })
        results.sort(key=lambda r: (-r["score"], -(r["year"] or 0), -r["abstract_len"], -r["has_doi"]))
        return results[:k]

    # -------- analytics --------
    def facets(self, by: str) -> pd.DataFrame:
        if by == "journal":
            series = self.df["Journal"]
        elif by == "year":
            series = self.df["Year"]
            series = series[series >= 2020]
        elif by == "author":
            series = self.df["Authors_list"].explode().apply(lambda x: x.split()[-1] if isinstance(x, str) else x)
        else:
            raise ValueError("Invalid facet")
        counts = series.value_counts().head(20)
        return pd.DataFrame({by: counts.index, "count": counts.values, "pct": counts.values / counts.sum() * 100})

    def yearly_counts(self, pmids: List[str] | None = None) -> pd.Series:
        df = self.df
        if pmids is not None:
            df = df[df["PMID"].isin(pmids)]
        return df[df["Year"] >= 2020]["Year"].value_counts().sort_index()

    def coauthor_network(self, pmids: List[str]):
        df = self.df[self.df["PMID"].isin(pmids)]
        G = nx.Graph()
        for authors in df["Authors_list"]:
            for i, a1 in enumerate(authors):
                for a2 in authors[i + 1:]:
                    if G.has_edge(a1, a2):
                        G[a1][a2]["weight"] += 1
                    else:
                        G.add_edge(a1, a2, weight=1)
        bet = nx.betweenness_centrality(G)
        metrics = pd.DataFrame([{ "author": n, "degree": G.degree(n), "betweenness": bet.get(n, 0)} for n in G.nodes])
        return G, metrics

# ------------------------------- CLI ----------------------------------

def parse_terms(value: str) -> List[str]:
    if not value:
        return []
    parts = [p.strip().strip('"') for p in value.split(';') if p.strip()]
    return parts

def add_search_args(p: argparse.ArgumentParser) -> None:
    p.add_argument('--query', type=str, default='', help='terms separated by ;')
    p.add_argument('--op', choices=['AND', 'OR', 'NOT'], default='AND')
    p.add_argument('--fields', choices=['ti', 'ab', 'tiab'], default='tiab')
    p.add_argument('--year-min', type=int, default=2020)
    p.add_argument('--year-max', type=int)
    p.add_argument('--journal-include', type=str, default='')
    p.add_argument('--journal-exclude', type=str, default='')
    p.add_argument('--author', type=str, default='')
    p.add_argument('--has-doi', action='store_true')
    p.add_argument('--exclude', type=str, default='')
    p.add_argument('--k', type=int, default=30)
    p.add_argument('--export-base', type=str)

def run_search(engine: SearchEngine, args) -> None:
    terms = parse_terms(args.query)
    journal_inc = parse_terms(args.journal_include)
    journal_exc = parse_terms(args.journal_exclude)
    exclude = parse_terms(args.exclude)
    results = engine.search(terms, op=args.op, fields=args.fields,
                            year_min=args.year_min, year_max=args.year_max,
                            journal_inc=journal_inc, journal_exc=journal_exc,
                            author=args.author, has_doi=args.has_doi,
                            exclude_terms=exclude, k=args.k)
    if not results:
        print('No results')
        return
    print(f"Showing {len(results)} results:\n")
    header = f"{'rank':<4} {'score':<6} year journal pmid doi title explain"
    print(header)
    for i, r in enumerate(results, 1):
        title = (r['title'][:117] + '...') if len(r['title']) > 120 else r['title']
        explain = ','.join(r['explanation'])
        print(f"{i:<4} {r['score']:.3f} {r['year']} {r['journal'][:20]} {r['pmid']} {r['doi'][:15]} {title} {explain}")
    if args.export_base:
        base = args.export_base
        os.makedirs('outputs', exist_ok=True)
        df = pd.DataFrame(results)
        df.to_csv(f'outputs/{base}_results.csv', index=False)
        with open(f'outputs/{base}_results.jsonl', 'w', encoding='utf-8') as f:
            for rec in results:
                f.write(json.dumps(rec, ensure_ascii=False) + '\n')
        with open(f'outputs/{base}_citations.txt', 'w', encoding='utf-8') as f:
            for rec in results:
                f.write(rec['citation_apa'] + '\n')
        print(f"Exported to outputs/{base}_results.csv/jsonl and citations.txt")


def run_facets(engine: SearchEngine, args) -> None:
    df = engine.facets(args.by)
    print(df.to_string(index=False))


def run_plot_yearly(engine: SearchEngine, args) -> None:
    pmids = None
    if args.filtered and args.query:
        terms = parse_terms(args.query)
        pmids = [r['pmid'] for r in engine.search(terms, op=args.op, fields=args.fields,
                                                  year_min=args.year_min, year_max=args.year_max,
                                                  journal_inc=parse_terms(args.journal_include),
                                                  journal_exc=parse_terms(args.journal_exclude),
                                                  author=args.author, has_doi=args.has_doi,
                                                  exclude_terms=parse_terms(args.exclude),
                                                  k=1000)]
    counts = engine.yearly_counts(pmids)
    os.makedirs('outputs', exist_ok=True)
    plt.figure()
    counts.sort_index().plot(kind='bar')
    plt.ylabel('count')
    plt.title('Yearly publications')
    plt.tight_layout()
    plt.savefig('outputs/yearly_counts.png')
    print('Saved outputs/yearly_counts.png')


def run_coauthors(engine: SearchEngine, args) -> None:
    terms = parse_terms(args.query)
    pmids = [r['pmid'] for r in engine.search(terms, op=args.op, fields=args.fields,
                                              year_min=args.year_min, year_max=args.year_max,
                                              journal_inc=parse_terms(args.journal_include),
                                              journal_exc=parse_terms(args.journal_exclude),
                                              author=args.author, has_doi=args.has_doi,
                                              exclude_terms=parse_terms(args.exclude),
                                              k=1000)]
    if not pmids:
        print('No data for coauthors')
        return
    G, metrics = engine.coauthor_network(pmids)
    os.makedirs('outputs', exist_ok=True)
    nx.write_gexf(G, 'outputs/coauthors.gexf')
    metrics.to_csv('outputs/coauthors_metrics.csv', index=False)
    print('Saved outputs/coauthors.gexf and coauthors_metrics.csv')


def main():
    parser = argparse.ArgumentParser(description='Local search over PubMed corpus')
    sub = parser.add_subparsers(dest='cmd', required=True)
    p_search = sub.add_parser('search')
    add_search_args(p_search)
    p_search.set_defaults(func=run_search)

    p_facets = sub.add_parser('facets')
    p_facets.add_argument('--by', choices=['journal', 'year', 'author'], required=True)
    p_facets.set_defaults(func=run_facets)

    p_plot = sub.add_parser('plot-yearly')
    add_search_args(p_plot)
    p_plot.add_argument('--filtered', action='store_true')
    p_plot.set_defaults(func=run_plot_yearly)

    p_co = sub.add_parser('coauthors')
    add_search_args(p_co)
    p_co.set_defaults(func=run_coauthors)

    args = parser.parse_args()
    df = load_corpus()
    engine = SearchEngine(df)
    args.func(engine, args)

if __name__ == '__main__':
    main()
