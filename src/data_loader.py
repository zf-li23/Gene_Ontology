"""Utilities for retrieving and loading interaction networks and GO annotations."""
from __future__ import annotations

import gzip
import io
import logging
import random
import shutil
from pathlib import Path
from typing import Dict, Iterable, List, Optional

import networkx as nx
import pandas as pd
import numpy as np
import requests

LOGGER = logging.getLogger(__name__)


def download_file(url: str, destination: Path, chunk_size: int = 1 << 20) -> Path:
    """Stream a file from ``url`` to ``destination``. Creates parents when needed."""
    destination.parent.mkdir(parents=True, exist_ok=True)
    LOGGER.info("Downloading %s -> %s", url, destination)
    with requests.get(url, stream=True, timeout=60) as response:
        response.raise_for_status()
        with open(destination, "wb") as fh:
            for chunk in response.iter_content(chunk_size=chunk_size):
                if chunk:
                    fh.write(chunk)
    return destination


def download_string_network(species_id: int, destination: Path, score_cutoff: int = 400) -> Path:
    """Download a STRING PPI edge list (gzipped) for ``species_id`` and filter by score."""
    base = "https://stringdb-static.org/download/protein.links.v12.0"
    url = f"{base}/{species_id}.protein.links.v12.0.txt.gz"
    gz_path = destination.with_suffix(destination.suffix + ".gz")
    download_file(url, gz_path)
    LOGGER.info("Filtering STRING links with combined score >= %s", score_cutoff)
    with gzip.open(gz_path, "rt") as src, open(destination, "w") as dst:
        header = src.readline().strip().split()
        dst.write("gene_a\tgene_b\tweight\n")
        for line in src:
            tokens = line.strip().split()
            if not tokens:
                continue
            score = int(tokens[2])
            if score >= score_cutoff:
                dst.write(f"{tokens[0]}\t{tokens[1]}\t{score / 1000:.3f}\n")
    gz_path.unlink(missing_ok=True)
    return destination


def load_edge_list(path: Path, weighted: bool = True) -> pd.DataFrame:
    """Load an interaction list with columns ``gene_a``, ``gene_b``, optional ``weight``."""
    df = pd.read_csv(path, sep="\t")
    required = {"gene_a", "gene_b"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Edge list missing columns: {missing}")
    if weighted and "weight" not in df.columns:
        LOGGER.warning("Weight column absent; injecting default weight = 1.0")
        df["weight"] = 1.0
    if not weighted:
        df["weight"] = 1.0
    return df[["gene_a", "gene_b", "weight"]]


def convert_scrin_csv_to_tsv(csv_path: Path, tsv_path: Optional[Path] = None,
                             pval_col: str = "pvalue", qval_col: str = "qvalue_BH",
                             enrich_col: str = "enrichment_ratio") -> Path:
    """Convert scrin CSV (gene_A,gene_B,...) to TSV with a more robust weight.

    Weight strategy (bounded, positive):
    w = log1p(enrichment_ratio_if_available) * (-log10(qvalue_BH_if_available else pvalue))
    Values are clipped to reasonable ranges to avoid exploding weights.
    """
    if tsv_path is None:
        tsv_path = csv_path.with_name(f"{csv_path.stem}_scrin_edges.tsv")
    df = pd.read_csv(csv_path)
    required = {"gene_A", "gene_B", pval_col}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"scrin CSV missing columns: {missing}")

    small = 1e-300
    # choose significance source: prefer qvalue if present
    if qval_col in df.columns:
        sig = df[qval_col].clip(lower=small)
    else:
        sig = df[pval_col].clip(lower=small)
    log_sig = -np.log10(sig)

    # effect size: enrichment_ratio when available; otherwise 1
    if enrich_col in df.columns:
        effect = np.log1p(df[enrich_col].clip(lower=0))
    else:
        effect = 1.0

    weight = (log_sig * effect).astype(float)
    # clip to prevent extreme outliers
    weight = np.clip(weight, 0, np.quantile(weight, 0.999))

    df_out = pd.DataFrame({
        "gene_a": df["gene_A"],
        "gene_b": df["gene_B"],
        "weight": weight,
    })
    tsv_path.parent.mkdir(parents=True, exist_ok=True)
    df_out.to_csv(tsv_path, sep="\t", index=False)
    LOGGER.info("Converted scrin CSV -> %s (%s edges)", tsv_path, len(df_out))
    return tsv_path


def load_go_annotations(path: Path) -> Dict[str, List[str]]:
    """Parse GO annotations as mapping gene -> list of GO terms."""
    df = pd.read_csv(path, sep="\t")
    required = {"gene_id", "go_term"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"GO file missing columns: {missing}")
    grouped = df.groupby("gene_id")["go_term"].apply(list)
    return grouped.to_dict()


def convert_string_links_to_tsv(input_path: Path, out_path: Optional[Path] = None, score_col: int = 2) -> Path:
    """Convert STRING/COG-style links file to a TSV with columns gene_a,gene_b,weight.

    The input is expected to be whitespace-separated with at least three columns: proteinA proteinB combined_score
    Weight is normalized to [0,1] by dividing combined_score by 1000.
    """
    if out_path is None:
        out_path = input_path.with_suffix(".tsv")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(input_path, "r") as src, open(out_path, "w") as dst:
        # try to detect and skip header lines that are not data
        first = src.readline()
        tokens = first.strip().split()
        if len(tokens) < 3 or not tokens[2].replace('.', '').isdigit():
            # assume first line is header; write header and continue
            pass
        else:
            # first line is data; reset file pointer
            src.seek(0)
        dst.write("gene_a\tgene_b\tweight\n")
        for line in src:
            if not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) <= score_col:
                continue
            a, b = parts[0], parts[1]
            try:
                score = float(parts[score_col])
            except Exception:
                # skip malformed score
                continue
            weight = max(0.0, min(1.0, score / 1000.0))
            dst.write(f"{a}\t{b}\t{weight:.3f}\n")
    return out_path


def generate_demo_network(edge_path: Path, go_path: Path, num_nodes: int = 30, seed: int = 7) -> None:
    """Write a quick synthetic network and GO table for offline experimentation."""
    rng = random.Random(seed)
    graph = nx.barabasi_albert_graph(num_nodes, 2, seed=seed)
    for u, v in graph.edges():
        graph.edges[u, v]["weight"] = round(rng.uniform(0.4, 1.0), 3)
    genes = [f"G{i:03d}" for i in range(num_nodes)]
    nx.relabel_nodes(graph, dict(zip(graph.nodes(), genes)), copy=False)
    edge_path.parent.mkdir(parents=True, exist_ok=True)
    with open(edge_path, "w") as fh:
        fh.write("gene_a\tgene_b\tweight\n")
        for u, v, data in graph.edges(data=True):
            fh.write(f"{u}\t{v}\t{data['weight']:.3f}\n")
    terms = ["GO:0008152", "GO:0009987", "GO:0006412"]
    go_path.parent.mkdir(parents=True, exist_ok=True)
    with open(go_path, "w") as fh:
        fh.write("gene_id\tgo_term\tdescription\n")
        for gene in genes:
            annotated = rng.sample(terms, k=rng.randint(1, len(terms)))
            for term in annotated:
                fh.write(f"{gene}\t{term}\tDemo term\n")


def ensure_inputs(edge_path: Path, go_path: Path) -> None:
    """Make sure both network and GO files exist, generating demo data when necessary."""
    if edge_path.exists() and go_path.exists():
        return
    LOGGER.warning("Input files missing. Generating demo network to keep the pipeline runnable.")
    generate_demo_network(edge_path=edge_path, go_path=go_path)
