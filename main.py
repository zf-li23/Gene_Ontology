"""End-to-end entry point for the intelligent bio-computation workflow.

Supports two modes:
- Default (no CLI edge file): runs the original config-driven pipeline with GO enrichment.
- CLI edge file provided: runs the unsupervised overlapping hierarchical detector and prints JSON.
"""
from __future__ import annotations

import json
import logging
import argparse
from collections import defaultdict
from pathlib import Path
from typing import Dict, List

import pandas as pd
import yaml

from src import data_loader
from src.algorithms.hierarchical_louvain import HierarchicalLouvain, HierarchyNode, flatten_hierarchy
from src.algorithms.louvain_baseline import run_louvain
from src import preprocess
from src import enrichment_analysis
from src import evaluation
from src import visualization
from src.overlap_detector import run_pipeline as run_overlap


LOGGER = logging.getLogger("intelligent_bio")


def configure_logging() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(name)s | %(message)s")


def load_config(config_path: Path) -> Dict:
    with open(config_path) as fh:
        return yaml.safe_load(fh)


def partition_to_dict(partition: Dict[str, int], prefix: str) -> Dict[str, List[str]]:
    grouped: Dict[int, List[str]] = defaultdict(list)
    for gene, community_id in partition.items():
        grouped[community_id].append(gene)
    labeled = {f"{prefix}{idx}": sorted(members) for idx, members in grouped.items()}
    return labeled


def hierarchy_to_dict(root: HierarchyNode) -> Dict[str, List[str]]:
    levels = flatten_hierarchy(root)
    labeled: Dict[str, List[str]] = {}
    for depth, communities in levels.items():
        for idx, members in enumerate(communities):
            labeled[f"H{depth}_{idx}"] = sorted(members)
    return labeled


def write_enrichment_table(results: Dict[str, List[Dict]], output_file: Path) -> None:
    rows = []
    for community_id, terms in results.items():
        for term in terms:
            rows.append({"community": community_id, **term})
    if not rows:
        LOGGER.info("No enrichments to write for %s", output_file)
        return
    df = pd.DataFrame(rows)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_file, index=False)


def resolve_dataset_paths(config: Dict, project_root: Path) -> Dict[str, Path]:
    paths_cfg = config.get("paths", {})
    dataset = paths_cfg.get("dataset")
    datasets = paths_cfg.get("datasets")
    if dataset and datasets and dataset in datasets:
        entry = datasets[dataset]
        return {
            "edge_path": project_root / entry["network_edges"],
            "go_path": project_root / entry["go_annotations"],
        }
    # legacy fallback
    return {
        "edge_path": project_root / paths_cfg["network_edges"],
        "go_path": project_root / paths_cfg["go_annotations"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Intelligent bio-computation pipeline")
    parser.add_argument("edge_file", nargs="?", help="Optional: run overlapping hierarchical detector on a custom edge list")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--prune", type=float, default=0.3)
    parser.add_argument("--resolution", type=float, default=1.0)
    args = parser.parse_args()

    if args.edge_file:
        configure_logging()
        output = run_overlap(args.edge_file, threads=args.threads, resolution=args.resolution, prune=args.prune)
        print(json.dumps(output, indent=2))
        return

    configure_logging()
    project_root = Path(__file__).resolve().parent
    config = load_config(project_root / "config.yaml")
    paths = config["paths"]
    dataset_paths = resolve_dataset_paths(config, project_root)
    edge_path = dataset_paths["edge_path"]
    go_path = dataset_paths["go_path"]
    output_dir = project_root / paths["output_dir"]
    figures_dir = project_root / paths["figures_dir"]
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    # Only auto-generate demo data when using the sample dataset; otherwise expect files to exist
    if config.get("paths", {}).get("dataset", "sample") == "sample":
        data_loader.ensure_inputs(edge_path, go_path)
    else:
        if not edge_path.exists() or not go_path.exists():
            raise FileNotFoundError(f"Expected data files not found: {edge_path} and/or {go_path}")
    edges = data_loader.load_edge_list(edge_path, weighted=config["preprocessing"]["weighted"])
    go_annotations = data_loader.load_go_annotations(go_path)

    graph = preprocess.prepare_graph(
        edges,
        keep_giant_component=config["preprocessing"]["keep_giant_component"],
        weighted=config["preprocessing"]["weighted"],
    )

    partition = run_louvain(
        graph,
        resolution=config["algorithms"]["louvain"]["resolution"],
        random_state=config["algorithms"]["louvain"].get("random_state"),
    )
    baseline_comm = partition_to_dict(partition, prefix="L")

    hierarchical = HierarchicalLouvain(
        min_size=config["algorithms"]["hierarchical"]["min_size"],
        max_depth=config["algorithms"]["hierarchical"]["max_depth"],
        modularity_threshold=config["algorithms"]["hierarchical"]["modularity_threshold"],
        resolution=config["algorithms"]["louvain"]["resolution"],
        random_state=config["algorithms"]["louvain"].get("random_state"),
    )
    hierarchy_root = hierarchical.fit(graph)
    hierarchical_comm = hierarchy_to_dict(hierarchy_root)

    background = enrichment_analysis.background_from_strategy(
        config["enrichment"]["background_strategy"],
        graph.nodes(),
        go_annotations,
    )
    baseline_enrich = enrichment_analysis.enrich_all_communities(
        baseline_comm,
        background,
        go_annotations,
        fdr_threshold=config["enrichment"]["fdr_threshold"],
        min_genes=config["enrichment"]["min_genes"],
    )
    hierarchical_enrich = enrichment_analysis.enrich_all_communities(
        hierarchical_comm,
        background,
        go_annotations,
        fdr_threshold=config["enrichment"]["fdr_threshold"],
        min_genes=config["enrichment"]["min_genes"],
    )

    consistency = evaluation.community_consistency(baseline_comm, go_annotations)
    metrics = evaluation.summarize_metrics(baseline_enrich, consistency)
    metrics_path = output_dir / "metrics.json"
    metrics_path.write_text(json.dumps(metrics, indent=2))
    LOGGER.info("Pipeline metrics: %s", metrics)

    write_enrichment_table(baseline_enrich, output_dir / "baseline_enrichment.csv")
    write_enrichment_table(hierarchical_enrich, output_dir / "hierarchical_enrichment.csv")

    visualization.plot_network(graph, partition, figures_dir / "network.png", layout=config["visualization"]["layout"])
    visualization.plot_hierarchy(hierarchy_root, figures_dir / "hierarchy.png")
    visualization.plot_enrichment_bar(
        baseline_enrich,
        figures_dir / "enrichment_bar.png",
        top_n=config["visualization"]["top_terms"],
    )


if __name__ == "__main__":
    main()
