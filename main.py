"""End-to-end entry point for the intelligent bio-computation workflow."""
from __future__ import annotations

import json
import logging
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


def main() -> None:
    configure_logging()
    project_root = Path(__file__).resolve().parent
    config = load_config(project_root / "config.yaml")
    paths = config["paths"]
    edge_path = project_root / paths["network_edges"]
    go_path = project_root / paths["go_annotations"]
    output_dir = project_root / paths["output_dir"]
    figures_dir = project_root / paths["figures_dir"]
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    data_loader.ensure_inputs(edge_path, go_path)
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
