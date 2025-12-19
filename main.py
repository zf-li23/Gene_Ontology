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
from src import visualization
from src import enrichment_runner
from src.overlap_detector import run_pipeline as run_overlap


LOGGER = logging.getLogger("intelligent_bio")


def configure_logging() -> None:
    # Keep library logs quiet; provide concise CLI output via prints
    logging.basicConfig(level=logging.WARNING, format="%(asctime)s | %(levelname)s | %(name)s | %(message)s")
    logging.getLogger("gseapy").setLevel(logging.WARNING)
    logging.getLogger("urllib3").setLevel(logging.WARNING)


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
            "go_path": (project_root / entry["go_annotations"]) if entry.get("go_annotations") else None,
        }
    # legacy fallback
    return {
        "edge_path": project_root / paths_cfg["network_edges"],
        "go_path": project_root / paths_cfg.get("go_annotations") if paths_cfg.get("go_annotations") else None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Intelligent bio-computation pipeline")
    parser.add_argument("edge_file", nargs="?", help="Optional: run overlapping hierarchical detector on a custom edge list")
    parser.add_argument("--organism", choices=["human", "mouse", "skip"], help="Optional: non-interactive organism choice for enrichment")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--prune", type=float, default=0.3)
    parser.add_argument("--resolution", type=float, default=1.0)
    args = parser.parse_args()

    configure_logging()
    project_root = Path(__file__).resolve().parent
    config = load_config(project_root / "config.yaml")
    paths = config["paths"]
    dataset_paths = resolve_dataset_paths(config, project_root)
    # interactive edge file choice when not provided via CLI
    if args.edge_file:
        chosen_edge = Path(args.edge_file)
    else:
        default_edge = dataset_paths["edge_path"]
        user_input = input(f"Edge file path [default: {default_edge}]: ").strip()
        chosen_edge = Path(user_input) if user_input else default_edge

    base_stem = chosen_edge.stem
    # handle scrin CSV -> TSV conversion automatically
    if chosen_edge.suffix.lower() == ".csv":
        base_stem = chosen_edge.stem  # keep original stem for naming outputs
        chosen_edge = data_loader.convert_scrin_csv_to_tsv(chosen_edge)

    go_path = dataset_paths.get("go_path")
    output_root = project_root / paths["output_dir"]
    run_dir = output_root / base_stem
    run_dir.mkdir(parents=True, exist_ok=True)

    if not chosen_edge.exists():
        raise FileNotFoundError(f"Edge file not found: {chosen_edge}")
    # Only auto-generate demo data when using the sample dataset; otherwise expect files to exist
    if config.get("paths", {}).get("dataset", "sample") == "sample" and not chosen_edge.exists():
        data_loader.ensure_inputs(chosen_edge, go_path)

    edges = data_loader.load_edge_list(chosen_edge, weighted=config["preprocessing"]["weighted"])

    # Note: GO annotation files are no longer used for hypergeometric tests.
    # Enrichment is performed via gseapy through `enrichment_runner` when enabled in config.

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

    # Optional Enrichr/GSEAPY enrichment for external databases (e.g., KEGG, GO)
    baseline_enrich = {}
    hierarchical_enrich = {}

    if config["enrichment"].get("enable_gseapy", False):
        try:
            # Prompt user for organism choice (Human or Mouse) for enrichment
            print(f"Detected {len(baseline_comm)} baseline communities and {len(hierarchical_comm)} hierarchical communities.")
            import sys
            # Determine choice: CLI flag > interactive prompt > config default
            if args.organism:
                choice = args.organism[0]
            else:
                if sys.stdin.isatty():
                    choice = input("Run enrichment for Human (h), Mouse (m), or skip (s) [h/m/s, default h]: ").strip().lower()
                else:
                    # non-interactive: use config default organism
                    default_org = config.get("enrichment", {}).get("gseapy_organism", "Human")
                    choice = default_org[0].lower() if default_org else "h"
            if choice == "m" or choice == "mouse":
                organism = "Mouse"
                gene_sets = [gs for gs in config["enrichment"].get("gseapy_gene_sets", []) if "Mouse" in gs or gs.startswith("GO_")]
            elif choice == "s":
                organism = None
                gene_sets = []
            else:
                organism = "Human"
                gene_sets = [gs for gs in config["enrichment"].get("gseapy_gene_sets", []) if "Human" in gs or gs.startswith("GO_")]

            if not organism:
                print("Skipping gseapy enrichment as requested.")
            else:
                print(f"Running enrichment for organism: {organism}. Gene sets: {gene_sets}")
                # Run enrichment on all communities (no hard limit); enrichment_runner shows per-community progress
                enrichr_baseline_df = enrichment_runner.enrich_communities(
                    baseline_comm,
                    gene_sets=gene_sets,
                    organism=organism,
                    cutoff=config["enrichment"].get("gseapy_cutoff", 0.05),
                    top_terms=config["enrichment"].get("gseapy_top_terms", 10),
                    out_file=run_dir / "enrichr_baseline.csv",
                    out_dir=run_dir / "per_community_enrich" / "baseline",
                )
                enrichr_hier_df = enrichment_runner.enrich_communities(
                    hierarchical_comm,
                    gene_sets=gene_sets,
                    organism=organism,
                    cutoff=config["enrichment"].get("gseapy_cutoff", 0.05),
                    top_terms=config["enrichment"].get("gseapy_top_terms", 10),
                    out_file=run_dir / "enrichr_hierarchical.csv",
                    out_dir=run_dir / "per_community_enrich" / "hierarchical",
                )

            def df_to_enrich_dict(df, cutoff=None):
                out = {}
                if df is None or df.empty:
                    return out
                # possible column names for adjusted p-value
                adj_candidates = ["adjusted_p_value", "Adjusted P-value", "Adjusted Pvalue", "adj_p", "Adj P", "Adjusted P-value "]
                term_candidates = ["term", "Term"]
                for cid, group in df.groupby("community_id"):
                    items = []
                    for _, row in group.iterrows():
                        term = None
                        for tcol in term_candidates:
                            if tcol in row.index and pd.notna(row.get(tcol)):
                                term = row.get(tcol)
                                break
                        # find adjusted p-value
                        q_val = None
                        for acol in adj_candidates:
                            if acol in row.index and pd.notna(row.get(acol)):
                                try:
                                    q_val = float(row.get(acol))
                                except Exception:
                                    q_val = None
                                break
                        # if cutoff provided, enforce significance threshold
                        if cutoff is not None:
                            if q_val is None or q_val > cutoff:
                                continue
                        # require a term name
                        if term is None:
                            continue
                        items.append({"term": term, "q_value": q_val if q_val is not None else 1.0})
                    if items:
                        out[str(cid)] = items
                return out

            cutoff_val = config["enrichment"].get("gseapy_cutoff", 0.05)
            baseline_enrich = df_to_enrich_dict(enrichr_baseline_df, cutoff=cutoff_val)
            hierarchical_enrich = df_to_enrich_dict(enrichr_hier_df, cutoff=cutoff_val)
            if not baseline_enrich:
                LOGGER.info("Enrichr returned no significant terms for baseline communities")

            # Ensure per-community CSVs exist (some runs may skip writing them inside runner)
            def _write_per_comm(df, subdir_name: str):
                if df is None or df.empty:
                    return
                outp = run_dir / "per_community_enrich" / subdir_name
                outp.mkdir(parents=True, exist_ok=True)
                for cid, group in df.groupby("community_id"):
                    group.to_csv(outp / f"{cid}_enrich.csv", index=False)

            try:
                _write_per_comm(enrichr_baseline_df, "baseline")
                _write_per_comm(enrichr_hier_df, "hierarchical")
            except Exception:
                pass

            # Summarize across communities (aggregate p-values) and write summary CSVs
            try:
                summary_base = enrichment_runner.summarize_enrichment(enrichr_baseline_df)
                if not summary_base.empty:
                    summary_base.to_csv(run_dir / "enrichment_summary_baseline.csv", index=False)
                summary_hier = enrichment_runner.summarize_enrichment(enrichr_hier_df)
                if not summary_hier.empty:
                    summary_hier.to_csv(run_dir / "enrichment_summary_hierarchical.csv", index=False)
            except Exception as exc:
                LOGGER.warning("Failed to create enrichment summary: %s", exc)
        except ImportError:
            LOGGER.warning("gseapy not installed; skipping Enrichr enrichment")
        except Exception as exc:  # pragma: no cover - runtime guard
            LOGGER.warning("Enrichr failed: %s", exc)

    # Export community assignments
    pd.DataFrame([
        {"community_id": cid, "members": " ".join(members)} for cid, members in baseline_comm.items()
    ]).to_csv(run_dir / "communities_louvain.csv", index=False)

    pd.DataFrame([
        {"community_id": cid, "members": " ".join(members)} for cid, members in hierarchical_comm.items()
    ]).to_csv(run_dir / "communities_hierarchical.csv", index=False)

    # Visualizations per run
    visualization.plot_network(graph, partition, run_dir / "network.png", layout=config["visualization"]["layout"])
    visualization.plot_hierarchy(hierarchy_root, run_dir / "hierarchy.png")
    # Plot enrichment bar if any Enrichr results exist
    if baseline_enrich:
        visualization.plot_enrichment_bar(
            baseline_enrich,
            run_dir / "enrichment_bar.png",
            top_n=config["visualization"]["top_terms"],
        )


if __name__ == "__main__":
    main()
