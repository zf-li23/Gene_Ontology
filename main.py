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
import shutil

import pandas as pd
import yaml

from src import data_loader
import re
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
    parser.add_argument("--dataset", choices=["sample", "scrin", "ppi"], help="Optional: dataset type to override config (sample/scrin/ppi)")
    parser.add_argument("--id-map", help="Optional TSV mapping file with columns 'source_id\tgene_symbol' to translate PPI IDs to gene symbols for enrichment")
    parser.add_argument("--organism", choices=["human", "mouse", "skip"], help="Optional: non-interactive organism choice for enrichment")
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--prune", type=float, default=0.3)
    parser.add_argument("--resolution", type=float, default=1.0)
    args = parser.parse_args()

    configure_logging()
    project_root = Path(__file__).resolve().parent
    config = load_config(project_root / "config.yaml")
    # allow overriding dataset type via CLI
    if args.dataset:
        config.setdefault("paths", {})["dataset"] = args.dataset
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
    # dataset-specific conversions
    if config.get("paths", {}).get("dataset") == "scrin" and chosen_edge.suffix.lower() == ".csv":
        base_stem = chosen_edge.stem  # keep original stem for naming outputs
        chosen_edge = data_loader.convert_scrin_csv_to_tsv(chosen_edge)
    if config.get("paths", {}).get("dataset") == "ppi":
        # support STRING/COG-style 'links' files (whitespace-separated combined_score)
        if chosen_edge.suffix.lower() in (".txt", ".gz") or "links" in chosen_edge.name:
            try:
                converted = data_loader.convert_string_links_to_tsv(chosen_edge)
                chosen_edge = converted
            except Exception:
                LOGGER.warning("Failed to convert PPI links file %s; expecting a TSV with gene_a/gene_b/weight", chosen_edge)

    go_path = dataset_paths.get("go_path")
    output_root = project_root / paths["output_dir"]
    run_dir = output_root / base_stem
    # If a previous run directory exists, remove it to ensure clean outputs for this run
    if run_dir.exists():
        try:
            shutil.rmtree(run_dir)
        except Exception:
            LOGGER.warning("Failed to remove existing run directory: %s", run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)

    if not chosen_edge.exists():
        raise FileNotFoundError(f"Edge file not found: {chosen_edge}")
    # Only auto-generate demo data when using the sample dataset; otherwise expect files to exist
    if config.get("paths", {}).get("dataset", "sample") == "sample" and not chosen_edge.exists():
        data_loader.ensure_inputs(chosen_edge, go_path)

    edges = data_loader.load_edge_list(chosen_edge, weighted=config["preprocessing"]["weighted"])

    # If user provided an ID mapping (e.g., COG -> gene symbol), apply it to edges
    if args.id_map:
        id_map_path = Path(args.id_map)
        if not id_map_path.exists():
            LOGGER.warning("ID map file not found: %s (skipping mapping)", id_map_path)
        else:
            try:
                map_df = pd.read_csv(id_map_path, sep="\t", header=None, names=["source", "target"] )
                id_map = dict(map_df.values)
                # map gene_a/gene_b
                edges["gene_a"] = edges["gene_a"].map(lambda x: id_map.get(x, x))
                edges["gene_b"] = edges["gene_b"].map(lambda x: id_map.get(x, x))
            except Exception as exc:
                LOGGER.warning("Failed to apply id map: %s", exc)

    # Note: GO annotation files are no longer used for hypergeometric tests.
    # Enrichment is performed via gseapy through `enrichment_runner` when enabled in config.

    # NOTE: algorithm integration will be invoked below as part of the standard run.

    # Prepare graph once for all algorithms
    graph = preprocess.prepare_graph(
        edges,
        keep_giant_component=config["preprocessing"]["keep_giant_component"],
        weighted=config["preprocessing"]["weighted"],
    )
    print(f"Prepared graph: {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges", flush=True)

    # Graph is already prepared above



    partition = run_louvain(
        graph,
        resolution=config["algorithms"]["louvain"]["resolution"],
        random_state=config["algorithms"]["louvain"].get("random_state"),
    )
    print(f"Baseline Louvain produced {len(set(partition.values()))} communities", flush=True)
    baseline_comm = partition_to_dict(partition, prefix="L")

    hierarchical = HierarchicalLouvain(
        min_size=config["algorithms"]["hierarchical"]["min_size"],
        max_depth=config["algorithms"]["hierarchical"]["max_depth"],
        modularity_threshold=config["algorithms"]["hierarchical"]["modularity_threshold"],
        resolution=config["algorithms"]["louvain"]["resolution"],
        random_state=config["algorithms"]["louvain"].get("random_state"),
    )
    hierarchy_root = hierarchical.fit(graph)
    print("Hierarchical detection complete", flush=True)
    hierarchical_comm = hierarchy_to_dict(hierarchy_root)

    # Run other algorithms
    other_results = {}
    try:
        from src.algorithms.algorithm_div import run_algorithm_div
        from src.algorithms.algorithm2 import run_algorithm2
        from src.algorithms.algorithm3 import run_algorithm3
        from src.algorithms.algorithm_spectral import run_algorithm_spectral
        
        # Div variants
        div_res = run_algorithm_div(graph)
        for name, part in div_res.items():
            # Convert partition to dict
            other_results[name] = partition_to_dict(part, prefix=f"{name}_")
            
        # Algo 2
        other_results["algorithm2"] = run_algorithm2(graph)
        
        # Algo 3
        other_results["algorithm3"] = run_algorithm3(graph)

        # Spectral
        try:
            other_results["spectral"] = run_algorithm_spectral(graph)
        except Exception as e:
            LOGGER.warning(f"Spectral failed: {e}")
        
    except Exception as e:
        LOGGER.warning(f"Failed to run other algorithms: {e}")

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
            # Always include GO gene sets (by basename starting with 'GO_'),
            # and include the organism-specific KEGG file depending on choice.
            available_sets = config["enrichment"].get("gseapy_gene_sets", [])
            go_sets = [gs for gs in available_sets if Path(gs).name.startswith("GO_")]
            if choice == "m" or choice == "mouse":
                organism = "Mouse"
                kegg_sets = [gs for gs in available_sets if "KEGG" in Path(gs).name and "Mouse" in Path(gs).name]
            elif choice == "s":
                organism = None
                kegg_sets = []
            else:
                organism = "Human"
                kegg_sets = [gs for gs in available_sets if "KEGG" in Path(gs).name and ("Human" in Path(gs).name or "Homo" in Path(gs).name)]
            # Combine GO sets + organism-specific KEGGs
            gene_sets = go_sets + kegg_sets

            if not organism:
                print("Skipping gseapy enrichment as requested.")
            else:
                print(f"Running enrichment for organism: {organism}. Gene sets: {gene_sets}")
                # Run enrichment on all communities (no hard limit); enrichment_runner shows per-community progress
                # detect whether gene sets refer to local GMT files -> disable progress and silence noisy gseapy logs
                local_gmt = any(Path(gs).exists() for gs in gene_sets)
                
                # Helper to run enrichment and save results
                def _process_enrichment(name, comms):
                    print(f"Enriching {name} ({len(comms)} communities)...")
                    df = enrichment_runner.enrich_communities(
                        comms,
                        gene_sets=gene_sets,
                        organism=organism,
                        cutoff=config["enrichment"].get("gseapy_cutoff", 0.05),
                        top_terms=config["enrichment"].get("gseapy_top_terms", 10),
                        out_file=run_dir / f"enrichr_{name}.csv",
                        out_dir=run_dir / "per_community_enrich" / name,
                        disable_progress=local_gmt,
                        silence_gseapy=local_gmt,
                    )
                    
                    # Summarize
                    try:
                        summ = enrichment_runner.summarize_enrichment(df)
                        if not summ.empty:
                            summ.to_csv(run_dir / f"enrichment_summary_{name}.csv", index=False)
                        concord = enrichment_runner.compute_concordance_scores(df, comms, cutoff=config["enrichment"].get("gseapy_cutoff", 0.05))
                        if not concord.empty:
                            concord.to_csv(run_dir / f"enrichment_concordance_{name}.csv", index=False)
                            overall = concord[concord["gene_set"].isna()]
                            if not overall.empty:
                                row = overall.iloc[0]
                                # Calculate gene coverage
                                all_genes_in_comms = set()
                                for genes in comms.values():
                                    all_genes_in_comms.update(genes)
                                coverage = len(all_genes_in_comms) / graph.number_of_nodes() if graph.number_of_nodes() > 0 else 0.0
                                print(f"{name} concordance normalized score: {row['normalized_score']:.3f}, gene coverage: {coverage:.3f}")
                    except Exception as e:
                        LOGGER.warning(f"Failed summary for {name}: {e}")
                    
                    return df

                # Collect all algorithms
                all_algorithms = {
                    "baseline": baseline_comm,
                    "hierarchical": hierarchical_comm,
                    **other_results
                }

                enrichment_results = {}
                for name, comms in all_algorithms.items():
                    enrichment_results[name] = _process_enrichment(name, comms)

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
            
            # Note: Per-community CSVs and summaries are already written by _process_enrichment
            
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

    # Export other algorithms community assignments
    for name, comms in other_results.items():
        pd.DataFrame([
            {"community_id": cid, "members": " ".join(members)} for cid, members in comms.items()
        ]).to_csv(run_dir / f"communities_{name}.csv", index=False)

    # Re-collect all algorithms if not already done (in case enrichment was skipped)
    all_algorithms = {
        "baseline": baseline_comm,
        "hierarchical": hierarchical_comm,
        **other_results
    }

    # Visualizations per run
    print("Generating visualizations for all algorithms...", flush=True)
    for name, comms in all_algorithms.items():
        # 1. Network Plot
        # Convert comms dict (cid->members) to partition dict (node->cid)
        partition = {}
        for cid, members in comms.items():
            for m in members:
                partition[m] = cid
        # Ensure all nodes are covered (assign to '0' or similar if missing, though they shouldn't be)
        # visualization.plot_network handles missing nodes gracefully
        visualization.plot_network(graph, partition, run_dir / f"network_{name}.png", layout=config["visualization"]["layout"])
        
        # 2. Enrichment Bar Plot (if available)
        # We need to reconstruct the enrichment dict from the results if they exist
        # Or if we are in the same scope, use enrichment_results
        # Since enrichment_results is local to the if block, we might need to reload or pass it out.
        # However, we can just check if the CSV exists or if we have the data.
        # Actually, let's use the enrichment_results dict if it exists in locals
        if 'enrichment_results' in locals() and name in enrichment_results:
             df = enrichment_results[name]
             if df is not None and not df.empty:
                 # We need the df_to_enrich_dict helper, which is also local.
                 # Let's redefine or move the helper to module level? 
                 # For now, I'll just duplicate the logic or assume it's available if I move the loop inside.
                 # But I want visualizations even if enrichment is disabled.
                 pass

    # Hierarchy plot is specific to hierarchical algo
    visualization.plot_hierarchy(hierarchy_root, run_dir / "hierarchy.png")
    
    # Plot enrichment bars if available
    if 'enrichment_results' in locals():
        cutoff_val = config["enrichment"].get("gseapy_cutoff", 0.05)
        for name, df in enrichment_results.items():
            if df is not None and not df.empty:
                # Convert DF to dict format expected by plot_enrichment_bar
                # We need to re-implement df_to_enrich_dict logic here or make it accessible
                enrich_dict = {}
                # ... (logic from df_to_enrich_dict) ...
                # Let's copy the logic briefly
                adj_candidates = ["adjusted_p_value", "Adjusted P-value", "adj_p"]
                term_candidates = ["term", "Term"]
                for cid, group in df.groupby("community_id"):
                    items = []
                    for _, row in group.iterrows():
                        term = None
                        for tcol in term_candidates:
                            if tcol in row.index and pd.notna(row.get(tcol)):
                                term = row.get(tcol)
                                break
                        q_val = None
                        for acol in adj_candidates:
                            if acol in row.index and pd.notna(row.get(acol)):
                                try: q_val = float(row.get(acol))
                                except: pass
                                break
                        if cutoff_val is not None and (q_val is None or q_val > cutoff_val): continue
                        if term: items.append({"term": term, "q_value": q_val if q_val is not None else 1.0})
                    if items: enrich_dict[str(cid)] = items
                
                if enrich_dict:
                    visualization.plot_enrichment_bar(
                        enrich_dict,
                        run_dir / f"enrichment_bar_{name}.png",
                        top_n=config["visualization"]["top_terms"],
                    )

    # Generate HTML report from template file (results/report_template.html)

        # Generate HTML report from template file (results/report_template.html)
        def _write_html_report(run_dir: Path):
                rep_path = run_dir / "report.html"
                per_comm_root = run_dir / "per_community_enrich"
                
                if not per_comm_root.exists():
                    return

                def _list_csv(folder: Path):
                    if not folder.exists():
                        return []
                    names = [p.name for p in folder.glob("*_enrich.csv")]
                    # natural sort by extracting integer components (e.g., L10 before L2 handled correctly)
                    def keyfn(name: str):
                        nums = re.findall(r"\d+", name)
                        if nums:
                            return tuple(int(x) for x in nums)
                        return (name,)
                    return sorted(names, key=keyfn)

                # Scan for all algorithm subdirectories
                algo_files = {}
                for subdir in per_comm_root.iterdir():
                    if subdir.is_dir():
                        algo_files[subdir.name] = _list_csv(subdir)

                # read template from repository results/ folder
                tpl_path = project_root / "results" / "report_template.html"
                try:
                        tpl = tpl_path.read_text()
                except Exception:
                        LOGGER.warning("Report template not found (%s); skipping HTML report", tpl_path)
                        return

                import json
                algo_files_json = json.dumps(algo_files)
                
                # Generate summaries HTML
                summaries_html = []
                for name in sorted(algo_files.keys()):
                    summaries_html.append(f'<li><strong>{name}</strong>: ')
                    links = []
                    # Check for summary files
                    if (run_dir / f"enrichr_{name}.csv").exists():
                        links.append(f'<a href="enrichr_{name}.csv">Full Results</a>')
                    if (run_dir / f"enrichment_summary_{name}.csv").exists():
                        links.append(f'<a href="enrichment_summary_{name}.csv">Summary</a>')
                    if (run_dir / f"enrichment_concordance_{name}.csv").exists():
                        links.append(f'<a href="enrichment_concordance_{name}.csv">Concordance</a>')
                    summaries_html.append(" | ".join(links))
                    summaries_html.append('</li>')
                
                html = tpl.replace('__ALGO_FILES__', algo_files_json).replace('__SUMMARIES_HTML__', "\n".join(summaries_html)).replace('__RUNNAME__', run_dir.name)

                try:
                        with open(rep_path, 'w') as fh:
                                fh.write(html)
                except Exception:
                        LOGGER.warning('Failed to write HTML report: %s', rep_path)

        _write_html_report(run_dir)


if __name__ == "__main__":
    main()
