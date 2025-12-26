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
    parser.add_argument("--downstream", action="store_true", help="Analyze CSV outputs in data/downstream and generate evaluations")
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
    downstream_dir = project_root / "data" / "downstream"
    # If downstream mode requested, handle it immediately (avoid prompting for edge file)
    if args.downstream:
        output_root = project_root / paths["output_dir"]
        summary_records = []
        if not downstream_dir.exists():
            raise FileNotFoundError(f"Downstream directory not found: {downstream_dir}")
        print(f"Running downstream evaluations for files in {downstream_dir}", flush=True)
        for csvf in sorted(downstream_dir.glob("*.csv")):
            try:
                print(f"Processing downstream file: {csvf.name}", flush=True)
                df_out = pd.read_csv(csvf)
                if "community_id" not in df_out.columns or "members" not in df_out.columns:
                    LOGGER.warning("Skipping %s: missing required columns 'community_id' and 'members'", csvf)
                    continue
                comms = {}
                for _, row in df_out.iterrows():
                    cid = str(row["community_id"])
                    mems = str(row["members"]).strip()
                    if "\t" in mems:
                        parts = mems.split("\t")
                    elif ";" in mems:
                        parts = [x.strip() for x in mems.split(";") if x.strip()]
                    else:
                        parts = [x for x in mems.split() if x]
                    comms[cid] = parts

                ds_out = output_root / "downstream" / csvf.stem
                ds_out.mkdir(parents=True, exist_ok=True)

                stats = {
                    "file": csvf.name,
                    "n_communities": len(comms),
                    "mean_community_size": float(pd.Series([len(v) for v in comms.values()]).mean() if comms else 0.0),
                }

                if config["enrichment"].get("enable_gseapy", False):
                    available_sets = config["enrichment"].get("gseapy_gene_sets", [])
                    go_sets = [gs for gs in available_sets if Path(gs).name.startswith("GO_")]
                    default_org = config.get("enrichment", {}).get("gseapy_organism", "Human")
                    organism = default_org
                    gene_sets = go_sets + [gs for gs in available_sets if "KEGG" in Path(gs).name and (organism in Path(gs).name or organism == "Human")]
                    print(f"Running enrichment for downstream file {csvf.name} (organism={organism})", flush=True)
                    enrich_df = enrichment_runner.enrich_communities(
                        comms,
                        gene_sets=gene_sets,
                        organism=organism,
                        cutoff=config["enrichment"].get("gseapy_cutoff", 0.05),
                        top_terms=config["enrichment"].get("gseapy_top_terms", 10),
                        out_file=ds_out / "enrichr_all.csv",
                        out_dir=ds_out / "per_community_enrich",
                        disable_progress=True,
                        silence_gseapy=True,
                    )
                    try:
                        summ = enrichment_runner.summarize_enrichment(enrich_df)
                        if not summ.empty:
                            summ.to_csv(ds_out / "enrichment_summary.csv", index=False)
                        concord = enrichment_runner.compute_concordance_scores(enrich_df, comms, cutoff=config["enrichment"].get("gseapy_cutoff", 0.05))
                        if not concord.empty:
                            concord.to_csv(ds_out / "enrichment_concordance.csv", index=False)
                            # find overall row where gene_set is NA
                            if "gene_set" in concord.columns:
                                overall = concord[concord["gene_set"].isna()]
                            else:
                                overall = concord.tail(1)
                            if not overall.empty:
                                orow = overall.iloc[0]
                                stats.update({
                                    "normalized_score": float(orow.get("normalized_score", float("nan"))),
                                    "frac_significant_communities": float(orow.get("frac_significant_communities", float("nan"))),
                                    "mean_score": float(orow.get("mean_score", float("nan"))),
                                    "weighted_score": float(orow.get("weighted_score", float("nan"))),
                                })
                    except Exception as exc:
                        LOGGER.warning("Downstream enrichment summarization failed for %s: %s", csvf, exc)
                else:
                    print("Enrichment disabled in config; skipping enrichment for downstream files", flush=True)

                pd.DataFrame([stats]).to_csv(ds_out / "summary_stats.csv", index=False)
                summary_records.append(stats)
            except Exception as exc:
                LOGGER.warning("Failed to process downstream file %s: %s", csvf, exc)
        # write comparison summary for all downstream files
        if summary_records:
            comp_out = output_root / "downstream" / "summary_comparison.csv"
            pd.DataFrame(summary_records).to_csv(comp_out, index=False)
            print(f"Wrote downstream comparison summary: {comp_out}", flush=True)
        print("Downstream evaluations complete", flush=True)
        return
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

    # If downstream mode requested, process all CSV outputs in data/downstream
    if args.downstream:
        if not downstream_dir.exists():
            raise FileNotFoundError(f"Downstream directory not found: {downstream_dir}")
        print(f"Running downstream evaluations for files in {downstream_dir}", flush=True)
        for csvf in sorted(downstream_dir.glob("*.csv")):
            try:
                print(f"Processing downstream file: {csvf.name}", flush=True)
                df_out = pd.read_csv(csvf)
                # Expect columns: community_id, members (space-separated or semicolon/comma)
                if "community_id" not in df_out.columns or "members" not in df_out.columns:
                    LOGGER.warning("Skipping %s: missing required columns 'community_id' and 'members'", csvf)
                    continue
                comms = {}
                for _, row in df_out.iterrows():
                    cid = str(row["community_id"])
                    mems = str(row["members"]).strip()
                    # support different separators
                    if "\t" in mems:
                        parts = mems.split("\t")
                    elif ";" in mems:
                        parts = [x.strip() for x in mems.split(";") if x.strip()]
                    else:
                        parts = [x for x in mems.split() if x]
                    comms[cid] = parts

                # prepare per-file output dir
                ds_out = output_root / "downstream" / csvf.stem
                ds_out.mkdir(parents=True, exist_ok=True)

                # basic stats
                stats = {
                    "file": csvf.name,
                    "n_communities": len(comms),
                    "mean_community_size": float(pd.Series([len(v) for v in comms.values()]).mean() if comms else 0.0),
                }
                pd.DataFrame([stats]).to_csv(ds_out / "summary_stats.csv", index=False)

                # If enrichment enabled, run enrichment and summarize
                if config["enrichment"].get("enable_gseapy", False):
                    # choose gene sets as in normal flow
                    available_sets = config["enrichment"].get("gseapy_gene_sets", [])
                    go_sets = [gs for gs in available_sets if Path(gs).name.startswith("GO_")]
                    # select organism default
                    default_org = config.get("enrichment", {}).get("gseapy_organism", "Human")
                    organism = default_org
                    gene_sets = go_sets + [gs for gs in available_sets if "KEGG" in Path(gs).name and (organism in Path(gs).name or organism == "Human")]
                    print(f"Running enrichment for downstream file {csvf.name} (organism={organism})", flush=True)
                    enrich_df = enrichment_runner.enrich_communities(
                        comms,
                        gene_sets=gene_sets,
                        organism=organism,
                        cutoff=config["enrichment"].get("gseapy_cutoff", 0.05),
                        top_terms=config["enrichment"].get("gseapy_top_terms", 10),
                        out_file=ds_out / "enrichr_all.csv",
                        out_dir=ds_out / "per_community_enrich",
                        disable_progress=True,
                        silence_gseapy=True,
                    )
                    # write summary and concordance if possible
                    try:
                        summ = enrichment_runner.summarize_enrichment(enrich_df)
                        if not summ.empty:
                            summ.to_csv(ds_out / "enrichment_summary.csv", index=False)
                        concord = enrichment_runner.compute_concordance_scores(enrich_df, comms, cutoff=config["enrichment"].get("gseapy_cutoff", 0.05))
                        if not concord.empty:
                            concord.to_csv(ds_out / "enrichment_concordance.csv", index=False)
                    except Exception as exc:
                        LOGGER.warning("Downstream enrichment summarization failed for %s: %s", csvf, exc)
                else:
                    print("Enrichment disabled in config; skipping enrichment for downstream files", flush=True)
            except Exception as exc:
                LOGGER.warning("Failed to process downstream file %s: %s", csvf, exc)
        print("Downstream evaluations complete", flush=True)
        return

    graph = preprocess.prepare_graph(
        edges,
        keep_giant_component=config["preprocessing"]["keep_giant_component"],
        weighted=config["preprocessing"]["weighted"],
    )
    print(f"Prepared graph: {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges", flush=True)

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
                enrichr_baseline_df = enrichment_runner.enrich_communities(
                    baseline_comm,
                    gene_sets=gene_sets,
                    organism=organism,
                    cutoff=config["enrichment"].get("gseapy_cutoff", 0.05),
                    top_terms=config["enrichment"].get("gseapy_top_terms", 10),
                    out_file=run_dir / "enrichr_baseline.csv",
                    out_dir=run_dir / "per_community_enrich" / "baseline",
                    disable_progress=local_gmt,
                    silence_gseapy=local_gmt,
                )
                enrichr_hier_df = enrichment_runner.enrich_communities(
                    hierarchical_comm,
                    gene_sets=gene_sets,
                    organism=organism,
                    cutoff=config["enrichment"].get("gseapy_cutoff", 0.05),
                    top_terms=config["enrichment"].get("gseapy_top_terms", 10),
                    out_file=run_dir / "enrichr_hierarchical.csv",
                    out_dir=run_dir / "per_community_enrich" / "hierarchical",
                    disable_progress=local_gmt,
                    silence_gseapy=local_gmt,
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
                # Compute concordance / evaluation scores comparing communities vs gene-sets
                try:
                    concord_base = enrichment_runner.compute_concordance_scores(enrichr_baseline_df, baseline_comm, cutoff=cutoff_val)
                    if not concord_base.empty:
                        concord_base.to_csv(run_dir / "enrichment_concordance_baseline.csv", index=False)
                    concord_hier = enrichment_runner.compute_concordance_scores(enrichr_hier_df, hierarchical_comm, cutoff=cutoff_val)
                    if not concord_hier.empty:
                        concord_hier.to_csv(run_dir / "enrichment_concordance_hierarchical.csv", index=False)
                    # Print brief summary
                    if not concord_base.empty:
                        overall = concord_base[concord_base["gene_set"].isna()].iloc[0]
                        print(f"Baseline concordance normalized score: {overall['normalized_score']:.3f}, frac significant communities: {overall['frac_significant_communities']:.3f}")
                    if not concord_hier.empty:
                        overall_h = concord_hier[concord_hier["gene_set"].isna()].iloc[0]
                        print(f"Hierarchical concordance normalized score: {overall_h['normalized_score']:.3f}, frac significant communities: {overall_h['frac_significant_communities']:.3f}")
                except Exception as exc:
                    LOGGER.warning("Failed to compute concordance scores: %s", exc)
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

        # Generate HTML report from template file (results/report_template.html)
        def _write_html_report(run_dir: Path):
                rep_path = run_dir / "report.html"
                per_base = run_dir / "per_community_enrich" / "baseline"
                per_hier = run_dir / "per_community_enrich" / "hierarchical"

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

                baseline_files = _list_csv(per_base)
                hier_files = _list_csv(per_hier)

                # read template from repository results/ folder
                tpl_path = project_root / "results" / "report_template.html"
                try:
                        tpl = tpl_path.read_text()
                except Exception:
                        LOGGER.warning("Report template not found (%s); skipping HTML report", tpl_path)
                        return

                import json
                baseline_json = json.dumps(baseline_files)
                hier_json = json.dumps(hier_files)
                html = tpl.replace('__BASELINE__', baseline_json).replace('__HIER__', hier_json).replace('__RUNNAME__', run_dir.name)

                try:
                        with open(rep_path, 'w') as fh:
                                fh.write(html)
                except Exception:
                        LOGGER.warning('Failed to write HTML report: %s', rep_path)

        _write_html_report(run_dir)


if __name__ == "__main__":
    main()
