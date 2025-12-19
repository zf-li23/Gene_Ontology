# Intelligent Bio Computation Pipeline

Unsupervised community detection and functional validation on gene/protein interaction networks. Supports standard Louvain + GO enrichment as well as a fast, overlapping, multi-level detector tuned for large PPIs.

## Features
- Data loaders for STRING/BioGRID-style edge lists and GO annotations (tab-delimited).
- Preprocessing: weighted graphs, giant-component trimming, weight normalization.
- Baseline Louvain and hierarchical Louvain (python-louvain).
- Overlapping, multi-level detector (isolate-partition + parallel Louvain + edge-seeded expansion) with topology self-checks.
- GO enrichment (hypergeometric + FDR), evaluation metrics, and matplotlib visualizations.
- Config-driven defaults plus CLI mode for custom edge lists.

## Project Layout
```
intelligent_bio_pipeline/
├── README.md              # You are here
├── config.yaml            # Paths and parameters for the config-driven pipeline
├── data/                  # Sample toy PPI and GO files
├── results/               # Metrics and figures (created at runtime)
├── src/
│   ├── algorithms/        # Louvain baseline + hierarchical wrappers
│   ├── data_loader.py     # Download/load/generate interaction + GO data
│   ├── preprocess.py      # Graph assembly and cleanup helpers
│   ├── enrichment_analysis.py # Hypergeometric + FDR
│   ├── evaluation.py      # Coverage and GO-consistency metrics
│   ├── visualization.py   # Network/hierarchy/enrichment plots
│   └── overlap_detector.py# Overlapping, multi-level detector + topology checks
└── main.py                # Entry point (config mode or CLI edge-file mode)
```

## Quickstart
Install core deps (example):
```
pip install networkx numpy numba community pandas pyyaml matplotlib scipy statsmodels requests
```

### A) Config-driven pipeline (with GO enrichment)
Runs on the dataset selected in `config.yaml` (`paths.dataset`). Presets:
- `sample`: data/sample/*.tsv (auto-generated if missing)
- `ppi`: data/PPI/ppi_network_edges.tsv, ppi_go_annotations.tsv
- `scrin`: data/scrin/scrin_edges.tsv, scrin_go_annotations.tsv (mRNA共定位)

Execute and write metrics/figures to `results/`:
```
python main.py
```
Outputs include `results/metrics.json`, enrichment CSVs, and figures in `results/figures/`.

### B) Overlapping hierarchical detector (unsupervised, GO-free)
Run on a custom edge list (TSV: gene_a gene_b [weight]) and emit JSON to stdout:
```
python main.py data/sample_network_edges.tsv --threads 8 --prune 0.3 --resolution 1.0 > results/hierarchy.json
```
- Supports overlapping memberships and multi-level tree output.
- Uses isolate-sets for parallel Louvain, edge-centrality seeding, and topology fingerprints (C(k) slope, clustering, modularity).

## Data Notes
- Sample files live in `data/sample/`. Real datasets go into `data/PPI/` or `data/scrin/`; update `config.yaml` via the `paths.datasets` map.
- `data_loader.ensure_inputs` will synthesize a demo network only when `paths.dataset=sample`; for `ppi` or `scrin`, the files must exist.

## Collaboration (GitHub)
1) Add remote (once):
```
git remote add origin https://github.com/<org_or_user>/<repo>.git
```
2) Pull latest (to sync with teammates):
```
git pull --rebase origin main
```
3) After edits, commit and push:
```
git add .
git commit -m "<message>"
git push origin main
```

## Performance Tips
- For ≥10k-node PPIs, use the overlapping detector with `--threads` to parallelize isolate-sets.
- Ensure `numba` is installed to accelerate isolate-set construction; adjust `--prune` for denser/sparser communities.
- If you see Numba JIT errors, upgrade it (e.g., `pip install -U numba`) and rerun; the isolate-set kernel avoids generators to stay nopython-safe.
