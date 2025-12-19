# PPI Dataset Placeholder

Place your PPI edge list and GO annotations here, then point `config.yaml -> paths.datasets.ppi` to the filenames.
Expected:
- `ppi_network_edges.tsv`: tab-separated `gene_a gene_b [weight]`
- `ppi_go_annotations.tsv`: tab-separated `gene_id go_term description`
