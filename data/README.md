# Data Layout

- `sample/` toy PPI + GO (auto-generated if missing when `paths.dataset=sample`).
- `PPI/` real protein-protein interactions + GO annotations. Expect TSV edge list (`gene_a gene_b [weight]`).
- `scrin/` mRNA共定位相互作用数据与对应注释。

Update `config.yaml -> paths.datasets` to point to your files when switching datasets.
