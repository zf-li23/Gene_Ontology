# Intelligent Bio Computation Pipeline

Unsupervised community detection and functional validation on gene/protein interaction networks. Supports standard Louvain + GO enrichment as well as a fast, overlapping, multi-level detector tuned for large PPIs.

## Features
- Data loaders for STRING/BioGRID-style edge lists and GO annotations (tab-delimited).
- Preprocessing: weighted graphs, giant-component trimming, weight normalization.
- Baseline Louvain and hierarchical Louvain (python-louvain).
- Overlapping, multi-level detector (isolate-partition + parallel Louvain + edge-seeded expansion) with topology self-checks.
- GO enrichment (hypergeometric + FDR) and matplotlib visualizations.
- Config-driven defaults plus CLI mode for custom edge lists.

## Project Layout
# Intelligent Bio Pipeline

这是一个集成了网络社区检测与基因集富集分析的流水线示例，包含：

- Louvain 与层次化 Louvain 社区检测
- 基于 `gseapy` / Enrichr 的 GO / KEGG 富集（支持本地 GMT）
- 每次运行会在 `results/<edge_stem>/` 产出表格、图像，并生成交互式 HTML 报告用于算法评价

下面是从刚 clone 仓库到查看报告的手把手说明（假设你已经 cd 到 `intelligent_bio_pipeline` 目录）：

## 1) 创建 & 激活 Python 环境（推荐 conda）

推荐使用 conda（与测试环境一致）：

```bash
conda create -n genomics python=3.11 -y
conda activate genomics
```

如果不使用 conda，也可以用 venv：

```bash
python -m venv .venv
source .venv/bin/activate
```

## 2) 安装依赖

```bash
pip install --upgrade pip
pip install -r requirements.txt
```

如果你使用 conda 并偏好 conda 包，可以按需用 `conda install` 安装对应包。

## 3) 下载 GMT 数据库（推荐，本地化可显著加速富集）

仓库含 `download_GMT.py`（用于从 Enrichr/MAayan 下载并生成 GMT 文件），运行：

```bash
# 将数据库写到 data/gmt
python download_GMT.py --outdir data/gmt
```

下载后 `data/gmt/` 下会包含 KEGG 与 GO 的 GMT 文件，主流程会优先使用本地 GMT。

## 4) 准备要分析的边表（edge list）

把你要检测的边表拷贝到 `data/scrin/` 下（示例使用 `test.csv`）：

```bash
cp /path/to/my_edges.csv data/scrin/test.csv
```

注意：如果输入是特定 SCRIN 风格 CSV，脚本会自动做必要的格式转换。

## 5) 配置（可选）

编辑 `config.yaml` 来修改：
- `paths.dataset`（scrin / ppi / sample）
- `enrichment.gseapy_gene_sets`（可填本地 GMT 路径 `data/gmt/*.gmt`）
- `enrichment.gseapy_cutoff`（显著性阈值，例如 0.05）

## 6) 运行主流水线（示例命令与参数）

基本运行（交互选择或按 config 默认）：

```bash
python main.py
```

使用指定边文件：

```bash
python main.py data/scrin/test.csv
```

常用 CLI 参数：

- `edge_file`（位置参数）: 指定边表文件路径
- `--organism` : `human` | `mouse` | `skip`（决定是否附加相应 KEGG 文件；GO GMT 始终包含）
- `--threads`  : 并行线程数（默认 8）
- `--prune`    : 剔除阈值（默认 0.3）
- `--resolution`: Louvain 分辨率（默认 1.0）

示例（非交互、指定 organism 与线程）：

```bash
python main.py data/scrin/test.csv --organism mouse --threads 4 --prune 0.25 --resolution 1.0
```

运行行为摘要：

- 会删除并重新创建 `results/<edge_stem>/`（保证输出干净）
- 执行 Louvain 與 Hierarchical Louvain 社区检测
- 使用 `gseapy`（本地 GMT 或 Enrichr API）对每个社区进行富集
- 写出 per-community CSV、合并富集表、汇总（Fisher）及 concordance 评分
- 生成交互式报告 `results/<edge_stem>/report.html`

## 7) 查看报告（建议用本地静态服务器）

浏览器直接打开 `file://.../report.html` 可能因跨域/安全限制导致 fetch 失败，建议在报告目录启动临时 HTTP 服务器：

```bash
cd results/test
python -m http.server 8000
# 在浏览器打开：http://localhost:8000/report.html
```

或在项目根直接指定目录：

```bash
python -m http.server --directory results/test 8000
```

## 8) 主要输出文件说明

- `results/<stem>/communities_louvain.csv` — baseline 社区成员
- `results/<stem>/communities_hierarchical.csv` — 层级社区成员
- `results/<stem>/per_community_enrich/baseline/*.csv` — baseline 每社区富集结果
- `results/<stem>/per_community_enrich/hierarchical/*.csv` — hierarchical 每社区富集结果
- `results/<stem>/enrichr_baseline.csv`, `enrichr_hierarchical.csv` — 合并富集表
- `results/<stem>/enrichment_summary_*.csv` — Fisher 合并 p 值汇总
- `results/<stem>/enrichment_concordance_*.csv` — 算法与基因库一致性指标
- `results/<stem>/report.html` — 交互式 HTML 报告

## 常见问题

- 如果页面加载 per-community CSV 失败，优先检查是否以 `file://` 打开：请使用 `python -m http.server`。 
- 若缺少 `gseapy`，请安装：

```bash
pip install gseapy
```

- 若需要把报告打包为单文件（内嵌图片与 CSV），或需要自定义前端展示（更多统计图），我可以帮你生成自包含 HTML 或改进前端脚本。

---

如需我替你把当前 `results/<stem>/` 打包为单个 HTML 或 ZIP，或把 README 翻译成英文版，告诉我即可。祝使用顺利！
