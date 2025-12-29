# Intelligent Bio Pipeline (智能生物计算流程)

这是一个集成了多种社区检测算法与基因集富集分析的综合生物信息学流水线。旨在从蛋白质相互作用网络（PPI）或其他生物网络中识别功能模块，并自动进行生物学意义注释。

## 主要功能

*   **多算法集成**：
    *   **Baseline Louvain**: 标准 Louvain 算法。
    *   **Hierarchical Louvain**: 层次化 Louvain 算法。
    *   **Divisive Louvain**: 包含 4 种加权变体 (`div_0_5`, `div_1`, `div_2`, `div_avg`)，模拟 C++ 版本的节点度加权逻辑，并包含微小社区过滤。
    *   **Spectral Clustering**: 自动估计 $k$ 值的谱聚类算法，包含 Eigengap Heuristic。
    *   **Algorithm 2 & 3**: 基于 DFS 的社区检测算法，包含全局合并（Global Merge）策略以消除冗余社区，并具备未分配节点的后处理机制以提升覆盖率。
*   **自动富集分析**：集成 `gseapy`，支持 GO (BP, CC, MF) 和 KEGG 通路富集分析。
*   **可视化报告**：自动生成网络结构图（节点大小代表加权度，颜色代表社区）、富集条形图，并汇总为交互式 HTML 报告。
*   **灵活配置**：支持通过 `config.yaml` 配置或命令行参数覆盖。

## 安装指南

推荐使用 Conda 环境：

```bash
conda create -n genomics python=3.11 -y
conda activate genomics
pip install -r requirements.txt
```

## 快速开始

使用默认配置运行示例数据：

```bash
python main.py
```

或者指定一个边缘列表文件：

```bash
python main.py data/scrin/test.csv --dataset scrin
```

## 命令行使用详解

`main.py` 是程序的统一入口，支持以下参数：

### 基础用法
```bash
python main.py [EDGE_FILE] [OPTIONS]
```

### 参数说明

*   **`EDGE_FILE`** (可选位置参数):
    *   指定输入的网络边缘列表文件路径（支持 CSV, TSV, TXT）。
    *   如果不提供，程序将使用 `config.yaml` 中定义的默认路径，或进入交互模式询问路径。

*   **`--dataset`** (可选):
    *   指定数据集类型，用于通过特定逻辑解析输入文件。
    *   可选值:
        *   `sample`: 默认格式。
        *   `scrin`: 针对 SCRIN 数据集的格式（自动处理 CSV 转 TSV）。
        *   `ppi`: 针对 PPI 网络（如 STRING 数据库）的格式。

*   **`--organism`** (可选):
    *   指定富集分析的目标物种。
    *   可选值:
        *   `human` (或 `h`): 人类 (Homo sapiens)。
        *   `mouse` (或 `m`): 小鼠 (Mus musculus)。
        *   `skip` (或 `s`): 跳过富集分析步骤。
    *   如果不指定，程序将在运行时交互式询问。

*   **`--id-map`** (可选):
    *   指定 ID 映射文件路径（TSV 格式，两列：`source_id`, `gene_symbol`）。
    *   用于将网络中的 ID（如 COG ID）转换为基因符号，以便进行富集分析。

*   **`--threads`** (可选):
    *   指定并行处理的线程数（默认为 8）。

*   **`--resolution`** (可选):
    *   Louvain 算法的分辨率参数（默认为 1.0）。值越大，社区越小；值越小，社区越大。

### 示例

1.  **运行 SCRIN 数据集并指定小鼠富集**：
    ```bash
    python main.py data/scrin/test.csv --dataset scrin --organism mouse
    ```

2.  **运行 PPI 数据并提供 ID 映射**：
    ```bash
    python main.py data/ppi/network.txt --dataset ppi --id-map data/ppi/id_mapping.tsv
    ```

## 输出结果

运行完成后，结果将保存在 `results/<文件名>/` 目录下：

*   **`report.html`**: **[核心]** 交互式 HTML 报告，汇总了所有算法的富集结果和链接。
*   **`communities_*.csv`**: 各算法生成的社区划分列表。
*   **`network_*.png`**: 各算法的社区网络可视化图。
*   **`enrichment_bar_*.png`**: 富集分析结果的前 10 个术语条形图。
*   **`enrichr_*.csv`**: 详细的富集分析统计表。
*   **`enrichment_summary_*.csv`**: 富集结果的摘要统计。
*   **`enrichment_concordance_*.csv`**: 社区与功能术语的一致性评分，包含 **Normalized Score** (归一化一致性得分) 和 **Gene Coverage** (基因覆盖率) 等关键指标。

## 目录结构

*   `src/`: 源代码目录。
    *   `algorithms/`: 包含所有社区检测算法实现。
*   `data/`: 存放输入数据和 GMT 基因集文件。
*   `results/`: 存放运行结果。
*   `config.yaml`: 全局配置文件（路径、预处理参数等）。

## 查看交互式报告（report.html）

生成的 `report.html` 位于某次运行的输出目录 `results/<run_name>/report.html` 中。推荐使用 HTTP 服务器来查看（避免直接用 `file://` 导致某些浏览器安全限制）。

本地查看（在本地机器上）：

```bash
# 切换到该次运行的结果目录
cd results/<run_name>
# 启动一个简单的 HTTP 服务器（Python 3）
python -m http.server 8000
# 在浏览器中打开：http://localhost:8000/report.html
```

远程服务器（在远程无头机上运行并通过端口转发查看）：

方法 A — 使用 SSH 端口转发（推荐）：

```bash
# 在本地终端建立隧道（将远程 8000 转发到本地 8000）
ssh -L 8000:localhost:8000 user@remote.host
# 在远程服务器上（ssh 后）运行：
cd /path/to/repo/results/<run_name>
python -m http.server 8000 --bind 127.0.0.1
# 在本地浏览器打开 http://localhost:8000/report.html
```

方法 B — 在远程机器上对外绑定（不推荐在不受信任网络上使用）：

```bash
cd /path/to/repo/results/<run_name>
# 绑定到 0.0.0.0，允许远程主机直接访问（注意防火墙/安全）
python -m http.server 8000 --bind 0.0.0.0
# 然后在本机浏览器打开 http://remote.host:8000/report.html
```

提示与注意事项：
- 如果 `report.html` 引用了本地 GMT 文件或其它数据资源，确保这些文件在 `results/<run_name>/` 的相对路径下或使用绝对路径可访问。
- 在无图形界面的服务器上运行富集（`gseapy`）可能需要网络访问；若无网络，请在 `config.yaml` 中配置本地 GMT 并将 `enable_gseapy` 设置为 `true` 并指定 GMT 路径。

## 常见问题与故障排查

- 如果网页在加载数据时显示空表或图像缺失，检查 `results/<run_name>/per_community_enrich/` 是否包含每个社区对应的 CSV 文件（`<community>_enrich.csv`）。
- 若出现 `ModuleNotFoundError`（例如 `gseapy`、`community` 等），请确认已按 `requirements.txt` 安装依赖，或使用推荐的 Conda 环境：

```bash
conda activate genomics
pip install -r requirements.txt
```

- 如果谱聚类输出警告 “Graph is not fully connected”，这是因为输入图有多个连通分量；建议在大型或非连通图上采用分量逐个聚类或先进行连通分量分解。

## 联系与支持

如需帮助，请在项目的 issue 中描述重现步骤或将运行日志贴上来，我们会协助定位问题。
