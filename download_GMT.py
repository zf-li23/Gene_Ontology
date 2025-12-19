import gseapy as gp
import os

os.makedirs('./data/gmt', exist_ok=True)

libraries = [
    'GO_Biological_Process_2025',
    'GO_Cellular_Component_2025', 
    'GO_Molecular_Function_2025',
    'KEGG_2019_Mouse',
    'KEGG_2021_Human'
]

for lib in libraries:
    print(f'正在下载: {lib}')
    try:
        # 1. 获取数据 (新版本API)
        gmt_data = gp.get_library(name=lib)
        # 2. 手动写入GMT文件
        output_path = f'./data/gmt/{lib}.gmt'
        with open(output_path, 'w') as f:
            # gmt_data 是一个字典，键是基因集名，值是基因列表
            for gene_set, genes in gmt_data.items():
                # GMT格式：基因集名[TAB]描述[TAB]基因1[TAB]基因2...
                # 这里用占位符作为描述，也可以写 lib 或其他信息
                line = f"{gene_set}\t{lib}\t" + "\t".join(genes) + "\n"
                f.write(line)
        print(f'已保存到: {output_path}')
        
    except Exception as e:
        print(f'下载 {lib} 失败: {e}')