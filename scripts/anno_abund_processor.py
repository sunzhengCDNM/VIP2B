#!/data/software/install/miniconda3/envs/python.3.7.0/bin/python3
import pandas as pd
import os, argparse
from pathlib import Path
from typing import Tuple, Dict, List

__doc__ = '''
  This script processes raw viral data(e.g., extracting target columns,
cleaning fields like viralverify_prediction, and linking protein annotations)
to build a standardized viral annotation database for downstream abundance analysis.
'''
__author__ = 'Liu Jiang'
__mail__ = 'liujiang9201@163.com'
__date__ = '2026/02/08 15:58:19'
__version__ = '1.0.0'

def init_output_dir(output_dir: str) -> Tuple[str, str, str]:
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    viral_tax_dir = os.path.join(output_dir, "viral_taxonomy")
    host_tax_dir = os.path.join(output_dir, "host_taxonomy")
    viral_func_dir = os.path.join(output_dir, "viral_function")
    Path(viral_tax_dir).mkdir(exist_ok=True)
    Path(host_tax_dir).mkdir(exist_ok=True)
    Path(viral_func_dir).mkdir(exist_ok=True)
    print(f"Output directory initialized successfully:")
    print(f"  Main directory: {output_dir}")
    print(f"  Viral taxonomy directory: {viral_tax_dir}")
    print(f"  Host taxonomy directory: {host_tax_dir}")
    print(f"  Viral function directory: {viral_func_dir}")
    return viral_tax_dir, host_tax_dir, viral_func_dir

def load_annotation_db(anno_path: str) -> Tuple[pd.DataFrame, set]:
    print("\nLoading annotation database...")
    required_cols = [
        "uhgv_votu", "is_jumbo_phage", "lifestyle", 
        "viralverify_prediction", "viral_taxonomy", "host_taxonomy", 
        "UniRef90_annotations", "info_annotations"
    ]
    df_anno = pd.read_csv(anno_path, sep='\t', compression='gzip')
    missing_cols = [col for col in required_cols if col not in df_anno.columns]
    if missing_cols:
        raise ValueError(f"Annotation database missing required columns: {missing_cols}")

    valid_jumbo = {"Yes", "No"}
    df_anno["is_jumbo_phage"] = df_anno["is_jumbo_phage"].apply(
        lambda x: x if x in valid_jumbo else "unknown"
    )
    valid_lifestyle = {"lytic", "temperate"}
    df_anno["lifestyle"] = df_anno["lifestyle"].apply(
        lambda x: x if x in valid_lifestyle else "unknown"
    )
    valid_viralverify = {"Chromosome", "Plasmid", "Uncertain", "Virus"}
    df_anno["viralverify_prediction"] = df_anno["viralverify_prediction"].apply(
        lambda x: x if x in valid_viralverify else "unknown"
    )

    anno_votu_set = set(df_anno["uhgv_votu"].dropna())
    print(f"Annotation database loaded successfully:")
    print(f"  Number of rows: {len(df_anno)}")
    print(f"  Number of unique uhgv_votu: {len(anno_votu_set)}")
    return df_anno[required_cols], anno_votu_set

def load_abundance_table(abund_path: str) -> Tuple[pd.DataFrame, List[str]]:
    print("\nLoading abundance table...")
    df_abund = pd.read_csv(abund_path, sep='\t')
    df_abund["uhgv_votu"] = df_abund["Taxonomy"].apply(
        lambda x: _extract_votu(x)
    )
    invalid_mask = df_abund["uhgv_votu"] == "missing_votu"
    if invalid_mask.sum() > 0:
        raise ValueError(
            f"Failed to extract vOTU from {invalid_mask.sum()} rows in abundance table (no 'vOTU-xxx' format in Taxonomy column)\n"
            f"Example rows:\n{df_abund[invalid_mask][['Taxonomy', 'uhgv_votu']].head(3).to_string(index=False)}"
        )
    sample_cols = [col for col in df_abund.columns if col not in ["Taxonomy", "uhgv_votu"]]
    if not sample_cols:
        raise ValueError("No sample columns found in abundance table (only Taxonomy column exists)")
    print(f"Abundance table loaded successfully:")
    print(f"  Number of rows: {len(df_abund)}")
    print(f"  Number of samples: {len(sample_cols)}")
    print(f"  Example sample columns: {sample_cols[:3]}")
    return df_abund[["uhgv_votu"] + sample_cols], sample_cols

def check_votu_match_and_merge(df_abund: pd.DataFrame, df_anno: pd.DataFrame, anno_votu_set: set) -> pd.DataFrame:
    print("\nChecking vOTU matching status...")
    abund_votu_set = set(df_abund["uhgv_votu"])
    unmatched_votu = abund_votu_set - anno_votu_set
    if unmatched_votu:
        unmatched_count = len(unmatched_votu)
        unmatched_ratio = (unmatched_count / len(abund_votu_set)) * 100
        raise ValueError(
            f"{unmatched_count} vOTUs in abundance table do not match annotation database ({unmatched_ratio:.1f}%)\n"
            f"Example unmatched vOTUs (first 10): {list(unmatched_votu)[:10]}\n"
            f"Please check: 1. Whether the annotation database contains these vOTUs; 2. Whether vOTU formats are consistent (e.g., vOTU-0001 vs vOTU_0001)"
        )
    df_merged = pd.merge(df_abund, df_anno, on="uhgv_votu", how="inner")
    print(f"vOTU matching completed successfully:")
    print(f"  Number of rows after merging: {len(df_merged)} (consistent with abundance table)")
    return df_merged

def generate_anno_and_func_files(df_merged: pd.DataFrame, sample_cols: list, output_dir: str, viral_func_dir: str) -> None:
    print("\nGenerating Phenotype.tsv...")
    fixed_anno_metrics = {
        "is_jumbo_phage": {
            "prefix": "is_jumbo_phage:",
            "values": ["Yes", "No", "unknown"]
        },
        "lifestyle": {
            "prefix": "lifestyle:",
            "values": ["lytic", "temperate", "unknown"]
        },
        "viralverify_prediction": {
            "prefix": "viralverify_prediction:",
            "values": ["Chromosome", "Plasmid", "Uncertain", "Virus", "unknown"]
        }
    }
    anno_result: Dict[str, Dict[str, float]] = {}

    for col, config in fixed_anno_metrics.items():
        prefix = config["prefix"]
        fixed_values = config["values"]
        for val in fixed_values:
            mask = df_merged[col] == val
            if mask.any():
                abund_sum = df_merged[mask][sample_cols].sum()
            else:
                abund_sum = pd.Series(0, index=sample_cols)
            metric_name = f"{prefix}{val}"
            anno_result[metric_name] = abund_sum.to_dict()

    df_annotation = pd.DataFrame.from_dict(anno_result, orient="index").sort_index()
    df_annotation.index.name = "Metric"
    pheno_path = os.path.join(output_dir, "Phenotype.tsv")
    df_annotation.to_csv(pheno_path, sep='\t', na_rep=0)
    print(f"Phenotype.tsv saved successfully: {pheno_path}")
    print(f"  Number of metrics: {len(df_annotation)}")
    print(f"  Metrics list: {list(df_annotation.index)}")

    print("\nGenerating Uniref90.tsv...")
    gene_abund_list = []
    for _, row in df_merged.iterrows():
        uniref_str = row["UniRef90_annotations"]
        if uniref_str in ["unknown", "-", ""] or pd.isna(uniref_str):
            continue
        genes = list(set([g.strip() for g in uniref_str.split("___") if g.strip()]))
        for gene in genes:
            gene_abund = {"gene": gene}
            gene_abund.update(row[sample_cols].to_dict())
            gene_abund_list.append(gene_abund)
    
    if not gene_abund_list:
        df_gene_norm = pd.DataFrame({"gene": []} | {col: [] for col in sample_cols})
        print("Warning: No UniRef90 gene annotations found, Uniref90.tsv will be empty")
    else:
        df_gene_abund = pd.DataFrame(gene_abund_list)
        df_gene_abund = df_gene_abund.groupby("gene")[sample_cols].sum().reset_index()
        df_gene_norm = df_gene_abund.copy()
        for col in sample_cols:
            total = df_gene_norm[col].sum()
            df_gene_norm[col] = df_gene_norm[col] / total if total > 0 else 0
    
    uniref_path = os.path.join(viral_func_dir, "Uniref90.tsv")
    df_gene_norm.to_csv(uniref_path, sep='\t', index=False, float_format="%.6f")
    print(f"Uniref90.tsv saved successfully: {uniref_path}")
    print(f"  Number of genes: {len(df_gene_norm)}")

    print("\nGenerating cluster.tsv...")
    cluster_abund_list = []
    for _, row in df_merged.iterrows():
        info_str = row["info_annotations"]
        if info_str in ["unknown", "-", ""] or pd.isna(info_str):
            continue
        clusters = list(set([c.strip() for c in info_str.split("___") if c.strip()]))
        for cluster in clusters:
            cluster_abund = {"cluster": cluster}
            cluster_abund.update(row[sample_cols].to_dict())
            cluster_abund_list.append(cluster_abund)
    
    if not cluster_abund_list:
        df_cluster_norm = pd.DataFrame({"cluster": []} | {col: [] for col in sample_cols})
        print("Warning: No info_annotations found, cluster.tsv will be empty")
    else:
        df_cluster_abund = pd.DataFrame(cluster_abund_list)
        df_cluster_abund = df_cluster_abund.groupby("cluster")[sample_cols].sum().reset_index()
        df_cluster_norm = df_cluster_abund.copy()
        for col in sample_cols:
            total = df_cluster_norm[col].sum()
            df_cluster_norm[col] = df_cluster_norm[col] / total if total > 0 else 0
    
    cluster_path = os.path.join(viral_func_dir, "cluster.tsv")
    df_cluster_norm.to_csv(cluster_path, sep='\t', index=False, float_format="%.6f")
    print(f"cluster.tsv saved successfully: {cluster_path}")
    print(f"  Number of clusters: {len(df_cluster_norm)}")

def generate_tax_abund_files(df_merged: pd.DataFrame, sample_cols: list, viral_tax_dir: str, host_tax_dir: str) -> None:
    tax_config = {
        "viral_taxonomy": [
            ['k__', 'p__', 'c__', 'o__', 'f__', 'g__'],
            viral_tax_dir
        ],
        "host_taxonomy": [
            ['d__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__'],
            host_tax_dir
        ]
    }

    for tax_col, (prefixes, save_dir) in tax_config.items():
        print(f"\nGenerating {tax_col} files (saved to: {save_dir})...")
        for level_idx, target_prefix in enumerate(prefixes):
            df_merged["current_tax_level"] = df_merged[tax_col].apply(
                lambda x: _get_tax_level(x, prefixes, level_idx, target_prefix)
            )
            df_level_abund = df_merged.groupby("current_tax_level")[sample_cols].sum().reset_index()
            df_level_abund.rename(columns={"current_tax_level": "Taxonomy"}, inplace=True)
            level_name_map = {"d": "domain", "k": "kingdom", "p": "phylum", "c": "class", 
                             "o": "order", "f": "family", "g": "genus", "s": "species"}
            level_name = level_name_map.get(target_prefix.rstrip('_'), target_prefix.rstrip('_'))
            file_name = f"{level_name}_abund.tsv"
            file_path = os.path.join(save_dir, file_name)
            df_level_abund.to_csv(file_path, sep='\t', index=False)
            print(f"  {file_name} saved successfully (number of taxa: {len(df_level_abund)})")
    if "current_tax_level" in df_merged.columns:
        df_merged.drop(columns="current_tax_level", inplace=True)

def _extract_votu(tax_str: str) -> str:
    if pd.isna(tax_str):
        return "missing_votu"
    parts = [p.strip() for p in str(tax_str).split(",") if "vOTU-" in p.strip()]
    return parts[-1] if parts else "missing_votu"

def _get_tax_level(tax_str: str, prefixes: list, level_idx: int, target_prefix: str) -> str:
    parts = [p.strip() for p in str(tax_str).split(";") if p.strip()]
    level_parts = []
    for part in parts:
        level_parts.append(part)
        if part.startswith(target_prefix):
            break
    while len(level_parts) < level_idx + 1:
        missing_prefix = prefixes[len(level_parts)]
        level_parts.append(f"{missing_prefix}unknown")
    return ";".join(level_parts)

def main():
    parser=argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter,
        epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(
            __author__,__mail__,__date__,__version__
            )
        )
    parser.add_argument('-i',help='Abundance.filtered.xls',dest='input',type=str,required=True)
    parser.add_argument('-d',help='metadata.tsv.gz',dest='metadata',type=str,required=True)
    parser.add_argument('-o',help='output dir',dest='output',type=str,required=True)
    args=parser.parse_args()

    print("=" * 60)
    print("Viral Annotation and Abundance Analysis Pipeline Started")
    print("=" * 60)

    viral_tax_dir, host_tax_dir, viral_func_dir = init_output_dir(args.output)

    df_anno, anno_votu_set = load_annotation_db(args.metadata)

    df_abund, sample_cols = load_abundance_table(args.input)

    df_merged = check_votu_match_and_merge(df_abund, df_anno, anno_votu_set)

    generate_anno_and_func_files(df_merged, sample_cols, args.output, viral_func_dir)

    generate_tax_abund_files(df_merged, sample_cols, viral_tax_dir, host_tax_dir)

    print("\n" + "=" * 60)
    print(f"All Files Generated Successfully! Final Results Saved to: {args.output}")
    print("=" * 60)

if __name__=="__main__":
    main()
