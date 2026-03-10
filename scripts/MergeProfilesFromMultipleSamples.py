#!/usr/bin/env python3
import os
import sys
import argparse
from pathlib import Path

def usage():
    """显示帮助信息"""
    help_msg = """
DESCRIPTION
    Merge the abundance profiles from multiple samples. 
    Can filter specific taxa and normalize remaining abundances when -f is specified.
USAGE
    python3 {script}
    
PARAMETERS
    -l <file> A list file indicating the sample_id and corresponding output files from the last step.
              e.g., sample_id<tab>qualitative/sample_id/sample_id.combine.xls
    -o <file>  The output file path
OPTIONS
    -f <str>   Taxon to filter (last part of Taxonomy column, comma-separated). 
               When provided, output file will contain filtered taxa removed and abundances normalized to 1.
    -h|--help   help
AUTHOR: Modified from LJ 2025.10.17
""".format(script=os.path.basename(sys.argv[0]))
    print(help_msg, file=sys.stderr)

def check_dir_for_file(path):
    """检查文件所在目录是否存在，不存在则创建"""
    dir_path = os.path.dirname(path)
    if dir_path and not os.path.isdir(dir_path):
        try:
            os.makedirs(dir_path, exist_ok=True)
        except OSError as e:
            print(f"Error: Directory {dir_path} cannot be created: {e}", file=sys.stderr)
            sys.exit(1)
    return True

def normalize_abundances(abundance_dict, samples):
    """将丰度值归一化，使每个样本的总和为1"""
    # 计算每个样本的当前总和
    sample_sums = {sample: 0.0 for sample in samples}
    for id_str, sample_abunds in abundance_dict.items():
        for sample in samples:
            sample_sums[sample] += sample_abunds.get(sample, 0.0)
    
    # 归一化每个样本的丰度值
    normalized = {}
    for id_str, sample_abunds in abundance_dict.items():
        normalized_abunds = {}
        for sample in samples:
            total = sample_sums[sample]
            if total == 0:
                normalized_abunds[sample] = 0.0
            else:
                normalized_abunds[sample] = sample_abunds.get(sample, 0.0) / total
        normalized[id_str] = normalized_abunds
    
    return normalized

def main():
    # 解析命令行参数
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-l", type=str, help="List file")
    parser.add_argument("-o", type=str, help="Output file path")
    parser.add_argument("-f", type=str, help="Taxon to filter (last part of Taxonomy)", default="human")
    parser.add_argument("-h", "--help", action="store_true", help="Show help")
    args = parser.parse_args()

    # 处理帮助信息
    if args.help:
        usage()
        sys.exit(0)

    # 检查必要参数
    if not all([args.l, args.o]):
        usage()
        print("Error: Parameters -l and -o are required", file=sys.stderr)
        sys.exit(1)

    # 初始化变量
    list_file = args.l
    outfile = args.o  # 输出文件路径
    filter_taxon = args.f  # 要过滤的分类单元（Taxonomy列最后一部分）

    # 确保输出目录存在
    check_dir_for_file(outfile)

    # 读取样品列表和数据
    hash_specie = {}  # 存储物种在各样本中的数量 {物种ID: {样本名: 数量}}
    hash_all = {}     # 存储每个样本的总数量 {样本名: 总数}
    sample_sort = []  # 保持样本顺序
    classify_col = -1 # 分类列索引
    head = ""         # 表头

    try:
        with open(list_file, 'r') as li:
            for line in li:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) != 2:
                    continue
                sample, path = parts
                path = os.path.abspath(path)
                
                if not os.path.exists(path):
                    print(f"Warning: {sample} {path} not exist, cannot calculate Abundance.", file=sys.stderr)
                    continue
                
                sample_sort.append(sample)
                
                try:
                    with open(path, 'r') as infile:
                        for line_in in infile:
                            line_in = line_in.strip()
                            if not line_in:
                                continue
                            tmp = line_in.split('\t')
                            
                            # 处理表头行
                            if tmp[0].startswith("Taxonomy"):
                                for i, col in enumerate(tmp):
                                    if col == "Theoretical_Tag_Num":
                                        classify_col = i - 1
                                        head = '\t'.join(tmp[:classify_col+1])  # 分类列部分的表头
                                        break
                            elif classify_col == -1:
                                continue  # 还没找到分类列，跳过
                            else:
                                # 处理数据行
                                if len(tmp) <= classify_col:
                                    continue  # 数据不完整，跳过
                                id_parts = tmp[:classify_col+1]
                                id_str = '\t'.join(id_parts)
                                
                                # 获取数量值（倒数第4列）
                                if len(tmp) < 4:
                                    continue
                                count = float(tmp[-4]) if tmp[-4] else 0.0
                                
                                # 更新物种数据
                                if id_str not in hash_specie:
                                    hash_specie[id_str] = {}
                                hash_specie[id_str][sample] = count
                                
                                # 更新样本总数
                                if sample in hash_all:
                                    hash_all[sample] += count
                                else:
                                    hash_all[sample] = count
                except IOError as e:
                    print(f"Error opening {path}: {e}", file=sys.stderr)
                    continue
    except IOError as e:
        print(f"Error opening {list_file}: {e}", file=sys.stderr)
        sys.exit(1)

    # 确定要输出的数据（原始或过滤后）
    if filter_taxon:
        # 筛选需要保留的物种（排除Taxonomy最后一部分匹配filter_taxon的）
        filtered_specie = {}
        for id_str in hash_specie:
            # 假设Taxonomy列是id_str的第一部分
            taxonomy_part = id_str.split('\t')[0]  # 获取Taxonomy列
            last_taxon = taxonomy_part.split(',')[-1].strip()
            if last_taxon != filter_taxon:
                filtered_specie[id_str] = hash_specie[id_str]
        
        # 对过滤后的丰度进行归一化（每个样本总和为1）
        output_data = normalize_abundances(filtered_specie, sample_sort)
    else:
        # 使用原始数据（计算百分比）
        output_data = {}
        for id_str, sample_abunds in hash_specie.items():
            percent_abunds = {}
            for sample in sample_sort:
                count = sample_abunds.get(sample, 0.0)
                total = hash_all.get(sample, 0)
                percent = count / total if total != 0 else 0.0
                percent_abunds[sample] = percent
            output_data[id_str] = percent_abunds

    # 输出结果到指定文件
    try:
        with open(outfile, 'w') as ou:
            sample_header = '\t'.join(sample_sort)
            ou.write(f"{head}\t{sample_header}\n")
            
            # 写入数据
            for id_str in sorted(output_data.keys()):
                judge = 0
                print_line = id_str
                for sample in sample_sort:
                    percent = output_data[id_str].get(sample, 0.0)
                    if percent == 0:
                        print_line += "\t0"
                    else:
                        print_line += f"\t{percent}"
                        judge += 1
                if judge != 0:
                    ou.write(f"{print_line}\n")
    except IOError as e:
        print(f"Error writing to {outfile}: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
