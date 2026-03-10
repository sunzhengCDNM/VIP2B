#!/data/software/install/miniconda3/envs/python.3.7.0/bin/python3

"""
功能：对多个输入的fa.gz文件进行酶切片段提取，去冗余后以fasta格式输出
修改说明：
1. 支持输入多个用逗号分隔的fa.gz文件
2. 去除冗余酶切片段（全局去重）
3. 直接以fasta格式存储，无需补充固定长度
4. 移除marisa_trie数据库相关功能
"""
########################################## import ################################################
import argparse, os, sys, gzip, re, collections
from datetime import datetime
############################################ ___ #################################################
__doc__ = 'This program is used for enzyme digestion fragment extraction and save as fasta'
__author__ = 'Liu Jiang (modified)'
__mail__ = 'jiang.liu@oebiotech.com'
__date__ = '2023/01/08 23:33:25 (modified)'
__version__ = '2.3'
############################################ main ##################################################
enzyme_pattern_dic = {
					'AlfI':r'(?=([AGCT]{10}GCA[AGCT]{6}TGC[AGCT]{10}))',
					'AloI':r'(?=([AGCT]{7}GAAC[AGCT]{6}TCC[AGCT]{7}))',
					'BaeI':r'(?=([AGCT]{10}AC[AGCT]{4}GTA[CT]C[AGCT]{7}))',
					'BcgI':r'(?=([AGCT]{10}CGA[AGCT]{6}TGC[AGCT]{10}))',
					'BplI':r'(?=([AGCT]{8}GAG[AGCT]{5}CTC[AGCT]{8}))',
					'BsaXI':r'(?=([AGCT]{9}AC[AGCT]{5}CTCC[AGCT]{7}))',
					'BslFI':r'(?=([AGCT]{6}GGGAC[AGCT]{14}))',
					'Bsp24I':r'(?=([AGCT]{8}GAC[AGCT]{6}TGG[AGCT]{7}))',
					'CjeI':r'(?=([AGCT]{8}CCA[AGCT]{6}GT[AGCT]{9}))',
					'CjePI':r'(?=([AGCT]{7}CCA[AGCT]{7}TC[AGCT]{8}))',
					'CspCI':r'(?=([AGCT]{11}CAA[AGCT]{5}GTGG[AGCT]{10}))',
					'FalI':r'(?=([AGCT]{8}AAG[AGCT]{5}CTT[AGCT]{8}))',
					'HaeIV':r'(?=([AGCT]{7}GA[CT][AGCT]{5}[AG]TC[AGCT]{9}))',
					'Hin4I':r'(?=([AGCT]{8}GA[CT][AGCT]{5}[GAC]TC[AGCT]{8}))',
					'PpiI':r'(?=([AGCT]{7}GAAC[AGCT]{5}CTC[AGCT]{8}))',
					'PsrI':r'(?=([AGCT]{7}GAAC[AGCT]{6}TAC[AGCT]{7}))',
					}

def report(level,info):
	date_now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
	if level == "ERROR":
		sys.stderr.write("{0} - {1} - ERROR - {2}\n".format(date_now,os.path.basename(__file__),info))
		sys.exit(1)
	elif level == "INFO":
		sys.stdout.write("{0} - {1} - INFO - {2}\n".format(date_now,os.path.basename(__file__),info))
	elif level == "DEBUG":
		sys.stdout.write("{0} - {1} - DEBUG - {2}\n".format(date_now,os.path.basename(__file__),info))
		sys.exit(1)
	return

def check_file(file):
	if os.path.exists(file):
		return file
	else:
		info = "file does not exist: {0}".format(file)
		report("ERROR",info)

def read_fa(genome_file):
	"""读取单个fa.gz文件中的序列"""
	seq_list = []
	seq = ''
	with gzip.open(genome_file, 'rt') as IN:
		for line in IN:
			line = line.strip()
			if line.startswith('>'):
				if seq:  # 避免添加空序列
					seq_list.append(seq)
				seq = ''
			else:
				seq += line.upper()
		if seq:  # 添加最后一条序列
			seq_list.append(seq)
	return seq_list

def extraction(fasta_files, enzyme_pattern_list):
	"""提取多个fasta文件的酶切片段并去冗余"""
	global_tag_set = set()  # 全局去重集合
	file_tag_map = []  # 存储(片段, 来源文件)
	
	# 遍历所有输入文件
	for file_idx, fasta_file in enumerate(fasta_files, 1):
		report("INFO", f"处理文件 {file_idx}/{len(fasta_files)}: {fasta_file}")
		seq_list = read_fa(fasta_file)
		
		# 处理每条序列的正向和反向互补链
		for seq in seq_list:
			# 正向链
			for pattern in enzyme_pattern_list:
				for tag in re.findall(pattern, seq):
					if tag and tag not in global_tag_set:
						global_tag_set.add(tag)
						file_tag_map.append((tag, os.path.basename(fasta_file)))
			
			# 反向互补链
			rev_comp_seq = seq[::-1].translate(str.maketrans('ACGTN', 'TGCAN'))
			for pattern in enzyme_pattern_list:
				for tag in re.findall(pattern, rev_comp_seq):
					if tag and tag not in global_tag_set:
						global_tag_set.add(tag)
						file_tag_map.append((tag, os.path.basename(fasta_file)))
	
	return file_tag_map

def progress_bar(p):
	sys.stdout.write('\r{}{} {}%'.format('#' * p, '-' * (50 - p), p * 2))
	sys.stdout.flush()

def save_fasta(output, tag_info):
	"""保存为fasta格式，标题包含片段ID和来源文件"""
	with open(output, 'w') as out_fh:
		for idx, (tag, filename) in enumerate(tag_info, 1):
			out_fh.write(f'>fragment_{idx}_from_{filename}\n')
			out_fh.write(f'{tag}\n')

def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawTextHelpFormatter,
		epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,__mail__,__date__,__version__))
	parser.add_argument('-i',help='输入的fa.gz文件，多个文件用逗号分隔',dest='input_files',type=str,required=True)
	parser.add_argument('-e',help='酶切位点组合，可选 AlfI/AloI/.../PsrI 中的一个或多个，多个酶之间用逗号间隔，选择全部酶可提供 all',dest='enzyme',type=str, required=True)
	parser.add_argument('-o',help='输出的fasta文件路径',dest='output',type=str,required=True)
	args=parser.parse_args()
	
	report("INFO", "开始运行酶切片段提取...")
	
	# 解析输入文件列表
	input_files = [check_file(f) for f in args.input_files.split(',')]
	report("INFO", f"共检测到 {len(input_files)} 个输入文件")
	
	# 解析酶切位点列表
	if args.enzyme == 'all':
		enzyme_pattern_list = list(enzyme_pattern_dic.values())
	else:
		enzyme_pattern_list = [enzyme_pattern_dic[enzyme] for enzyme in args.enzyme.split(',')]
	report("INFO", f"使用的酶切位点: {args.enzyme}")
	
	# 提取并去冗余酶切片段
	tag_info = extraction(input_files, enzyme_pattern_list)
	
	# 保存结果
	save_fasta(args.output, tag_info)
	report("INFO", f"共提取到 {len(tag_info)} 个独特的酶切片段")
	report("INFO", f"结果已保存至 {args.output}")

if __name__=="__main__":
	main()
	sys.stdout.write('\n')
	report("INFO", "处理完成!")
