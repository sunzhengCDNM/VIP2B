#!/data/software/install/miniconda3/envs/python.3.7.0/bin/python3

"""
数据库结构为
key： 8 位 genome_ID + 4 位 tag copy num(最大到 9999)
value: tag

v2.0 20230220
1. 修复基因组电子酶切时跳过重叠 tag 的 bug。

v2.1 20230223
1. 使用正则表达式的“正向前瞻”和“零宽度断言”功能进行重叠 tag 的识别。

"""
########################################## import ################################################
import argparse, os, sys, marisa_trie, gzip, re, collections, math
from datetime import datetime
############################################ ___ #################################################
__doc__ = 'This program is used for MAP2B database building'
__author__ = 'Liu Jiang'
__mail__ = 'jiang.liu@oebiotech.com'
__date__ = '2023/01/08 23:33:25'
__version__ = '2.1'
############################################ main ##################################################
enzyme_pattern_dic = {'BcgI':r'(?=([AGCT]{10}CGA[AGCT]{6}TGC[AGCT]{10}))',
					'CjePI':r'(?=([AGCT]{7}CCA[AGCT]{7}TC[AGCT]{8}))',
					'BsaXI':r'(?=([AGCT]{9}AC[AGCT]{5}CTCC[AGCT]{7}))'
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

def save_db(out, n, Tag_lst, ID_lst):
	trie = marisa_trie.BytesTrie(zip(ID_lst, Tag_lst))
	trie.save('{}.marisa{}'.format(out, n))

def read_fa(genome):
	seq_list = []
	seq = ''
	with gzip.open(genome, 'rt') as IN:
		for line in IN:
			line = line.strip()
			if line.startswith('>'):
				seq_list.append(seq)
				seq = ''
			else:
				seq += line.upper()
	seq_list.append(seq)
	return [seq for seq in seq_list if seq]

def extraction(genome, enzyme_pattern, tag_dic):
	for seq in read_fa(genome):
		for tag in re.findall(enzyme_pattern, seq):
			tag_dic.setdefault(tag, None)
		for tag in re.findall(enzyme_pattern, seq[::-1].translate(str.maketrans('ACGTN', 'TGCAN'))):
			tag_dic.setdefault(tag, None)
	return tag_dic

def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawTextHelpFormatter,
		epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,__mail__,__date__,__version__))
	parser.add_argument('-i',help='宿主基因组列表',dest='genome_list',type=str,required=True)
	parser.add_argument('-e',help='酶切位点，可选 BcgI / CjePI / BsaXI',dest='enzyme',type=str,choices=['BcgI', 'CjePI', 'BsaXI'], required=True)
	parser.add_argument('-o',help='输出数据库前缀, 例如 BcgI.human',dest='output',type=str,required=True)
	args=parser.parse_args()
	info = "runing..."
	report("INFO",info)

	enzyme_pattern = enzyme_pattern_dic[args.enzyme]
	tag_dic = {}
	with open(check_file(args.genome_list), 'r') as IN:
		for line in IN:
			line = line.strip()
			if line.startswith('#') or not line:continue
			tag_dic = extraction(check_file(line), enzyme_pattern, tag_dic)

	marisa_key_lst = tag_dic.keys()
	marisa_value_lst = ['0'.encode('utf-8')] * len(marisa_key_lst)
	trie = marisa_trie.BytesTrie(zip(marisa_key_lst, marisa_value_lst))
	trie.save('{}.marisa'.format(args.output))

if __name__=="__main__":
	main()
	sys.stdout.write('\n')
	info = "finish!"
	report("INFO",info)
