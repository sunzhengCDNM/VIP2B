#!/data/software/install/miniconda3/envs/python.3.7.0/bin//home/liujiang/software/install/mamba/envs/map2b/bin/python3

"""
数据库结构为
key： 酶切出来的序列
value: 序列id__序列方向__位置索引（序列位置索引从 0 开始，反向序列的位置从最后一个字符往前计算）

"""
########################################## import ################################################
import argparse, os, sys, marisa_trie, gzip, re, collections, math
from datetime import datetime
############################################ ___ #################################################
__doc__ = '多线程、多存储结构的电子酶切脚本'
__author__ = 'Liu Jiang'
__mail__ = 'jiang.liu@oebiotech.com'
__date__ = '2024/06/01 23:33:25'
__version__ = '1.0'
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

def open_file(sequence_file):
	if sequence_file.endswith('.gz'):
		return gzip.open(sequence_file, 'rt')
	else:
		return open(sequence_file, 'r')

def read_fafq(sequence_file):
	with open_file(sequence_file) as IN:
		try:
			first_line = IN.readline().strip()
		except UnicodeDecodeError:
			report('ERROR', '请给定 fa 或者 fq 格式的序列文件，不要提供 marisa 或者其他二进制的格式！')
	if first_line.startswith('>'):
		id = None
		seq = ''
		with open_file(sequence_file) as IN:
			for line in IN:
				line = line.strip()
				if line.startswith('>'):
					if id is not None:
						yield id, seq
					id = line.lstrip('>').split()[0]
					seq = ''
				else:
					seq += line.upper()
			if id is not None:
				yield id, seq
	elif first_line.startswith('@'):
		with open_file(sequence_file) as IN:
			while True:
				id = IN.readline().split()
				if not id:
					break  # End of file
				seq = IN.readline().strip()
				_ = IN.readline()  # Skip the separator line
				_ = IN.readline()  # Skip the quality score line
				yield id[0], seq
	else:
		report('ERROR', 'Unknown file format (not fastq or fasta)!')

def extraction(sequence_file, value_len, enzyme_pattern_sele_dic):
	global ori_seq_ct, dig_seq_ct
	ori_seq_ct = 0
	dig_seq_ct = 0
	for id, seq in read_fafq(sequence_file):
		cr_seq = seq[::-1].translate(str.maketrans('ACGTN', 'TGCAN'))
		for enzyme, enzyme_pattern in enzyme_pattern_sele_dic.items():
			ori_seq_ct += 1
			for matches in re.finditer(enzyme_pattern, seq):
				dig_seq_ct += 1
				new_id = '{}_{}_{}_{}'.format(id, '+', (matches.start() + 1), enzyme).rjust(value_len, '.')[-value_len:]
				yield matches.group(1), new_id
			for matches in re.finditer(enzyme_pattern, cr_seq):
				dig_seq_ct += 1
				new_id = '{}_{}_{}_{}'.format(id, '-', (matches.start() + 1), enzyme).rjust(value_len, '.')[-value_len:]
				yield matches.group(1), new_id

def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawTextHelpFormatter,
		epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,__mail__,__date__,__version__))
	parser.add_argument('-i',help='fa/fq文件或者数据库文件',dest='input',type=str,required=True)
	parser.add_argument('-e',help='酶切位点组合，可选 AlfI/AloI/BaeI/BcgI/BplI/BsaXI/BslFI/Bsp24I/CjeI/CjePI/CspCI/FalI/HaeIV/Hin4I/PpiI/PsrI 中的一个或多个，多个酶之间用逗号间隔，如果选择全部酶可直接提供 all',dest='enzyme',type=str, required=True)
	parser.add_argument('-o',help='输出文件前缀，建议使用样本名或者物种名',dest='output',type=str,required=False)
	parser.add_argument('-of',help='输出文件的压缩格式 marisa/gzip，默认为 marisa',dest='outfmt',type=str,choices=['marisa', 'gzip'], default='marisa')
	parser.add_argument('-l',help='marisa 值（id）长度，默认 50，建议不要更改，除非你对本程序及 marisa_trie 非常熟悉',dest='value_len',type=int,default=50)
	parser.add_argument('--dump',help='将 marisa 转换成 fasta 格式（标准输出）',dest='dump',action='store_true')
	args=parser.parse_args()
#	enzyme_length = enzyme_pattern_dic[args.enzyme][0]
#	fmt = '{}c'.format(enzyme_length)

	fmt = '{}c'.format(args.value_len)
	if not args.output and not args.dump:
		report('ERROR', '缺少参数：给定 -o 将结果输出到文件，或者给定 --dump 将 marisa 转换成 fasta')

	if args.dump:
		trie = marisa_trie.RecordTrie(fmt).mmap(args.input)
		for tag in set(trie.keys()):
			for id in trie[tag]:
				sys.stdout.write('>{}\n{}\n'.format(''.join([str(i, 'utf-8') for i in list(id)]).lstrip('.'), tag))
	else:
#		enzyme_list = args.enzyme.split(',')
		if args.enzyme == 'all':
			enzyme_pattern_sele_dic = enzyme_pattern_dic
		else:
			enzyme_pattern_sele_dic = {enzyme: enzyme_pattern_dic[enzyme] for enzyme in args.enzyme.split(',')}
#		[enzyme_pattern_dic[enzyme] for enzyme in args.enzyme.split(',')]
#		enzyme_pattern = enzyme_pattern_dic[enzyme_list[0]][1]
		if args.outfmt == 'marisa':
			marisa_key_lst, marisa_value_lst = [], []
			for seq, id in extraction(args.input, args.value_len, enzyme_pattern_sele_dic):
				marisa_key_lst.append(seq)
				marisa_value_lst.append(id.encode('utf-8'))
			trie = marisa_trie.BytesTrie(zip(marisa_key_lst, marisa_value_lst))
			trie.save('{}.fa.marisa'.format(args.output))
		else:
#			with gzip.open('{}.{}.fa.gz'.format(args.output, args.enzyme), 'wt') as OUT:
			with gzip.open('{}.fa.gz'.format(args.output), 'wt') as OUT:
				for seq, id in extraction(args.input, args.value_len, enzyme_pattern_sele_dic):
					OUT.write('>{}\n{}\n'.format(id.lstrip('.'), seq))
#		with open('{}.{}.dige.stat.xls'.format(args.output, args.enzyme), 'w') as STAT:
		with open('{}.dige.stat.xls'.format(args.output), 'w') as STAT:
			STAT.write('sample\tenzyme\tinput_sequence_num\tenzyme_reads_num\tpercent\n')
			STAT.write('{}\t{}\t{}\t{}\t{}\t'.format(args.output, args.enzyme, ori_seq_ct, dig_seq_ct, round((dig_seq_ct/ori_seq_ct), 2)))



if __name__=="__main__":
	main()
