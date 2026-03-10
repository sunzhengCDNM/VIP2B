 #!/usr/bin/env python3
########################################## import ################################################
import argparse, os, sys, re, random, glob, gzip, marisa_trie, collections
from datetime import datetime
############################################ ___ #################################################
__doc__ = '本程序用于MAP2B流程的定量建库'
__author__ = 'Liu Jiang'
__mail__ = 'jiang.liu@oebiotech.com'
__date__ = '2022/11/24 22:10:41'
__version__ = '1.0.0'
############################################ main ##################################################
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
		return os.path.abspath(file)
	else:
		info = "file does not exist: {0}".format(file)
		report("ERROR",info)

def check_dir(dir):
	dir = os.path.abspath(dir)
	if not os.path.exists(dir):
		os.system("mkdir -p {0}".format(dir))
		info = "mkdir: {0}".format(dir)
		report("INFO",info)
	return(dir)


def ID_subdb(thd, ID_lst):
	subdb_ID_dic = {} # {classify:[ID]}
	for ID in sorted(ID_lst):
		subdb_ID_dic.setdefault((((int(ID) // thd) + 1) * thd), []).append(ID)
	return subdb_ID_dic

def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawTextHelpFormatter,
		epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,__mail__,__date__,__version__))
	parser.add_argument('-d',help='database prefix, like BcgI_CjePI',dest='database',type=str,required=True)
	parser.add_argument('-l',help='DB level, Class/Order/Family/Genus/Species/All',dest='level',
			choices=['Class', 'Order', 'Family', 'Genus', 'Species', 'All'], type=str,required=True)
	parser.add_argument('-c',help='classfy file',dest='classfy',type=str,required=True)
	parser.add_argument('-s',help='the size of the sub-databases when establishing the main database',dest='size',
			type=int,required=True)
	parser.add_argument('-o',help='outfile prefix, like BcgI.species.uniq',dest='prefix',type=str,required = True)
#	parser.add_argument('-e',help='enzyme, choose from BcgI / CjePI / BsaXI',dest='enzyme',
#			choices=['BcgI', 'CjePI', 'BsaXI'],type=str,required=True)
	parser.add_argument('-n',help='copy number for tag, choos from s[ingle] and m[ultiple]',dest='copy',
			choices=['s', 'm'],type=str,required=True)
	parser.add_argument('-p',help='pred.result',dest='pred',type=str,required=False)
	parser.add_argument('--intersection',help='intersection or union of tags between genomes, default union',
			dest='intersection',action='store_true')
	args=parser.parse_args()
	info = "runing..."
	report("INFO",info)

#	enzyme_dic = {
#			'CjePI':27,
#			'BcgI':32,
#			'BsaXI':27
#			}
#	outdir = check_dir(args.outdir)
	level_dic = {
			'Phylum':[3],
			'Class':[4],
			'Order':[5],
			'Family':[6],
			'Genus':[7],
			'Species':[8],
			'All':[3, 4, 5, 6, 7, 8]
			}
	level2_dic = {
			3:'Phylum',
			4:'Class',
			5:'Order',
			6:'Family',
			7:'Genus',
			8:'Species'
			}

	for level in level_dic[args.level]:
		report('INFO', 'Building a {}-level database...'.format(level2_dic[level]))
#		marisa_file = '{}/{}.{}.marisa'.format(outdir, args.enzyme, level2_dic[level])
#		stat_file = '{}/{}.{}.stat.xls'.format(outdir, args.enzyme, level2_dic[level])
#		fmt = '{}c'.format(enzyme_dic[args.enzyme])
		marisa_file = '{}.marisa'.format(args.prefix,  level2_dic[level])
		stat_file = '{}.stat.xls'.format(args.prefix, level2_dic[level])
		fmt = '40c'
		ID_taxo_dic = {} # {ID:taxo}
		spe_ID_dic = {} # {taxo:[ID]}

		if args.pred:
			taxo_pred_list = []
			with open(args.pred, 'r') as IN:
				for line in IN:
					line = line.strip()
					if line.startswith('Taxonomy') or not line:continue
					tmp = line.split('\t')
					if tmp[-2] == '1':
						taxo_pred_list.append(tmp[0])
			with gzip.open(check_file(args.classfy), 'rt') as IN:
				for line in IN:
					line = line.strip()
					if line.startswith('#') or not line:continue
					tmp = line.split('\t')
					taxo = ','.join(tmp[1:level])
					if taxo in taxo_pred_list:
						ID_taxo_dic.setdefault(tmp[0], taxo)
						spe_ID_dic.setdefault(taxo, []).append(tmp[0])
		else:
			with gzip.open(check_file(args.classfy), 'rt') as IN:
				for line in IN:
					line = line.strip()
					if line.startswith('#') or not line:continue
					tmp = line.split('\t')
					taxo = ','.join(tmp[1:level])
					ID_taxo_dic.setdefault(tmp[0], taxo)
					spe_ID_dic.setdefault(taxo, []).append(tmp[0])
		tag_sID_dic = {} # {tag:[sID]}
		subdb_ID_dic = ID_subdb(args.size, ID_taxo_dic.keys())
		ID_tagNum_dic = collections.defaultdict(int) # {ID:tag_num}
		for subdb_num, ID_lst in subdb_ID_dic.items():
			marisa = '{}.marisa{}'.format(args.database, subdb_num)
			trie = marisa_trie.RecordTrie(fmt).mmap(marisa)
			for ID in ID_lst:
				for sID in set(trie.keys(ID)):
					tag_tmp_list = trie[sID]
					ID_tagNum_dic[ID] += len(tag_tmp_list)
					for tag in tag_tmp_list:
						tag_sID_dic.setdefault(''.join([str(i, 'utf-8') for i in tag]), []).append(sID)

		ID_lst_trie, uniq_tag_lst_trie= [], []
		ID_theo_tag_num_dic = collections.defaultdict(int) # {ID: uniq_tag_num}
		if args.copy == 'm':
			for tag, sID_lst in tag_sID_dic.items():
				tmp_taxo_lst = [] # [taxo]
				for sID in sID_lst:
					tmp_taxo_lst.append(ID_taxo_dic[sID[:8]])
				if len(set(tmp_taxo_lst)) == 1:
					if args.intersection and len(list(set([sID[:8] for sID in sID_lst]))) != len(spe_ID_dic[tmp_taxo_lst[0]]):continue
					tmp_ID_lst_trie = []
					for sID in sID_lst:
						tmp_ID_lst_trie.append(sID[:8])
						ID_theo_tag_num_dic[sID[:8]] += int(sID[8:])
					ID_lst_trie += list(set(tmp_ID_lst_trie))
					uniq_tag_lst_trie += [tag] * len(set(tmp_ID_lst_trie))
		elif args.copy == 's':
			for tag, sID_lst in tag_sID_dic.items():
				tmp_taxo_lst = [] # [taxo]
				for sID in sID_lst:
					tmp_taxo_lst.append(ID_taxo_dic[sID[:8]])
				if len(set(tmp_taxo_lst)) == 1: # 该tag只出现在一个物种内
					if args.intersection and len(set([sID[:8] for sID in sID_lst])) != len(spe_ID_dic[tmp_taxo_lst[0]]):continue # tag出现的genome数不等于该物种下的genome数
					tmp_ID_lst_trie = [sID[:8] for sID in sID_lst if sID[8:] == '0001']
					for ID in tmp_ID_lst_trie:
						ID_theo_tag_num_dic[ID] += 1
					ID_lst_trie += list(set(tmp_ID_lst_trie))
					uniq_tag_lst_trie += [tag] * len(set(tmp_ID_lst_trie))

		del tag_sID_dic
		trie = marisa_trie.BytesTrie(zip(uniq_tag_lst_trie, [ID.encode('utf-8') for ID in ID_lst_trie]))
		del uniq_tag_lst_trie, ID_lst_trie
		trie.save(marisa_file)
		del trie

	## 统计文件
		taxo_theo_tag_num_dic = {} # {taxo: uniq_tag_num}
		for ID, sum_t in ID_theo_tag_num_dic.items():
			taxo_theo_tag_num_dic.setdefault(ID_taxo_dic[ID], []).append(sum_t)
		for taxo, t in taxo_theo_tag_num_dic.items():
			taxo_theo_tag_num_dic[taxo] = round(sum(t)/len(t), 4)
		with open(stat_file, 'w') as OUT:
			OUT.write('#Genome\tTaxonomy\tGenomeTheoUniqTag\tTaxoTheoUniqTag\tGenomeTheoTag\n')
			for ID in sorted(ID_theo_tag_num_dic.keys()):
				OUT.write('{}\t{}\t{}\t{}\t{}\n'.format(ID, ID_taxo_dic[ID], ID_theo_tag_num_dic[ID], taxo_theo_tag_num_dic[ID_taxo_dic[ID]], ID_tagNum_dic[ID]))

if __name__=="__main__":
	main()
	info = "finish!"
	report("INFO",info)
