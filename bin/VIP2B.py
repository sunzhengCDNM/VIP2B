#!python3
########################################## import ################################################
import argparse, os, sys, re, random, glob, gzip
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor
src_dir = os.path.abspath(os.path.dirname(__file__) + '/../scripts')
def_db = os.path.abspath(os.path.dirname(__file__) + '/../database/8Enzyme')
config_dir = os.path.abspath(os.path.dirname(__file__) + '/../config')
############################################ ___ #################################################
__doc__ = ''
__author__ = 'Zheng Sun, Jiang Liu'
__mail__ = 'sunzheng0618@gmail.com'
__date__ = '2026/3/6 22:02:47'
__version__ = '1.1'
"""
update 20250303 v1.0.1
1. rename VIP2B

updata 20250709 v1.0.3
1. upgrade the machine learning algorithm to XGBoost.

updata 20250904 v1.0.4
1. updated the XGBoost model.

updata 20251017 v1.0.5
1. updated database.

updata 20260306 v1.1
1. added function annotations and other extended annotations.
2. Fixed other minor bugs.
"""
############################################ main ##################################################
def report(level, info):
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
	return dir

def check_db(db_dir):
	db_lst = '{config_dir}/def_db.list'.format(config_dir = config_dir)
	t = 0
	if not os.path.exists(db_lst):
		report('ERROR', 'Unable to find {db_lst}, please download from github'.format(db_lst = db_lst))
	else:
		with open(db_lst, 'r') as IN:
			for line in IN:
				line = line.strip()
				if line.startswith('#') or not line:continue
				tmp = line.split('\t')
				db_file = '{}/{}'.format(db_dir, tmp[0])
				if not os.path.exists(db_file):
					t = 1
	if t == 1:
		report('INFO', 'Downloading database, due to network reasons, may take a long time, please wait')
		exe_shell(f'python3 {src_dir}/DownloadDB.py -l {db_lst} -d {db_dir}', 'DownloadDB')
	return

def exe_shell(cmd, desc=None):
	desc = None # Disable the 'desc'
	if desc:
		report('INFO', 'Running: {}'.format(desc))
	else:
		report('INFO', 'Running: {}'.format(cmd))
#	print(cmd)
	code = os.system(cmd)
	if code == 0:
		return
	else:
		if desc:
			report('ERROR', 'failed to run: {}'.format(desc))
		else:
			report('ERROR', 'failed to run: {}'.format(cmd))

def mkdb(db_dir, smp, quan_db, cutoff, level):
	exe_shell('python3 {src_dir}/CreatDB4AllLevel.py -d {db} -c {db_dir}/abfh_classify_with_speciename.txt.gz -p {O}/1.qual/{smp}/pred.result -o {quan_db}/{smp} -s {cutoff} -n m -l {level} {intersection}'.format(src_dir = src_dir, db = args.database, db_dir = db_dir, smp = smp, quan_db = quan_db, O = O, cutoff = cutoff, level = level, intersection = "--intersection" if args.intersection else ''), 'creatQuanDB')
	return

def check_data(data_file):
	data_dic = {}
	with open(data_file,'r') as IN:
		for line in IN:
			line = line.strip()
			if line.startswith('#') or not line:continue
			tmp = line.split('\t')
			data_dic.setdefault(tmp[0], [])
			for f in [check_file(i) for i in tmp[1:] if i]:
				data_dic.setdefault(tmp[0], []).append(f)
	for smp, data_list in data_dic.items():
			if len(data_list) > 2:
				report('ERROR', 'Sample {0} has more than 2 Reads, or sample {0} is named duplicate: {1}'.format(smp, ', '.join(data_list)))
			elif len(data_list) == 0:
				report('ERROR', 'Sample {} has no Reads'.format(smp))
	return data_dic

def cc_abd(db, abfh, enzyme_smp_file, O_, processes, level, threshold):
	exe_shell('python3 {src_dir}/CalculateRelativeAbundance_Single2bEnzyme.py -d {db} -c {abfh} -l {enzyme_smp_file} -o {O} -p {processes} -t {level} -ct {threshold}'.format(src_dir = src_dir, db = db, abfh = abfh, enzyme_smp_file = enzyme_smp_file, O = O_, processes = processes, level = level, threshold = threshold), 'CalculateRelativeAbundance_Single2bEnzyme')

def extra_tag(reads, enzyme_dir, smp, db_dir):
	if len(reads) == 2:
		check_dir(enzyme_dir)
		exe_shell('python3 {src_dir}/sequence_digestion.py -i {reads} -e {enzyme} -o {enzyme_dir}/{smp}_1 -of gzip'.format(src_dir = src_dir, reads = reads[0], enzyme = enzyme, enzyme_dir = enzyme_dir, smp = smp), '2bRADExtraction')
		exe_shell('python3 {src_dir}/sequence_digestion.py -i {reads} -e {enzyme} -o {enzyme_dir}/{smp}_2 -of gzip'.format(src_dir = src_dir, reads = reads[1], enzyme = enzyme, enzyme_dir = enzyme_dir, smp = smp), '2bRADExtraction')
		exe_shell('cat {enzyme_dir}/{smp}_1.fa.gz {enzyme_dir}/{smp}_2.fa.gz >{enzyme_dir}/{smp}.all.fa.gz'.format(enzyme = enzyme, enzyme_dir = enzyme_dir, smp = smp), 'mergePEReads')
		with open('{enzyme_dir}/{smp}_1.dige.stat.xls'.format(enzyme_dir = enzyme_dir, smp = smp), 'r') as IN:
			for line in IN:
				tmp = line.strip().split('\t')
				if tmp[0] == 'sample':continue
				rr = int(tmp[2])
				er = int(tmp[3])
		with open('{enzyme_dir}/{smp}_2.dige.stat.xls'.format(enzyme_dir = enzyme_dir, smp = smp), 'r') as IN:
			for line in IN:
				tmp = line.strip().split('\t')
				if tmp[0] == 'sample':continue
				rr += int(tmp[2])
				er += int(tmp[3])
		with open('{enzyme_dir}/{smp}.dige.stat.xls'.format(enzyme_dir = enzyme_dir, smp = smp), 'w') as OUT:
			OUT.write('sample\tenzyme\tinput_reads_num\tenzyme_reads_num\tpercent\n')
			OUT.write('{sample}\t{enzyme}\t{input_reads_num}\t{enzyme_reads_num}\t{percent}%\n'.format(sample = smp, enzyme = enzyme, input_reads_num = rr, enzyme_reads_num = er, percent = (round((er / rr), 4) * 100)))
		exe_shell('rm {enzyme_dir}/{smp}_1.* {enzyme_dir}/{smp}_2.*'.format(enzyme_dir = enzyme_dir, smp = smp), 'clean_tmp')
	else:
		check_dir(enzyme_dir)
		exe_shell('python3 {src_dir}/sequence_digestion.py -i {reads} -e {enzyme} -o {enzyme_dir}/{smp} -of gzip && mv {enzyme_dir}/{smp}.fa.gz {enzyme_dir}/{smp}.all.fa.gz'.format(src_dir = src_dir, reads = reads[0], enzyme = enzyme, enzyme_dir = enzyme_dir, smp = smp), '2bRADExtraction')
#	if args.host:
#		exe_shell('python3 {src_dir}/host_filter.py -i {enzyme_dir}/{smp}.all.fa.gz -d {db_dir}/ -e {enzyme} -s {host} -o {enzyme_dir}/{smp}'.format(src_dir = src_dir, enzyme = enzyme, enzyme_dir = enzyme_dir, smp = smp, host = args.host, db_dir = db_dir), 'host_filter')
#	else:
#		exe_shell('cd {enzyme_dir} && ln -sf {smp}.all.fa.gz {smp}.fa.gz'.format(enzyme = enzyme, enzyme_dir = enzyme_dir, smp = smp))
	exe_shell('cd {enzyme_dir} && ln -sf {smp}.all.fa.gz {smp}.fa.gz'.format(enzyme = enzyme, enzyme_dir = enzyme_dir, smp = smp))
	return

def run_dige(db_dir, data_dic):
	enzyme_dir = check_dir(O + '/0.dige')
	done_file = O + '/0.dige/done'
	enzyme_smp_file = O + '/0.dige/enzyme_smp.list'
	if os.path.exists(done_file):
		report('INFO', 'The data digestion has been completed, go to the next step')
	else:
		report('INFO', 'Start of data digestion')
		check_dir(O + '/0.dige')
		executor = ProcessPoolExecutor(args.processes)
		pool = []
		with open(enzyme_smp_file, 'w') as OUT:
			for smp, reads in data_dic.items():
				OUT.write('{smp}\t{enzyme_dir}/{smp}/{smp}.fa.gz\n'.format(enzyme_dir = enzyme_dir, smp = smp))
				pool.append(executor.submit(extra_tag, reads, (enzyme_dir + '/' + smp), smp, db_dir))
		executor.shutdown()
		for res in pool:
			res.result()
		exe_shell('touch {}'.format(done_file), 'dige_done')
	return enzyme_smp_file

def run_qual(db_dir, enzyme_smp_file, data_dic, none_micro_smp_file, level):
	threshold = args.threshold
	if threshold.startswith('G'):
		t_m = 'Gscore'
		try:
			t_t = int(threshold.lstrip('G'))
		except:
			report('ERROR', f'The -t {threshold} parameter is wrong, please check')
	elif threshold.startswith('M'):
		t_m = 'ML'
		try:
			t_t = float(threshold.lstrip('M'))
		except:
			report('ERROR', f'The -t {threshold} parameter is wrong, please check')
	else:
		report('ERROR', f'The -t {threshold} parameter is wrong, please check')

	done_file = O + '/1.qual/done_q'
	if os.path.exists(done_file):
		report('INFO', 'The qualitative analysis has been completed, go to the next step')
	else:
		report('INFO', 'Start qualitative analysis')
		check_dir(O + '/1.qual')
#		cc_abd('{}/{}.{}.uniq'.format(db_dir, enzyme, level), '{}/abfh_classify_with_speciename.txt.gz'.format(db_dir), enzyme_smp_file, O + '/1.qual', args.processes, level)
		cc_abd('{db}.{level}.uniq'.format(db = args.database, level = level), '{}/abfh_classify_with_speciename.txt.gz'.format(db_dir), enzyme_smp_file, O + '/1.qual', args.processes, level, 0)
		exe_shell('touch {}'.format(done_file), '')

	done_file = O + '/1.qual/done_p'
	none_micro_smp_list = []
	if os.path.exists(done_file):
		report('INFO', 'The predictive analysis has been completed, go to the next step')
	else:
		report('INFO', 'Start predictive analysis')
		check_dir(O + '/1.qual')
		for smp in data_dic.keys():
			try:
				if t_m == 'Gscore':
					exe_shell('python3 {src_dir}/gscore_filter.py -i {O}/1.qual/{smp}/{smp}.xls -o {O}/1.qual/{smp}/pred.result -g {threshold}'.format(src_dir = src_dir, O = O, smp = smp, threshold = t_t), 'gscore_filter')
				elif t_m == 'ML':
					exe_shell('python3 {src_dir}/MAP2B_ML.py -i {O}/1.qual/{smp}/{smp}.xls -m {config_dir}/XGB_none_0238.pkl -n none -T {threshold} -o {O}/1.qual/{smp}/pred.result'.format(src_dir = src_dir, O = O, smp = smp, config_dir = config_dir, threshold = t_t), 'MAP2B_ML')
				else:
					report('ERROR', f'The -t {threshold} parameter is wrong, please check')
			except:
				report('INFO', 'If an error is reported there, it may be due to the absence of microorganisms, so please ignore it')
				none_micro_smp_list.append(smp)
			tmp_p_list = []
			try:
				with open('{O}/1.qual/{smp}/pred.result'.format(O = O, smp = smp), 'r') as IN:
					for line in IN:
						tmp = line.split()
						if tmp[0] == 'Taxonomy':continue
						if int(tmp[-2]) == 1:
							tmp_p_list.append(1)
			except:
				pass
#			print(tmp_p_list)
			if not tmp_p_list:
				none_micro_smp_list.append(smp)
		if none_micro_smp_list:
			with open(none_micro_smp_file, 'a') as OUT:
				OUT.write('\n'.join(none_micro_smp_list))
		if len(set(none_micro_smp_list)) == len(set(data_dic.keys())):
			report('ERROR', 'No microorganisms were found in all samples')
		exe_shell('touch {}'.format(done_file), 'filter_done')
	return

def get_datadic(enzyme_smp_file, none_micro_smp_file):
	data_dic = check_data(enzyme_smp_file)
	if not os.path.exists(none_micro_smp_file):
		exe_shell('touch {}'.format(none_micro_smp_file), 'prepare')
	else:
		with open(none_micro_smp_file, 'r') as IN:
			for line in IN:
				line = line.strip()
				if line.startswith('#') or not line:continue
				try:
					del data_dic[line.split('\t')[0]]
				except:
					pass
	return data_dic

def run_mkdb(enzyme_fa_dic, db_dir, level):
	done_file = O + '/2.mkdb/done'
	if os.path.exists(done_file):
		report('INFO', 'The quantitative database building has been completed, go to the next step')
	else:
		report('INFO', 'Start quantitative database building')
		check_dir(O + '/2.mkdb')
		executor = ProcessPoolExecutor(args.processes)
		pool = []
		for smp, fa in enzyme_fa_dic.items():
			quan_db = check_dir('{O}/2.mkdb/{smp}'.format(O = O, smp = smp))
			check_dir(quan_db)
			with open('{}/reads.list'.format(quan_db), 'w') as OUT:
				OUT.write('{smp}\t{fa}\n'.format(smp = smp, fa = fa[0]))
			pool.append(executor.submit(mkdb, db_dir, smp, quan_db, args.cutoff, level))
		executor.shutdown()
		for res in pool:
			res.result()
		exe_shell('touch {}'.format(done_file), 'mkdb_done')

def run_quan(data_dic, db_dir, level, cov_thresh):
	done_file = O + '/3.quan/done'
	abd_list = '{}/3.quan/abd.list'.format(O)
	quan_smp_set = set(data_dic.keys())
	if os.path.exists(done_file):
		report('INFO', 'The quantitative analysis has been completed, go to the next step')
	else:
		report('INFO', 'Start quantitative analysis')
		check_dir(O + '/3.quan')
		executor = ProcessPoolExecutor(args.processes)
		pool = []
		for smp in quan_smp_set:
			pool.append(executor.submit(cc_abd, '{O}/2.mkdb/{smp}/{smp}'.format(O = O, smp = smp), '{}/abfh_classify_with_speciename.txt.gz'.format(db_dir), '{}/2.mkdb/{}/reads.list'.format(O, smp), O + '/3.quan', 1, level, cov_thresh))
		executor.shutdown()
		for res in pool:
			res.result()
		with open(abd_list, 'w') as OUT:
			for smp in quan_smp_set: 
				OUT.write('{}\t{}/3.quan/{}/{}.xls\n'.format(smp, O, smp, smp))
		exe_shell('touch {}'.format(done_file), 'quan_done')
	return(abd_list)

def run_stat(abd_list, enzyme_smp_file, db_dir):
	done_file = O + '/4.stat/done'
	if os.path.exists(done_file):
		pass
	else:
		report('INFO', 'Start stat analysis')
		check_dir(O + '/4.stat')
		exe_shell('python3 {src_dir}/MergeProfilesFromMultipleSamples.py -l {abd_list} -o {O}/4.stat/Abundance.tsv'.format(src_dir = src_dir, abd_list = abd_list, O = O), 'MergeProfilesFromMultipleSamples')
		exe_shell('python3 {src_dir}/MergeCoverageFromMultipleSamples.py -i {abd_list} -o {O}/4.stat/Coverage.tsv'.format(src_dir = src_dir, abd_list = abd_list, O = O), 'MergeCoverageFromMultipleSamples')
		exe_shell('python3 {src_dir}/dige_stat.py -i {O}/0.dige -l {enzyme_smp_file} -e {enzyme} -o {O}/4.stat/Tags.tsv'.format(src_dir = src_dir, enzyme = enzyme, enzyme_smp_file = enzyme_smp_file, O = O), 'dige_stat')
		exe_shell('python3 {src_dir}/anno_abund_processor.py -i {O}/4.stat/Abundance.tsv -d {db_dir}/metadata.tsv.gz -o {O}/4.stat'.format(src_dir = src_dir, enzyme = enzyme, db_dir = db_dir, enzyme_smp_file = enzyme_smp_file, O = O, config_dir = config_dir), 'dige_stat')
		exe_shell('touch {}'.format(done_file), 'stat_done')

def main():
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawTextHelpFormatter,
		epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,__mail__,__date__,__version__))
	parser.add_argument('-i',help='The filepath of the sample list. Each line includes an input sample ID and the file path of corresponding DNA sequence data where each field should be separated by <tab>. A line in this file that begins with # will be ignored. like \n \
	sample <tab> shotgun.1.fq(.gz) (<tab> shotgun.2.fq.gz)',dest='input',type=str,required=True)
	parser.add_argument('-o',help='Output directory, default {}/VIP2B_result'.format(os.getcwd()),dest='output',type=str,default='{}/VIP2B_result'.format(os.getcwd()))
	parser.add_argument('-l',help='taxo level, choose from Class/Order/Family/Genus/Species, default Species',dest='level',choices=['Class', 'Order', 'Family', 'Genus', 'Species'],type=str,default='Species')
	parser.add_argument('-e',help='Enzyme, One or more of the following enzymes: AlfI/AloI/BaeI/BcgI/BplI/BsaXI/BslFI/Bsp24I/CjeI/CjePI/CspCI/FalI/HaeIV/Hin4I/PpiI/PsrI. Multiple enzymes should be separated by commas. If you want to select all enzymes, you can simply provide "all". It should be noted that the selection of enzymes needs to correspond with the database. The default combination is AlfI,BcgI,BslFI,CjeI,CjePI,FalI,HaeIV,Hin4I.',dest='enzyme',type=str,default='AlfI,BcgI,BslFI,CjeI,CjePI,FalI,HaeIV,Hin4I')
	parser.add_argument('-d',help=f'Database for pipeline, default {def_db}',dest='database',type=str,default=def_db)
	parser.add_argument('-p',help='Number of processes, note that more threads may require more memory, default 1',dest='processes',type=int,default=1)
	parser.add_argument('-t',help='Threshold for species identification, G5 means using gscore > 5 and M0.5 means using ML probability > 0.5 as a filtering parameter, G2/G5/M0.1/M0.5 are a few commonly used options, default M0.5',dest='threshold',type=str,default='M0.5')
	parser.add_argument('-c',help='cut off for database, default 30000',dest='cutoff',type=int,default=30000)
	parser.add_argument('-f',help='threshold for coverage filtering, default 0.5',dest='cov_thresh',type=float,default=0.5)
	parser.add_argument('--intersection',help='intersection or union of tags between genomes, default union',dest='intersection',action='store_true')
#	parser.add_argument('-n',help='negative control, use comma separation, default is the sample starting with OENC_',dest='nc',type=str,required=False)
#	parser.add_argument('-H',help='host, choose from human/mouse/gallus or use commas to link multiple hosts, default None',dest='host', type=str,default=None)
	global args, enzyme, O
	args=parser.parse_args()
	if args.enzyme == 'all':
		enzyme = 'AlfI,AloI,BaeI,BcgI,BplI,BsaXI,BslFI,Bsp24I,CjeI,CjePI,CspCI,FalI,HaeIV,Hin4I,PpiI,PsrI'
	else:
		enzyme = args.enzyme
#	if args.source == 'SelfBuilt':
#		db_dir = args.database
#	else:
#		db_dir = check_dir(args.database + '/' + args.source)
	db_dir = '/'.join(args.database.split('/')[:-1])

	# check database
#	check_db(db_dir, enzyme, args.source)
	if args.database == def_db:
		check_db(db_dir)
	# prepare data
	O = check_dir(args.output)
	with open(f'{O}/enzyme.log', 'a') as OUT:
		OUT.write('{}\t{}\n'.format(datetime.now().strftime('%Y-%m-%d %H:%M:%S'), enzyme))
	enzyme_smp_file = check_file(args.input)
	none_micro_smp_file = O + '/none_micro_smp.txt'
	data_dic = get_datadic(enzyme_smp_file, none_micro_smp_file)
	# dige
	enzyme_smp_file = run_dige(db_dir, data_dic)
	# qual
	run_qual(db_dir, enzyme_smp_file, data_dic, none_micro_smp_file, args.level)
	# mkdb
	data_dic = get_datadic(enzyme_smp_file, none_micro_smp_file)
	run_mkdb(data_dic, db_dir, args.level)
	# quan
	abd_list = run_quan(data_dic, db_dir, args.level, args.cov_thresh)
	# stat 
	run_stat(abd_list, enzyme_smp_file, db_dir)

	report('INFO', 'Congratulations, all work has been completed')

if __name__=="__main__":
	main()
