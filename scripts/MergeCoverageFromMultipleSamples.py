import argparse, os, sys, re, random, glob
from datetime import datetime
bindir = os.path.abspath(os.path.dirname(__file__))

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
	return()

def check_file(file):
	if os.path.exists(file):
		return(os.path.abspath(file))
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

def main():
	parser=argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
	parser.add_argument('-i',help='input file',dest='input',type=str,required=True)
	parser.add_argument('-o',help='output file',dest='output',type=str,required=True)
	parser.add_argument('-f',help='filter Taxonomy by last comma-separated part (default: human)',dest='filter',type=str,default="human")
	args=parser.parse_args()
	info = "running..."
#	report("INFO",info)
	info_dic = {}
	I = check_file(args.input)
	with open(I,'r') as IN:
		for line1 in IN:
			smp, abd = line1.split()
			with open(abd, 'r') as ABD:
				for line2 in ABD:
					if line2.startswith('Taxonomy'):continue
					tmp = line2.strip().split('\t')
					tmp_dic = info_dic.setdefault(smp, {})
					tmp_dic.setdefault(tmp[0], tmp[3])
	
	# 收集所有物种并根据过滤条件筛选
	all_spe_list = []
	for smp, spe_cov_dic in info_dic.items():
		all_spe_list.extend(spe_cov_dic.keys())
	all_spe_set = set(all_spe_list)
	
	# 应用过滤：移除最后一部分匹配的物种
	if args.filter:
		filtered_spe = []
		for spe in all_spe_set:
			parts = spe.split(',')
			# 检查是否有逗号分隔部分，且最后一部分匹配过滤条件
			if parts and parts[-1] == args.filter:
				continue  # 跳过需要过滤的物种
			filtered_spe.append(spe)
		all_spe_sorted = sorted(filtered_spe)
	else:
		all_spe_sorted = sorted(all_spe_set)

	smp_list = list(info_dic.keys())
	with open(args.output, 'w') as OUT:
		OUT.write('Taxonomy\t{}\n'.format('\t'.join(smp_list)))
		for spe in all_spe_sorted:
			cov_list = []
			for smp in smp_list:
				try:
					# 保持原脚本的丰度处理逻辑（>1则设为1，否则用原始值）
					if float(info_dic[smp][spe]) > 1:
						cov = '1'
					else:
						cov = info_dic[smp][spe]
				except KeyError:
					cov_list.append('0')
			OUT.write('{}\t{}\n'.format(spe, '\t'.join(cov_list)))

if __name__=="__main__":
	main()
	info = "finish!"
#	report("INFO",info)
