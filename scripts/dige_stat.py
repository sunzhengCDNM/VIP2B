#!/data/software/install/miniconda3/envs/python.3.7.0/bin/python3
########################################## import ################################################
import argparse, os, sys, re, random, glob
from datetime import datetime
bindir = os.path.abspath(os.path.dirname(__file__))
sys.path.append('/data/USER/liujiang/script/lib')
#import pandas as pd
############################################ ___ #################################################
__doc__ = ''
__author__ = 'Liu Jiang'
__mail__ = 'jiang.liu@oebiotech.com'
__date__ = '2023/06/24 14:03:32'
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
	parser=argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.RawTextHelpFormatter,
		epilog='author:\t{0}\nmail:\t{1}\ndate:\t{2}\nversion:\t{3}'.format(__author__,__mail__,__date__,__version__))
	parser.add_argument('-i',help='input dir',dest='input',type=str,required=True)
	parser.add_argument('-l',help='list',dest='list',type=str,required=True)
	parser.add_argument('-e',help='enzyme',dest='enzyme',type=str,required=True)
	parser.add_argument('-o',help='output file',dest='output',type=str,required=True)
	args=parser.parse_args()
	I = check_file(args.list)
	with open(I,'r') as IN, open(args.output, 'w') as OUT:
#		OUT.write('{}\n'.format('\t'.join('Sample,Enzyme,Raw_Reads,Enzyme_Tag,Enzyme_Tag_Ratio,Host_Tag,Host_Tag_Ratio,Other_Tag,Other_Tag_Ratio'.split(','))))
		OUT.write('{}\n'.format('\t'.join('Sample,Enzyme,Raw_Reads,Enzyme_Tag,Enzyme_Tag_Ratio'.split(','))))
		for line in IN:
			smp = line.strip().split('\t')[0]
			dige_stat = glob.glob('{}/{}/{}.dige.stat.xls'.format(args.input, smp, smp))
#			host_stat = glob.glob('{}/{}/{}.{}.host.stat.xls'.format(args.input, smp, smp, args.enzyme))
			with open(dige_stat[0], 'r') as DIGE:
				for line in DIGE:
					tmp = line.strip().split('\t')
					if tmp[0] == 'sample':continue
					smp_dige_stat = tmp
#			if host_stat:
#				with open(host_stat[0], 'r') as HOST:
#					for line in HOST:
#						tmp = line.strip().split('\t')
#						if tmp[0] != 'final':continue
#						OUT.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format('\t'.join(smp_dige_stat[:-1]), float((smp_dige_stat[-1]).strip('%')) / 100,  tmp[2], round((int(tmp[2]) / int(smp_dige_stat[3])), 4), tmp[3], round((int(tmp[3]) / int(smp_dige_stat[3])), 4)))
#			else:
#				OUT.write('{}\t{}\n'.format('\t'.join(smp_dige_stat[:-1]), float((smp_dige_stat[-1]).strip('%')) / 100, smp_dige_stat[3]))
				OUT.write('{}\t{}\n'.format('\t'.join(smp_dige_stat[:-1]), round((int(smp_dige_stat[3])/int(smp_dige_stat[2])), 4)))

if __name__=="__main__":
	main()
