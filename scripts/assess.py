#!/data/software/install/miniconda3/envs/python.3.7.0/bin/python3
########################################## import ################################################
import argparse, os, sys
import pandas as pd
from datetime import datetime
from scipy.spatial import distance
############################################ ___ #################################################
__doc__ = ''
__author__ = 'Liu Jiang'
__mail__ = 'jiang.liu@oebiotech.com'
__date__ = '2024/09/30 14:40:24'
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
	parser.add_argument('-g',help='ground truth file',dest='ground_truth',type=str,required=True)
	parser.add_argument('-i',help='abundance file',dest='abundance',type=str,required=True)
	parser.add_argument('-t',help='tag',dest='tag',type=str,required=False)
	parser.add_argument('-o',help='output file',dest='output',type=str,required=True)
	args=parser.parse_args()
	info = "runing..."
	report("INFO",info)

	if args.tag:
		tag = f'{args.tag}\t'
	else:
		tag = ''
	gt_df = pd.read_csv(check_file(args.ground_truth), sep = '\t', index_col = 0)
	smp_df = pd.read_csv(check_file(args.abundance), sep = '\t', index_col = 0)
	with open(args.output, 'w') as OUT:
		if tag:
			OUT.write(f'Tag\tSample\tTruth\tDetected\tTP\tFP\tFN\tPrecision\tRecall\tF1\tL2\tBC\n')
		else:
			OUT.write(f'Sample\tTruth\tDetected\tTP\tFP\tFN\tPrecision\tRecall\tF1\tL2\tBC\n')
		for smp in gt_df.columns:
			if smp in smp_df:
				gt_df_tmp = gt_df[smp][gt_df[smp] > 0]
				smp_df_tmp = smp_df[smp][smp_df[smp] > 0]
				TP = len(smp_df_tmp[smp_df_tmp.index.isin(gt_df_tmp.index)])
				FP = len(smp_df_tmp[~smp_df_tmp.index.isin(gt_df_tmp.index)])
				FN = len(gt_df_tmp[~gt_df_tmp.index.isin(smp_df_tmp.index)])
				GT = len(gt_df_tmp)
				DT = len(smp_df_tmp)
				merge_df_tmp = pd.merge(gt_df_tmp, smp_df_tmp, left_index=True, right_index=True, how='outer').fillna(0)
				col1 = merge_df_tmp[smp + '_x'].values
				col2 = merge_df_tmp[smp + '_y'].values
				try:
					L2 = round(1 - distance.euclidean(col1, col2), 4)
				except:
					La = 'nan'
				try:
					BC = round(1 - distance.braycurtis(col1, col2), 4)
				except:
					BC = 'nan'
				try:
					precision = round(TP / DT, 4)
				except:
					precision = 'nan'
				try:
					recall = round(TP / GT, 4)
				except:
					recall = 'nan'
				try:
					F1 = round(2 * ((precision * recall)/(precision + recall)), 4)
				except:
					F1 = 'nan'
				OUT.write(f'{tag}{smp}\t{GT}\t{DT}\t{TP}\t{FP}\t{FN}\t{precision}\t{recall}\t{F1}\t{L2}\t{BC}\n')
			else:
				report('INFO', f'{smp} not in ground truth file')

if __name__=="__main__":
	main()
	info = "finish!"
	report("INFO",info)
