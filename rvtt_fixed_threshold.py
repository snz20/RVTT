'''
Author: Sumaiya Nazeen (sumaiya_nazeen@hms.harvard.edu)
This program runs the fixed threshold version of RVTT.
Input:
- mutations file produced by preprocess_gzvcf.py file
- pathway gene list (txt file containing one gene per line)
- tab-separated phenotype file in .fam format 
- fixed minor allele frequency threshold
- number of permutations to run (suggested value = 10000)
- random seed 
- output file name
Output: 
- output file contains the RVTT results for the input genelist under each of the following 
categories of variants: "damaging","damaging_missense","missense_variant","LoF","synonymous","neutral"
'''

import sys
import numpy as np
import scipy as sp
import pandas as pd
import scipy.stats as spstats
from collections import Counter
from joblib import Parallel, delayed

categories = ["structural_interaction_variant","exon_variant","intiator_codon_variant","start_lost","frameshift_variant", "inframe_deletion" ,"inframe_insertion" ,"intron_variant" ,"missense_variant" ,"protein_altering_variant" ,"splice_acceptor_variant" ,"splice_donor_variant" ,"splice_region_variant" ,"stop_gained" ,"stop_lost" ,"synonymous_variant", "damaging", "neutral","High","Medium","Low","PD","PN","SD","SN","damaging_missense","LoF"]
sel_categories = ["damaging","damaging_missense","missense_variant","LoF","synonymous","neutral"]
def calc_test_statistic(features,summary_df):
	obs_z = np.zeros(len(features))
	print("inside_test")
	for i in range(len(features)):
		df = summary_df[i*2:(i+1)*2].iloc[:2,2:]
		df = df.loc[:, (df != 0).any(axis=0)]
		df.index = [0,1]
		df_pct = df/df[df.columns].sum()
		df_cols = df.sum(axis=0)
		df_rows = df.sum(axis=1)
		df_rows_pct = df_rows/df_rows.sum()
		s_bar = 0
		scores = df.columns.map(int)
		for j in range(df.shape[1]):
			s_bar += df_cols.iloc[j]*scores[j]
		N = df.sum().sum()
		s_bar /= N
		denom = 0
		for j in range(df.shape[1]):
			denom += df_cols.iloc[j]*(scores[j] - s_bar)**2
		b = 0
		for j in range(df.shape[1]):	
			b += df_cols.iloc[j]*(scores[j] - s_bar)*(df_pct.iloc[0,j]-df_rows_pct[0])
		b /= denom
		#print(b)
		b_sq = b**2
		z_sq = b_sq / (df_rows_pct[0] * df_rows_pct[1]) * denom
		z = z_sq ** 0.5
		obs_z[i] = z
	return obs_z
	
def calc_p_values(obs_z,features):
	obs_pone = np.ones(len(features))
	obs_ptwo = np.ones(len(features))
	for i in range(len(features)):
		obs_pone[i] = spstats.norm.sf(abs(obs_z[i]))
		obs_ptwo[i] = spstats.norm.sf(abs(obs_z[i]))*2
	return obs_pone, obs_ptwo

def summarize_matrix(df, case_control):
	np_val = df.values
	case_ind = [i for i, x in enumerate(case_control) if x == 2]
	control_ind = [i for i, x in enumerate(case_control) if x == 1]
	nrows = 2 * len(df.columns)
	ncols = 15
	zero_data = np.zeros(shape=(nrows,ncols))
	df2 = pd.DataFrame(zero_data, columns=["Type","Group","0","1","2","3","4","5","6","7","8","9","10","11","12"])
	for i in range(len(df.columns)):
		df2.iloc[2*i,0] = df.columns[i]
		df2.iloc[2*i+1,0] = df.columns[i]
		df2.iloc[2*i,1] = "Case"
		df2.iloc[2*i+1,1] = "Control"
	for i in range(len(sel_categories)):
		case_bins = np.bincount(np_val[case_ind,i])
		control_bins = np.bincount(np_val[control_ind,i])
		case_ii = np.nonzero(case_bins)[0]
		control_ii = np.nonzero(control_bins)[0]
		#print(case_bins,control_bins)
		for a,b in zip(case_ii, case_bins[case_ii]):
			if a<=6:
				df2.iloc[2*i,2+a] = b
			else:
				df2.iloc[2*i,2+7] += b	
		for a,b in zip(control_ii, control_bins[control_ii]):
			if a<=6:
				df2.iloc[2*i+1,2+a] = b
			else:
				df2.iloc[2*i+1,2+7] += b	
	df2 = df2.loc[:, (df2 != 0).any(axis=0)]
	return df2

#def permutation_test(df, case_control, h_maf, N, obs_z,seed):
#	np.random.seed(seed)
#	z_mat = np.zeros((N,len(categories)))
#	maf = h_maf
#	for i in range(N):
#		case_control_p = np.random.permutation(case_control)
#		h_maf = np.random.uniform(min(maf)-0.00001, max(maf),len(maf))
#		z_mat[i] = calc_vt_stat(df, case_control_p, h_maf)
#	p_count = np.zeros(len(categories))
#	for i in range(len(categories)):
#		p_count[i] = len(np.extract(z_mat[:,i]>=obs_z[i],z_mat[:,i]))
#	print(p_count)
#	p_count = (p_count+1)/(N+1)
#	return p_count

def core_task(arglist):
	df = arglist[0]
	names = arglist[1]
	case_control = arglist[2]
	cutoff = arglist[3]
	res = calc_vt_stat(df, names, case_control, cutoff)
	return res

def permutation_test(df, names, case_control, cutoff, N, obs_z,seed):
	np.random.seed(seed)
	#maf = h_maf
	args = []
	for i in range(N):
		case_control_p = np.random.permutation(case_control)
		args.append([df, names, case_control_p, cutoff])
	results = Parallel(n_jobs=20,verbose=0,backend='multiprocessing')(map(delayed(core_task), args))
	z_mat = np.array(results)
	p_count = np.zeros(len(sel_categories))
	for i in range(len(sel_categories)):
		p_count[i] = len(np.extract(z_mat[:,i]>=obs_z[i],z_mat[:,i]))
	print(z_mat)
	print(p_count)
	p_count = (p_count+1)/(N+1)
	print(p_count)
	return p_count
	
def read_files(matfile, genefile, famfile, cutoff):
	dfall = pd.read_table(matfile,sep='\t',header=0)
	d = dict()
	for k in dfall.columns:
		d[k] = k.strip('#')
	dfall = dfall.rename(columns=d)
	genes = [l.strip() for l in open(genefile)]
	df = dfall[dfall.Gene.isin(genes)]
	lines = [l.strip() for l in open(famfile)]
	case_control = []
	names = []
	for l in lines[1:]:
		x = l.split('\t')
		ind = x[1]
		val = int(x[5])
		#if x[1] in names:
		names.append(ind)
		case_control.append(val)
	return df, names, case_control

def create_indiv_count_matrix(df, names, case_control, cutoff):
	df_sel = df.loc[df['PopMAF'] <= cutoff]
	hdr = list(df_sel.columns)
	pind = hdr.index('polyphen')
	sind = hdr.index('sift')
	csqind = hdr.index('top_csq')
	csq2ind = hdr.index('csq2')
	impind = hdr.index('impact')
	afind = hdr.index('PopMAF')
	indv_ind = hdr.index('mutated_individuals')
	r,c = df_sel.shape   
	df_path = pd.DataFrame(data=np.zeros((len(names),len(sel_categories))),index=names,columns=sel_categories,dtype=int)	
	for i in range(r):
		cur_row = df_sel.iloc[i,:]
		gid = cur_row[0]
		mutid = cur_row[0]+'_'+cur_row[1]+'_'+cur_row[3]
		tot_len = len(cur_row)
		polyphen = cur_row[pind]
		sift = cur_row[sind]
		popmaf = float(cur_row[afind])
		category = cur_row[csqind]
		impact = cur_row[impind]
		csq = cur_row[csqind]
		csq2 = cur_row[csq2ind]
		csq3 = '.'		
		if csq2 == 'damaging':
			if csq == "frameshift_variant" or csq=="splice_acceptor_variant" or csq=="splice_donor_variant" or csq=="stop_gained" or csq=="stop_lost" or csq=="start_lost" or csq=="splice_region_variant" or csq=="structural_interaction_variant" or csq=="initiator_codon_variant":
				csq3  = 'LoF'
			elif csq == "missense_variant":
				csq3 = 'damaging_missense'
		gts = cur_row[indv_ind].split(';')
		for v in gts:
			i = -1
			if v in names:
				i = names.index(v)
			if i!=-1:
				if csq in sel_categories:
					df_path.loc[v,csq] += 1
				if csq2 in sel_categories:
					df_path.loc[v,csq2] += 1
				if csq3 in sel_categories:
					df_path.loc[v,csq3] += 1
	return df_path
	
def calc_vt_stat(df, names, case_control, cutoff):
	df_path = create_indiv_count_matrix(df, names, case_control, cutoff)
	summary_df = summarize_matrix(df_path, case_control)
    	#print(summary_df)
	features = df_path.columns
	obs_z = calc_test_statistic(features,summary_df)
	#print(obs_z)
	z_score = obs_z
	z_score_mod = np.nan_to_num(z_score)
	return z_score_mod
		
def main():
    if len(sys.argv) < 8:
        print("Usage: python rvtt_fixed_threshold.py <variantfile> <genelist> <famfile> <cutoff> <N> <seed> <outfile>")
        exit(1)	
    variantfile = sys.argv[1]
    genelist = sys.argv[2]
    famfile = sys.argv[3]
    cutoff = float(sys.argv[4])
    N = int(sys.argv[5])
    seed = int(sys.argv[6])
    outfile = sys.argv[7]
    df, names, case_control= read_files(variantfile, genelist, famfile, cutoff)
    obs_z = calc_vt_stat(df, names, case_control, cutoff)
    perm_p = permutation_test(df, names, case_control, cutoff, N, obs_z,seed)
    resc = pd.DataFrame({'Type': sel_categories, 'vt-fixed': cutoff, 'vt-z-score': list(obs_z), 'perm-p': list(perm_p)}, columns=['Type', 'vt-fixed', 'vt-z-score','perm-p'])
	#print(obs_hmaf, obs_z)
    print(resc)
    resc.to_csv(outfile, sep='\t',index=False)
	
if __name__=="__main__":
	main()
