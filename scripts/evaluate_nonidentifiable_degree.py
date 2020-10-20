#!/bin/python

import sys
from pathlib import Path
import numpy as np
from matplotlib import pyplot as plt
import seaborn
import pandas as pd
from matplotlib_venn import venn2


def GetFeature(line, key):
	s=line.index(key)
	t=line.index(";", s+1)
	return line[(s+len(key)+2):(t-1)]


def map_gene_trans(gtffile):
	GeneTransMap = {}
	TransGeneMap = {}
	with open(gtffile, 'r') as fp:
		for line in fp:
			if line[0] == '#':
				continue
			strs = line.strip().split("\t")
			if strs[2] == "transcript":
				gene_id = GetFeature(line, "gene_id")
				trans_id = GetFeature(line, "transcript_id")
				if gene_id in GeneTransMap:
					GeneTransMap[gene_id].append( trans_id )
				else:
					GeneTransMap[gene_id] = [trans_id]
				TransGeneMap[trans_id] = gene_id
	return GeneTransMap, TransGeneMap


def read_both_bounds_without_normalization(filename, eps = 1e-8):
	exp_bounds = {}
	with open(filename, 'r') as fp:
		for line in fp:
			if line[0] == '#':
				continue
			strs = line.strip().split("\t")
			lb = max(0, float(strs[1]))
			ub = max(0, float(strs[2]))
			if lb > ub - eps:
				ub = lb
			exp_bounds[strs[0]] = (lb, ub)
	return exp_bounds


def read_bootstrap_map(filename):
	bootstrap_map = {}
	tmp_mat = np.loadtxt(filename, skiprows=1)
	trans_names = []
	with open(filename, 'r') as fp:
		line = fp.readline()
		trans_names = line.strip().split("\t")
	assert( len(trans_names) == tmp_mat.shape[1] )
	for i in range(len(trans_names)):
		bootstrap_map[trans_names[i]] = tmp_mat[:, i]
	return bootstrap_map


def read_express_solvable_flag(filename):
	unsolvable_transcripts = []
	with open(filename, 'r') as fp:
		linecount = 0;
		for line in fp:
			linecount += 1
			if linecount == 1:
				continue
			strs = line.strip().split("\t")
			assert( strs[13] == "F" or strs[13] == "T" )
			if strs[13] == "F":
				unsolvable_transcripts.append( strs[1] )
	return unsolvable_transcripts


def count_isoforms_with_permuted_ranking_bounds(exp_bounds, GeneTransMap):
	total_list = []
	for g,v in GeneTransMap.items():
		trans_list = [x for x in v if x in exp_bounds]
		reverted_trans_list = []
		for idx1 in range(len(trans_list)):
			for idx2 in range(idx1+1, len(trans_list)):
				t1_range = exp_bounds[trans_list[idx1]]
				t2_range = exp_bounds[trans_list[idx2]]
				# check whether overlap
				if max(t1_range[0], t2_range[0]) < min(t1_range[1], t2_range[1]):
					reverted_trans_list += [trans_list[idx1], trans_list[idx2]]
		reverted_trans_list = list(set(reverted_trans_list))
		total_list += reverted_trans_list
	return total_list


def count_isoforms_with_permuted_ranking_bootstrap(bootstrap_map, GeneTransMap):
	total_list = []
	for g,v in GeneTransMap.items():
		trans_list = [x for x in v if x in bootstrap_map]
		reverted_trans_list = []
		for idx1 in range(len(trans_list)):
			for idx2 in range(idx1+1, len(trans_list)):
				t1_bootstrap = bootstrap_map[trans_list[idx1]]
				t2_bootstrap = bootstrap_map[trans_list[idx2]]
				if np.sum(t1_bootstrap > t2_bootstrap) > 0 and np.sum(t2_bootstrap > t1_bootstrap) > 0:
					reverted_trans_list += [trans_list[idx1], trans_list[idx2]]
		reverted_trans_list = list(set(reverted_trans_list))
		total_list += reverted_trans_list
	return total_list


def process_one_sample(salmon_folder, GeneTransMap):
	# read bootstrap matrix
	bootstrap_map = read_bootstrap_map(salmon_folder + "/quant_bootstraps.tsv")
	# read uncertain bounds under complete reference assumption
	lp_exp_bounds = read_both_bounds_without_normalization(salmon_folder + "/graphsalmon/salmon_lp_bound.txt")
	maxflow_exp_bounds = read_both_bounds_without_normalization(salmon_folder + "/graphsalmon/gs_maxflow_bound.txt")
	assert(set(lp_exp_bounds.keys()) == set(maxflow_exp_bounds.keys()))
	# a extra info
	num_total_trans = np.sum([len(v) for v in GeneTransMap.values()])
	sample_id = salmon_folder.split("/")[-1].split("_")[-1]
	# for each lambda values, compute the number of isoforms such that its ranking can be permuted with another isoform if uncertainty is considered
	df = None
	for lambd in np.arange(0, 1.01, 0.1):
		tmp_exp_bounds = {k:(lambd*maxflow_exp_bounds[k][0] + (1-lambd)*v[0], lambd*maxflow_exp_bounds[k][1] + (1-lambd)*v[1]) for k,v in lp_exp_bounds.items()}
		trans_list_bound = count_isoforms_with_permuted_ranking_bounds(tmp_exp_bounds, GeneTransMap)
		count_bound = len(trans_list_bound)
		tmp = pd.DataFrame( {"sample id":sample_id, "lambda":lambd, "Percent transcripts":1.0*count_bound/num_total_trans, "source":"non-identifiability"}, index=[0] )
		if df is None:
			df = tmp
		else:
			df = df.append( tmp )
	trans_list_bootstrap = count_isoforms_with_permuted_ranking_bootstrap(bootstrap_map, GeneTransMap)
	count_bootstrap = len(trans_list_bootstrap)
	df = df.append( pd.DataFrame( {"sample id":sample_id, "lambda":2, "Percent transcripts":1.0*count_bootstrap/num_total_trans, "source":"bootstrapping"}, index=[0] ) )
	print( sample_id, df.shape )
	return df


def plot_percentage_flipped_transcripts(hubmap_path, data_path):
	# folders = [str(x) for x in Path(hubmap_path).glob("salmon*")]
	folders = [str(x) for x in Path(hubmap_path).glob("ERR*")]
	folders.sort()
	GeneTransMap, TransGeneMap = map_gene_trans(data_path + "/gencode.v26.annotation.gtf")
	df = None
	for salmon_folder in folders:
		if df is None:
			df = process_one_sample(salmon_folder, GeneTransMap)
		else:
			df = df.append( process_one_sample(salmon_folder, GeneTransMap) )
	# plotting
	plt.rcParams.update({'font.size': 14, 'text.usetex': True})
	fig, axes = plt.subplots(1, 2, figsize = (8,3.5), gridspec_kw={'width_ratios': [1.5, 1]})
	seaborn.boxplot(data=df, x="lambda", y="Percent transcripts", linewidth=1, ax=axes[0])
	axes[0].set_xticklabels(["{:.2f}".format(x) for x in np.arange(0, 1.01, 0.1)] + ["bootstrapping"], rotation=30)
	axes[0].set_xlabel(r'\lambda')
	# venn diagram to compare the overlapping flipped isoforms under lambda = 0.9 and bootstrapping
	bootstrap_map = read_bootstrap_map(folders[5] + "/quant_bootstraps.tsv")
	lp_exp_bounds = read_both_bounds_without_normalization(folders[5] + "/graphsalmon/salmon_lp_bound.txt")
	maxflow_exp_bounds = read_both_bounds_without_normalization(folders[5] + "/graphsalmon/gs_maxflow_bound.txt")
	tmp_exp_bounds = {k:(0.7*maxflow_exp_bounds[k][0] + 0.3*v[0], 0.7*maxflow_exp_bounds[k][1] + 0.3*v[1]) for k,v in lp_exp_bounds.items()}
	trans_list_bound = set(count_isoforms_with_permuted_ranking_bounds(tmp_exp_bounds, GeneTransMap))
	trans_list_bootstrap = set(count_isoforms_with_permuted_ranking_bootstrap(bootstrap_map, GeneTransMap))
	plt.rcParams.update({'font.size': 12, 'text.usetex': True})
	v = venn2(subsets = {'10':len(trans_list_bound-trans_list_bootstrap), '01':len(trans_list_bootstrap-trans_list_bound), '11':len(trans_list_bound & trans_list_bootstrap)}, set_labels=('', ''), ax=axes[1])
	v.get_patch_by_id('10').set_color('yellow')
	v.get_patch_by_id('01').set_color('red')
	v.get_patch_by_id('11').set_color('orange')
	v.get_label_by_id('10').set_text("non-identifiability\n{}".format( len(trans_list_bound-trans_list_bootstrap) ))
	v.get_label_by_id('10').set_y(0.1)
	v.get_label_by_id('01').set_text("bootstrapping\n{}".format( len(trans_list_bootstrap-trans_list_bound) ))
	v.get_label_by_id('01').set_y(0.1)
	v.get_label_by_id('11').set_text("shared\n{}".format( len(trans_list_bound & trans_list_bootstrap) ))
	v.get_label_by_id('11').set_y(-0.1)
	fig.tight_layout()
	fig.savefig(hubmap_path + "/percentage_flipped_transcripts.pdf", transparent = True, bbox_inches='tight')


if __name__ == "__main__":
	if len(sys.argv) == 1:
		print("python evaluate_nonidentifiable_degree.py <path to Human Body Map quantification output> <path to data downloading output>")
	else:
		hubmap_path = sys.argv[1]
		data_path = sys.argv[2]
		plot_percentage_flipped_transcripts()