#!/bin/python

import sys
import numpy as np

floatpoint_error = 1e-15

def ReadSalmonQuant(filename):
	NumReads = {}
	Expression = {}
	fp = open(filename, 'r')
	linecount = 0
	for line in fp:
		linecount += 1
		if linecount == 1:
			continue
		strs = line.strip().split("\t")
		NumReads[strs[0]] = float(strs[4])
		Expression[strs[0]] = float(strs[4]) / float(strs[2])
	fp.close()
	return NumReads, Expression


def ReadAbundanceGap(filename):
	AbundanceGap = {}
	fp = open(filename, 'r')
	for line in fp:
		if line[0] == '#':
			continue
		strs = line.strip().split("\t")
		lb = float(strs[1])
		ub = float(strs[2])
		if lb < 0:
			lb = 0
		if ub < 0:
			ub = 0
		if lb > ub:
			ub = lb
		AbundanceGap[strs[0]] = (lb, ub)
	fp.close()
	return AbundanceGap


def NormalizeAbundanceGap(AbundanceGap, NumReads):
	'''
	normalization of salmon flow:
		the flow of salmon is calculated by transcript numreads / effective length (sum of flow * efflen over all transcripts = numreads)
	normalization of graph salmon flow:
		sum of (flow * efflen) = numreads
	different samples have different numreads
	post-normalize lb and ub is equivalent to pre-normalize the flow and bound using the normalized flows
	normalize by the numreads so that different samples have the same number of numreads
	'''
	sumreads = np.sum(list(NumReads.values()))
	AbundanceGap = {t:(v[0] / sumreads * 1e6, v[1] / sumreads * 1e6) for t,v in AbundanceGap.items()}
	return AbundanceGap


def NormalizeExpression(Expression, NumReads):
	sumreads = np.sum(list(NumReads.values()))
	Expression = {t:v / sumreads * 1e6 for t,v in Expression.items()}
	return Expression


def ReadCelltypeMetadata(filename):
	'''
	this is for Immune29 dataset
	'''
	CelltypeSampleMap = {}
	fp = open(filename, 'r')
	for line in fp:
		strs = line.strip().split(",")
		if strs[30] in CelltypeSampleMap:
			CelltypeSampleMap[strs[30]].append(strs[0])
		else:
			CelltypeSampleMap[strs[30]] = [strs[0]]
	fp.close()
	return CelltypeSampleMap


def ProcessUncertaintyMatrix(GroupSampleMap, folder, id_prefix = "", bound_suffix = "graphsalmon/gs_maxflow_bound.txt"):
	transnames = None
	Expression_pergroup = {}
	Lb_pergroup = {}
	Ub_pergroup = {}
	for c,samples in GroupSampleMap.items():
		mat_exp = None
		mat_lb = None
		mat_ub = None
		for s in samples:
			salmonfile = folder + "/" + id_prefix + s + "/quant.sf"
			boundfile = folder + "/" + id_prefix + s + "/" + bound_suffix
			NumReads, Expression = ReadSalmonQuant(salmonfile)
			Expression = NormalizeExpression(Expression, NumReads)
			Bounds = ReadAbundanceGap(boundfile)
			Bounds = NormalizeAbundanceGap(Bounds, NumReads)
			# Bounds = AdjustBounds(Bounds, Expression)
			print("finish reading {}:{}".format(c, s))
			# order the transcript names by dictionary order
			tmptransnames = list(Expression.keys())
			tmptransnames.sort()
			if transnames is None:
				transnames = tmptransnames
			else:
				assert( len(transnames) == len(tmptransnames) and np.all([transnames[i] == tmptransnames[i] for i in range(len(transnames))]) )
			# process matrix
			vec_exp = np.array([Expression[t] for t in transnames]).reshape( (len(transnames), 1) )
			vec_lb = np.array([Bounds[t][0] for t in transnames]).reshape( (len(transnames), 1) )
			vec_ub = np.array([Bounds[t][1] for t in transnames]).reshape( (len(transnames), 1) )
			if mat_exp is None:
				mat_exp = vec_exp
				mat_lb = vec_lb
				mat_ub = vec_ub
			else:
				mat_exp = np.hstack( (mat_exp, vec_exp) )
				mat_lb = np.hstack( (mat_lb, vec_lb) )
				mat_ub = np.hstack( (mat_ub, vec_ub) )
		Expression_pergroup[c] = mat_exp
		Lb_pergroup[c] = mat_lb
		Ub_pergroup[c] = mat_ub
	return transnames, Expression_pergroup, Lb_pergroup, Ub_pergroup


def MeanLowerBounds_varylambda(salmon_exps, lbs):
	assert(len(salmon_exps) == len(lbs))
	interception = np.mean(salmon_exps)
	slope = np.mean([lbs[i] - salmon_exps[i] for i in range(len(lbs))])
	return interception, slope


def MeanUpperBounds_varylambda(salmon_exps, ubs):
	assert(len(salmon_exps) == len(ubs))
	interception = np.mean(salmon_exps)
	slope = np.mean([ubs[i] - salmon_exps[i] for i in range(len(ubs))])
	return interception, slope


def CalculateIValue_mean(interception1_min, slope1_min, interception2_min, slope2_min, interception1_max, slope1_max, interception2_max, slope2_max, quantile = 0):
	# lb and ub at the region start
	lb1_start = interception1_min
	ub1_start = interception1_max
	lb2_start = interception2_min
	ub2_start = interception2_max
	# lb and ub at the region end
	lb1_end = interception1_min + slope1_min
	ub1_end = interception1_max + slope1_max
	lb2_end = interception2_min + slope2_min
	ub2_end = interception2_max + slope2_max
	# adjust for quantile for start
	q1_start = (ub1_start - lb1_start) * quantile
	q2_start = (ub2_start - lb2_start) * quantile
	lb1_start += q1_start
	ub1_start -= q1_start
	assert(lb1_start <= ub1_start)
	lb2_start += q2_start
	ub2_start -= q2_start
	assert(lb2_start <= ub2_start)
	# adjust for quantile for end
	q1_end = (ub1_end - lb1_end) * quantile
	q2_end = (ub2_end - lb2_end) * quantile
	lb1_end += q1_end
	ub1_end -= q1_end
	assert(lb1_end <= ub1_end)
	lb2_end += q2_end
	ub2_end -= q2_end
	assert(lb2_end <= ub2_end)
	# adjust for quantile for slope
	q_slope1_min = lb1_end - lb1_start
	q_slope1_max = ub1_end - ub1_start
	q_slope2_min = lb2_end - lb2_start
	q_slope2_max = ub2_end - ub2_start
	# calculate IV
	IV = 1
	# overlap in the beginning
	if max(lb1_start, lb2_start) <= min(ub1_start, ub2_start):
		IV = 0
	elif ub1_start < lb2_start:
		# not overlap at all through all lambda values
		if ub1_end < lb2_end:
			IV = 2
		# overlap in the middle
		else:
			x = (lb2_start - ub1_start) / (q_slope1_max - q_slope2_min)
			assert(x >= 0 and x <= 1)
			IV = x
	elif ub2_start < lb1_start:
		# not overlap at all through all lambda values
		if ub2_end < lb1_end:
			IV = 2
		# overlap in the middle
		else:
			x = (ub2_start - lb1_start) / (q_slope1_min - q_slope2_max)
			assert(x >= 0 and x <= 1)
			IV = x
	return IV


def WriteIV_mean(outputfile, transnames, Expression_pergroup, Lb_pergroup, Ub_pergroup):
	assert(len(Expression_pergroup) == 2 and len(Lb_pergroup) == 2 and len(Ub_pergroup) == 2)
	groups = list(Expression_pergroup.keys())
	groups.sort()
	# process and write IV for each transcript
	fp = open(outputfile, 'w')
	fp.write("# Name\t{}Mean\t{}Mean\tIV\tIV25\n".format(groups[0], groups[1]))
	for i in range(len(transnames)):
		t = transnames[i]
		c1_exp = Expression_pergroup[groups[0]][i,:]
		c1_lb = Lb_pergroup[groups[0]][i,:]
		c1_ub = Ub_pergroup[groups[0]][i,:]
		c2_exp = Expression_pergroup[groups[1]][i,:]
		c2_lb = Lb_pergroup[groups[1]][i,:]
		c2_ub = Ub_pergroup[groups[1]][i,:]
		# min lb of c1
		interception1_min, slope1_min = MeanLowerBounds_varylambda(c1_exp, c1_lb)
		# max ub for c1
		interception1_max, slope1_max = MeanUpperBounds_varylambda(c1_exp, c1_ub)
		# min lb for c2
		interception2_min, slope2_min = MeanLowerBounds_varylambda(c2_exp, c2_lb)
		# max ub for c2
		interception2_max, slope2_max = MeanUpperBounds_varylambda(c2_exp, c2_ub)
		# calculate IV
		IV = CalculateIValue_mean(interception1_min, slope1_min, interception2_min, slope2_min, interception1_max, slope1_max, interception2_max, slope2_max)
		IV25 = CalculateIValue_mean(interception1_min, slope1_min, interception2_min, slope2_min, interception1_max, slope1_max, interception2_max, slope2_max, quantile = 0.25)
		fp.write("{}\t{}\t{}\t{}\t{}\n".format(t, np.mean(c1_exp), np.mean(c2_exp), IV, IV25))
	fp.close()


def WriteExampleCurve(outputfile, indexes, GroupSampleMap, transnames, Expression_pergroup, Lb_pergroup, Ub_pergroup):
	fp = open(outputfile, 'w')
	fp.write("# Name\tGroup\tSample\treference_proportion\tlb\tub\n")
	for idx in indexes:
		for c,samples in GroupSampleMap.items():
			tname = transnames[idx]
			exp = Expression_pergroup[c][idx,:]
			lb = Lb_pergroup[c][idx,:]
			ub = Ub_pergroup[c][idx,:]
			assert( len(exp) == len(samples) )
			# writing the begining and ending of line of individual samples
			for i in range(len(samples)):
				s = samples[i]
				fp.write("{}\t{}\t{}\t0\t{}\t{}\n".format(tname, c, s, exp[i], exp[i]))
				fp.write("{}\t{}\t{}\t1\t{}\t{}\n".format(tname, c, s, lb[i], ub[i]))
			# writing the line of the mean of the group
			fp.write("{}\t{}\t{}\t0\t{}\t{}\n".format(tname, c, "Mean", np.mean(exp), np.mean(exp)))
			fp.write("{}\t{}\t{}\t1\t{}\t{}\n".format(tname, c, "Mean", np.mean(lb), np.mean(ub)))
	fp.close()


def tmp():
	print(Expression_pergroup[groups[0]][i,:])
	print(Expression_pergroup[groups[1]][i,:])
	print(Lb_pergroup[groups[0]][i,:])
	print(Lb_pergroup[groups[1]][i,:])
	print(Ub_pergroup[groups[0]][i,:])
	print(Ub_pergroup[groups[1]][i,:])


if __name__=="__main__":
	if len(sys.argv) == 1:
		print("python ComputeIValue.py <output folder as in metarun.sh>")
	else:
		script_folder = sys.argv[0]
		folder = sys.argv[1]

		# get the meta file directory from script_folder
		script_folder_split = script_folder.split("/")
		if len(script_folder_split) > 2:
			script_folder = "/".join(script_folder_split[:-2])
		elif len(script_folder_split) == 2:
			script_folder = "./"
		elif len(script_folder_split) == 1:
			script_folder = "../"
		metafile = script_folder + "data/Metadata_immune29.txt"
		print(metafile)

		CelltypeSampleMap = ReadCelltypeMetadata(metafile)
		print(CelltypeSampleMap)

		transnames, Expression_pergroup, Lb_pergroup, Ub_pergroup = ProcessUncertaintyMatrix(CelltypeSampleMap, folder + "/Immune29/", "", "graphsalmon/gs_maxflow_bound.txt")
		WriteIV_mean(folder + "/Immune29/IValue_mean_ext.txt", transnames, Expression_pergroup, Lb_pergroup, Ub_pergroup)

		# most unreliable ones, suceptible to <10% unannotated transcript expression
		tnames = ["ENST00000424085.6", "ENST00000258412.7", "ENST00000445061.5", "ENST00000273258.3", "ENST00000296255.7", "ENST00000503065.1",\
		 "ENST00000252725.9", "ENST00000339121.9", "ENST00000396466.5", "ENST00000394936.7", "ENST00000531791.1", "ENST00000545638.2", "ENST00000439220.6",\
		 "ENST00000396410.8", "ENST00000564400.5", "ENST00000420290.6", "ENST00000359680.9", "ENST00000476532.1", "ENST00000369762.6"]
		# reliable ones, the comparison between conditions using salmon expression and ub agrees with each other
		tnames += ["ENST00000619423.4", "ENST00000435064.5", "ENST00000371733.7"]
		# reliable oens, the comparison between conditions using salmon and ub are opposite
		tnames += ["ENST00000314289.12"]
		indexes = [transnames.index(tname) for tname in tnames]
		WriteExampleCurve(folder + "/Immune29/IValue_curve_example.txt", indexes, CelltypeSampleMap, transnames, Expression_pergroup, Lb_pergroup, Ub_pergroup)
