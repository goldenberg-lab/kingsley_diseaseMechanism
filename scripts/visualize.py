import numpy as np
import pandas as pd
from scipy import stats, integrate
import matplotlib.pyplot as plt
import time
import sys
import pyBigWig
import os

import seaborn as sns
sns.set(color_codes=True)

def	drawDistanceDist():

	distance_list = []
	with open('data/sorted_enhancers.txt') as fp:
		firstline = fp.readline().strip()
		
		for line in fp:
			tmp = line.strip().split()

			chromosome = tmp[0]
			start = int(tmp[1])
			end = int(tmp[2])

			distance_list.append(end-start)

			if end-start > 100000:
				print('enhancer: %s, %d, %d, %d' % (chromosome, start, end, end-start))

	drawFigure(distance_list, 'enhancer_distance', 20000, 0, 50000)

def	drawInterEnhancersDist():
	distance_list = []
	with open('data/sorted_enhancers.txt') as fp:
		# Ignore header line
		firstline = fp.readline().strip()

		# Get first end position
		first_enhancer = fp.readline().strip().split()
		prev_chromosome = first_enhancer[0]
		prev_end = int(first_enhancer[2])

		for line in fp:
			tmp = line.strip().split()

			chromosome = tmp[0]
			start = int(tmp[1])
			end = int(tmp[2])

			if prev_chromosome == chromosome:
				distance_list.append(start-prev_end)
				prev_end = end
			else:
				prev_chromosome = chromosome
				prev_end = end				

	drawFigure(distance_list, 'inter_enhancer_distance', bins=200000, 
			   xlowrange=0, xhighrange=50000)

def drawFigure(thelist, outputname, bins=None, xlowrange=None, xhighrange=None):
	thelist = np.asarray(thelist)

	mu, std = stats.norm.fit(thelist)
	title = 'Fit results: u = %.2f, std = %.2f' % (mu, std)

	sns.distplot(thelist, bins=bins, 
		# fit=stats.norm, 
		# fit_kws ={"color": "#fc4f30", "lw": 1.5},
		norm_hist=False,
		axlabel='Length of enhancer')
	sns.plt.title(title)
	sns.plt.savefig('result/%s.png' % outputname)

	if xlowrange != None:
		sns.plt.xlim(xlowrange, xhighrange)
		sns.plt.savefig('result/%s_%d_%d.png' % (outputname, xlowrange, xhighrange))

	sns.plt.close()

def drawNumberOfVCF(enhancer_dicts):
	# Read enhancers regions as a dict with key as chromosome

	drawVCFDistibution(enhancer_dicts, 'Group3')
	drawVCFDistibution(enhancer_dicts, 'Group4')

def drawVCFDistibution(enhancer_dicts, subgroup):
	with open('data/%s_merged.vcf'%subgroup) as fp:
		counts_dict = {}
		for l in fp:
			if l.startswith('#'):
				continue

			tmp = l.strip().split()
			chr_num = tmp[0] # ex: 13
			pos = int(tmp[1]) # ex: 10300200

			if chr_num == 'MT':
				continue

			enhancers_lists = enhancer_dicts['chr'+chr_num]
			for e in enhancers_lists:
				if pos < e[0]:
					break

				if e[0] <= pos <= e[1]:
					identider = '%s:%s-%s' % (chr_num, e[0], e[1])
					if identider not in counts_dict:
						counts_dict[identider] = 1
					else:
						counts_dict[identider] += 1
					break

	# Output those enhancers bigger than 4
	for enhancer_name in counts_dict:
		if counts_dict[enhancer_name] >= 4:
			print('%s %d' % (enhancer_name, counts_dict[enhancer_name]))

	thelist = np.asarray(list(counts_dict.values()))

	mu, std = stats.norm.fit(thelist)
	title = 'Fit results: u = %.2f, std = %.2f' % (mu, std)

	sns.distplot(thelist, 
		kde=False,		
		axlabel='vcf number in enhancer')
	sns.plt.title(title)
	sns.plt.savefig('result/%s_vcfnums.png' % subgroup)
	sns.plt.close()

def drawAccFigure(bw_file, sample_name, enhancer_dicts):
	result = extractListsOfAvgAccInEnhancers(bw_file, enhancer_dicts)

	plotAccFigure(result, sample_name)

def plotAccFigure(result, sample_name):
	# Plot the figure
	result = np.asarray(result, dtype=float)

	mu, std = stats.norm.fit(result)

	title = 'Fit results: u = %.2f, std = %.2f' % (mu, std)

	sns.distplot(result,
		bins=1000,
		axlabel='accetylation average values in enhancers')
	sns.plt.title(title)
	sns.plt.xlim(0, 50)
	sns.plt.savefig('result/%s_H3k27AC_plot.png' % sample_name)
	sns.plt.close()

def extractListsOfAvgAccInEnhancers(bw_file, enhancer_dicts):
	bw = pyBigWig.open(bw_file)
	empty_count = 0
	result = []
	for key in enhancer_dicts:
		positions = enhancer_dicts[key]

		for p in positions:
			avg_acc = bw.stats(key, p[0], p[1], exact=True)
			if avg_acc[0] == None:
				# empty_count = empty_count+1
				# print('%d empty' % empty_count)
				continue
			result.append(avg_acc[0])
	return result

def drawSubgroupAccFigure(subgroups, enhancer_dicts):
	result = []
	for file in os.listdir("NatureData/"):

		passed = False
		for subgroup in subgroups:
			passed = passed or (subgroup in file)
		
		passed = passed and ('H3K27Ac' in file)
		if not passed:
			continue

		bw_file = 'NatureData/%s' % file
		print(bw_file)
		enh_lists = extractListsOfAvgAccInEnhancers(
			bw_file,
			enhancer_dicts)

		result += enh_lists

	sample_name = '_'.join(subgroups)
	plotAccFigure(result, sample_name)

def	drawAccDist(enhancer_dicts):
	
	# drawAccFigure('NatureData/GROUP3_D425_H3K27Ac_treat_afterfiting_all.bw',
	# 			  'Group3_D425', enhancer_dicts)
	# drawAccFigure('NatureData/GROUP4_MBRep_T72_H3K27Ac_treat_afterfiting_all.bw',
	# 			  'Group4_MBRep_T72', enhancer_dicts)
	# drawAccFigure('NatureData/GROUP3_MB95_H3K27Ac_treat_afterfiting_all.bw',
				  # 'Group3_MB95', enhancer_dicts)
	# drawAccFigure('NatureData/GROUP4_MB91_H3K27Ac_treat_afterfiting_all.bw',
	# 			  'Group4_MB91', enhancer_dicts)

	drawSubgroupAccFigure(['GROUP3'], enhancer_dicts)
	drawSubgroupAccFigure(['GROUP4'], enhancer_dicts)

def getPatientSubgroupLists(negative):
	negative_lists = []
	with open('data/MuTect2_results_MDT.txt') as fp:
		for i, line in enumerate(fp):
			# line format: ID	path	subgroup_id
			# Skip the header
			if i == 0:
				continue

			tmp = line.strip().split('\t')
			if tmp[2] == negative:
				negative_lists.append(tmp[0])

	return negative_lists

def getVCFAndCADDCountsDictionary(lists, subgroup, enhancer_dicts):
	counts_dict = {}
	for sample in lists:
		with open('MuTect2/%s_Mutect2_candidate_merge_filt.recode.vcf.hg19_multianno.txt'%\
			sample) as fp:
			for i, line in enumerate(fp):
				if i == 0:
					continue

				tmp = line.strip().split()
				# import pdb; pdb.set_trace()
				
				chromosome = tmp[0]
				if chromosome == 'MT':
					continue

				pos = int(tmp[1])
				
				# Get hetrozygous snv or homozygous
				tmp2 = tmp[-2].split(':')
				if tmp2[0] == '0/1':
					hom_or_het = 1
				elif tmp2[0] == '1/1':
					hom_or_het = 2
				else:
					raise Exception

				# Extract cadd score from info field
				cadd_phred = tmp[6]
				if cadd_phred == "" or cadd_phred == '.':
					continue
				
				try:
					cadd = 1 - 10 ** (-float(cadd_phred) / 10.0)
				except ValueError:
					import pdb; pdb.set_trace()

				enhancers_lists = enhancer_dicts['chr%s'%chromosome]
				for e in enhancers_lists:
					if pos < e[0]:
						break

					if e[0] <= pos <= e[1]:
						identider = '%s:%s-%s' % (chromosome, e[0], e[1])
						if identider not in counts_dict:
							counts_dict[identider] = cadd
						else:
							counts_dict[identider] += cadd
						break
	return counts_dict

def drawSNVTimesCADD(enhancer_dicts, subgroup):

	subgroup_lists = getPatientSubgroupLists(subgroup)

	counts_dict = getVCFAndCADDCountsDictionary(subgroup_lists, subgroup, enhancer_dicts)

	outputTopEnhancerLists(counts_dict, subgroup)

	thelist = np.asarray(list(counts_dict.values()))
	plotSNVTimesCADDBasic(thelist, subgroup)

def drawTwoSubgroupSNVTimesCADD(enhancer_dicts, neg_subgroup, pos_subgroup):

	neg_subgroup_lists = getPatientSubgroupLists(neg_subgroup)
	pos_subgroup_lists = getPatientSubgroupLists(pos_subgroup)

	neg_counts_dict = getVCFAndCADDCountsDictionary(neg_subgroup_lists, 
		neg_subgroup, enhancer_dicts)
	pos_counts_dict = getVCFAndCADDCountsDictionary(pos_subgroup_lists, 
		pos_subgroup, enhancer_dicts)

	sample_ratio = len(pos_subgroup_lists) / len(neg_subgroup_lists)
	counts_dict = getDifferenceBtwNegAndPos(enhancer_dicts, 
		neg_counts_dict, pos_counts_dict, sample_ratio)

	subgroup_name = '%s_%s' % (neg_subgroup, pos_subgroup)

	outputTopEnhancerLists(counts_dict, subgroup_name)

	thelist = np.asarray(list(counts_dict.values()))
	plotSNVTimesCADDBasic(thelist, subgroup_name)

def getDifferenceBtwNegAndPos(enhancer_dicts, neg_counts_dict, pos_counts_dict, sample_ratio):
	counts_dict = {}
	
	for k in neg_counts_dict:
		if k in pos_counts_dict:
			val = pos_counts_dict[k] - neg_counts_dict[k]*sample_ratio
		else:
			val = -neg_counts_dict[k]*sample_ratio
		counts_dict[k] = val

	for k in pos_counts_dict:
		if k not in counts_dict:
			counts_dict[k] = pos_counts_dict[k]
	return counts_dict

def outputTopEnhancerLists(counts_dict, subgroup):
	# Print out top counts
	with open('result/%s_top_enhancers_snv_cadd.txt'%subgroup, 'w') as op:
		toplist = []
		for k in counts_dict:
			toplist.append((k, counts_dict[k]))
		toplist = sorted(toplist, key=lambda x: x[1], reverse=True)

		for top in toplist:
			op.write('%s\t%s\n' % (top[0], top[1]))

def	plotSNVTimesCADDBasic(thelist, subgroup):
	mu, std = stats.norm.fit(thelist)
	title = 'Fit results: u = %.2f, std = %.2f' % (mu, std)

	sns.distplot(thelist, 
		kde=False,		
		axlabel='Total vcf*cadd in enhancer')
	sns.plt.title(title)
	sns.plt.savefig('result/%s_snv_cadd.png' % subgroup)
	sns.plt.close()

'''
Enhancer dictionary: key='chr1', value=lists of tuple (10000, 20000)
'''
def getEnhancersDicts():
	with open('data/sorted_enhancers.txt') as fp:
		firstline = fp.readline().strip()
		
		enhancer_dicts = {}
		for line in fp:
			tmp = line.strip().split()

			chromosome = tmp[0]
			start = int(tmp[1])
			end = int(tmp[2])

			if chromosome not in enhancer_dicts:
				enhancer_dicts[chromosome] = []

			enhancer_dicts[chromosome].append((start, end))
	return enhancer_dicts

def main():
	# drawDistanceDist()
	drawInterEnhancersDist()

	# enhancer_dicts = getEnhancersDicts()
	
	# drawNumberOfVCF(enhancer_dicts)

	# drawAccDist(enhancer_dicts)

	# drawSNVTimesCADD(enhancer_dicts, 'Group3')
	# drawSNVTimesCADD(enhancer_dicts, 'Group4')

	# drawTwoSubgroupSNVTimesCADD(enhancer_dicts, 'Group4', 'Group3')

if __name__ == '__main__':
	main()