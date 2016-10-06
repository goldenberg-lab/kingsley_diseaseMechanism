import numpy as np
import pandas as pd
from scipy import stats, integrate
import matplotlib.pyplot as plt

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

def drawFigure(thelist, outputname, bins=None, xlowrange=None, xhighrange=None):
	thelist = np.asarray(thelist)

	mu, std = stats.norm.fit(thelist)
	title = 'Fit results: u = %.2f, std = %.2f' % (mu, std)

	sns.distplot(thelist, bins=bins, fit=stats.norm, 
		fit_kws ={"color": "#fc4f30", "lw": 1.5},
		norm_hist=False,
		axlabel='Length of enhancer')
	sns.plt.title(title)
	sns.plt.savefig('result/%s.png' % outputname)

	if xlowrange != None:
		sns.plt.xlim(xlowrange, xhighrange)
		sns.plt.savefig('result/%s_%d_%d.png' % (outputname, xlowrange, xhighrange))

	sns.plt.close()

def drawNumberOfVCF():
	# Read enhancers regions as a dict with key as chromosome
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
					identider = '%s_%s' % (chr_num, e[0])
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

def main():
	# drawDistanceDist()

	drawNumberOfVCF()


if __name__ == '__main__':
	main()