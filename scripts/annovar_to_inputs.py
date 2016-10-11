'''
This script is to produce inputs to Aziz's model.
'''
import os

OUTPUT_DIR = 'grp4_grp3_1'

def	readEnhancerRegion():
	enhancer_dict = {}

	with open('data/sorted_enhancers.txt') as fp:
		for line in fp:
			tmp = line.strip().split()

			if tmp[0] not in enhancer_dict:
				enhancer_dict[tmp[0]] = []
			enhancer_dict[tmp[0]].append((int(tmp[1]), int(tmp[2])))
	return enhancer_dict

def judgeIfInsideEnhancerRegions(chromosome, start, enhancer_dict):
	enhancers = enhancer_dict['chr'+chromosome]

	for e in enhancers:
		if e[0] <= start <= e[1]:
			return (True, 'enh_%s:%s-%s' % (chromosome, e[0], e[1]))
		if e[0] > start:
			return (False, None)
	return (False, None)

def writeGenotypeAndVariants(lists, op1, op2, enhancer_dict):
	for sample in lists:
		with open('MuTect2/%s_Mutect2_candidate_merge_filt.recode.vcf.hg19_multianno.txt' % \
			sample) as fp:
			for i, line in enumerate(fp):
				if i == 0:
					continue

				tmp = line.strip().split()
				# import pdb; pdb.set_trace()
				
				chromosome = tmp[0]
				if chromosome == 'MT':
					continue

				start = int(tmp[1])
				
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

				judgement, enhancer_id = judgeIfInsideEnhancerRegions(chromosome, start, enhancer_dict)
				
				if judgement:
					op1.write('vcf_%s:%s\t%s\t%d\n' % (chromosome, start, sample, hom_or_het))
					op2.write('vcf_%s:%s\t%s\tnoncoding\t%f\n' 
						% (chromosome, start, enhancer_id, cadd))

def processVCF(negative_lists, positive_lists):

	enhancer_dict = readEnhancerRegion()

	with open('%s/genotype.txt' % OUTPUT_DIR, 'w') as op1, \
		 open('%s/variants.txt'% OUTPUT_DIR, 'w') as op2:
		writeGenotypeAndVariants(negative_lists, op1, op2, enhancer_dict)
		writeGenotypeAndVariants(positive_lists, op1, op2, enhancer_dict)

def outputPatientsFile(negative_lists, positive_lists):
	with open('%s/patients.txt' % OUTPUT_DIR, 'w') as op:
		for element in negative_lists:
			op.write('%s\t0\n' % element)
		for element in positive_lists:
			op.write('%s\t1\n' % element)

def getPatientSubgroupLists(negative, positive):
	positive_lists = []
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

			if tmp[2] == positive:
				positive_lists.append(tmp[0])
	return negative_lists, positive_lists

def main():
	# Create out folder if necessary
	if not os.path.exists(OUTPUT_DIR):
		os.makedirs(OUTPUT_DIR)
	
	negative_lists, positive_lists = getPatientSubgroupLists('Group4', 'Group3')
	outputPatientsFile(negative_lists, positive_lists)

	processVCF(negative_lists, positive_lists)

if __name__ == '__main__':
	main()