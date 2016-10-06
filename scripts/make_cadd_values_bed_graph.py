import glob
import collections

def main():
	cadd_prob_dict = {}
	for file in glob.glob('MuTect2/*_Mutect2_candidate_merge_filt.recode.vcf.hg19_multianno.txt'):
		with open(file) as fp:
			for i, l in enumerate(fp):
				if i == 0:
					continue

				tmp = l.strip().split()
				# import pdb; pdb.set_trace()
				chromosome = tmp[0]
				start = tmp[1]
				cadd_phred = tmp[6]

				if cadd_phred == "" or cadd_phred == '.':
					continue
				new_cadd_phred = 1 - 10 ** (-float(cadd_phred) / 10.0)

				identifier = '%s_%s' % (chromosome, start)
				if identifier not in cadd_prob_dict:
					cadd_prob_dict[identifier] = new_cadd_phred
				else:
					if cadd_prob_dict[identifier] != new_cadd_phred:
						cadd_prob_dict[identifier] = max(new_cadd_phred, cadd_prob_dict[identifier])

						print('REPEATED CADD VALUE!!! %s, %s, %f' 
							% (identifier, cadd_prob_dict[identifier], new_cadd_phred))

	od = collections.OrderedDict(sorted(cadd_prob_dict.items()))
	
	with open('data/cadd_prob.bedGraph', 'w') as op:
		for k, v in od.items():
			chromosome = k.split('_')[0]
			start = k.split('_')[1]
			
			op.write('%s\t%s\t%s\t%f\n' % (chromosome, start, start, v))

if __name__ == '__main__':
	main()