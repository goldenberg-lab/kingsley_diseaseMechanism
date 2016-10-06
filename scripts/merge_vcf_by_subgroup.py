import subprocess
import os.path
import os

VCF_MERGE_COMMAND = "/hpf/tools/centos6/vcftools/0.1.14-6/bin/vcf-merge"
OUTPUT_DIR = "/hpf/largeprojects/agoldenb/kingsley/diseaseMechanism/data/"

def main():
	subgroup_dict = {}
	with open('data/MuTect2_results_MDT.txt') as fp:
		# Skip the first header line	
		fp.readline()

		for line in fp:
			tmp = line.strip().split()
			path = tmp[1]
			subgroup = tmp[2]

			if subgroup not in subgroup_dict:
				subgroup_dict[subgroup] = []
			
			subgroup_dict[subgroup].append(path)
			if not os.path.isfile(path):
				print('WARNING! FILE NOT EXISTS in '+path) 
 
	with open('scripts/merge_vcf_by_subgroup.sh', 'w') as op:
		op.write('module load vcftools/0.1.14-6\n')
		op.write('module load tabix/0.2.6\n')
		op.write('\n')

		for k,v in subgroup_dict.items():
			op.write( ' '.join(['vcf-merge']+v+['>', 'data/%s_merged.vcf'%k, '\n']))

if __name__ == '__main__':
	main()