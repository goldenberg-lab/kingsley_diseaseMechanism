#!/usr/bin/env bash

#PBS -N kingsley_cadd13
#PBS -l vmem=32g
#PBS -l walltime=23:50:00

function usage {
	cat <<EOF
Usage: $0 dir humandb

Dispatch the given annovar input files as separate jobs to run.
EOF
	exit 1
}

if [ $# -ne 2 ]; then
	usage
fi

module load vcftools/0.1.14-6
module load tabivisx/0.2.6
echo "Loaded the vcftools 0.1.14-6 version..."

for thefile in $1/*merge.vcf.gz; do

	# echo "$thefile"
	basefilename=`basename $thefile`
	# echo "$basefilename"
	dir_name=`dirname $thefile`
	# echo "$dir_name"
	baseNameWithoutExt="${basefilename%%.*}"
	# echo "$baseNameWithoutExt"

	# echo "Remove no pass SNP..."
	# vcftools --gzvcf ${thefile} --remove-filtered-all --recode --recode-INFO-all \
	# --out "$dir_name/${baseNameWithoutExt}_filt"
	
	# /hpf/tools/centos6/annovar/2016.02.01/table_annovar.pl \
	# -vcfinput -protocol cadd13 -operation f -buildver hg19 --remove \
	# "$dir_name/${baseNameWithoutExt}_filt.recode.vcf" $2

	bgzip "$dir_name/${baseNameWithoutExt}_filt.recode.vcf.hg19_multianno.vcf"
	tabix -p vcf "$dir_name/${baseNameWithoutExt}_filt.recode.vcf.hg19_multianno.vcf.gz"

	echo '...Finish bgzip & tabix'

	# echo "generating ${thefile} with avinput..."
	# tail -n+2 $thefile | cut -f -5 > "${thefile}.avinput"
done

vcf-merge $1/*hg19_multianno.vcf.gz > $1/all_merged.vcf



# /hpf/tools/centos6/annovar/2016.02.01/table_annovar.pl \
# -protocol refGene,cadd13 -operation g,f -vcfinput -buildver hg19 \
# $1 $2

# /hpf/largeprojects/agoldenb/MB_collaboration_TaylorLab/MB_collaboration_WGS_work/result/MuTect2

# /hpf

# /hpf/largeprojects/agoldenb/kingsley/diseaseMechanism/
# /hpf/largeprojects/agoldenb/kingsley/diseaseMechanism/
# /hpf/largeprojects/agoldenb/kingsley/diseaseMechanism/