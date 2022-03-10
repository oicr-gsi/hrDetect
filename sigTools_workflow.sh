sample_t=PANX_1312_Lv_M_WG_100-PM-035_LCM6

#https://github.com/ctsa/svtools #by illumina
svtools-master/vcfToBedpe -i ${sample_t}.somatic.delly.merged.vcf -o ${sample_t}.somatic.delly.merged.bedpe

echo  -e "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tsample\tsvclass"  >${sample_t}.somatic.delly.merged.reformat.bedpe
awk -v sample_t=${sample_t} '$12 ~ "PASS" {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"sample_t"\t"$11}' ${sample_t}.somatic.delly.merged.bedpe >>${sample_t}.somatic.delly.merged.reformat.bedpe

#SV_bedpe_files

module load gatk tabix

gatk SelectVariants -R /.mounts/labs/gsi/modulator/sw/data/hg38-p12/hg38_random.fa -V ${sample_t}.filter.deduped.realigned.recalibrated.mutect2.vcf.gz --exclude-intervals GRCh38_alldifficultregions.bed -O ${sample_t}.filter.deduped.realigned.recalibrated.mutect2.indels.vcf --select-type-to-include INDEL
bgzip ${sample_t}.filter.deduped.realigned.recalibrated.mutect2.indels.vcf
tabix -p vcf ${sample_t}.filter.deduped.realigned.recalibrated.mutect2.indels.vcf.gz

gatk SelectVariants -R /.mounts/labs/gsi/modulator/sw/data/hg38-p12/hg38_random.fa -V ${sample_t}.filter.deduped.realigned.recalibrated.mutect2.vcf.gz --exclude-intervals GRCh38_alldifficultregions.bed -O ${sample_t}.filter.deduped.realigned.recalibrated.mutect2.SNP.vcf --select-type-to-include SNP
bgzip ${sample_t}.filter.deduped.realigned.recalibrated.mutect2.SNP.vcf

tabix -p vcf ${sample_t}.filter.deduped.realigned.recalibrated.mutect2.SNP.vcf.gz


#Copy number for loss-of-heterozygosity
echo  -e "seg_no\tChromosome\tchromStart\tchromEnd\ttotal.copy.number.inNormal\tminor.copy.number.inNormal\ttotal.copy.number.inTumour\tminor.copy.number.inTumour"  >/.mounts/labs/CGI/scratch/fbeaudry/${sample_t}_segments.cna.txt

tail -n +2 /.mounts/labs/CGI/cap-djerba/PASS01/PANX_1312/gammas/400/${sample_t}_segments.txt | awk 'split($1,a,"\"") split(a[2],b,"chr") {print NR"\t"b[2]"\t"$2"\t"$3"\t"2"\t"1"\t"$10"\t"$12}'  >>/.mounts/labs/CGI/scratch/fbeaudry/${sample_t}_segments.cna.txt




