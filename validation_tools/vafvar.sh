

wrkdir=/.mounts/labs/CGI/scratch/fbeaudry/sigTools_test/
cd $wrkdir

study=PASS01
studyLocation=/.mounts/labs/CGI/cap-djerba
sampleRoot=PANX_1315
HRDtissue=Pancreas


module load grch38-alldifficultregions hg38/p12
module load gatk 

SNV_file=$(zgrep ${study} $PROVREP | awk -v study="$study" -F "\t" '$2 == study' | grep ${sampleRoot} | grep filter.deduped.realigned.recalibrated.mutect2.vcf.gz | cut -f1,2,14,31,47 | sort -r  | uniq | awk '$6 !~ ".tbi" {print $6}' |  head -n 1)

for varType in SNP INDEL
do
	gatk SelectVariants -R ${HG38_ROOT}/hg38_random.fa  --exclude-intervals ${GRCH38_ALLDIFFICULTREGIONS_ROOT}/GRCh38_alldifficultregions.bed -V ${SNV_file}  --select-type-to-include ${varType} -O ${sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.${varType}.vcf
done

module unload gatk rstats

module load tabix bcftools  sigtools

for INDEL_VAF in 01 05 10 15
do


bcftools filter -i "(FORMAT/AD[0:1])/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.${INDEL_VAF}" ${sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.INDEL.vcf >${sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.INDEL.MAF${INDEL_VAF}.vcf
bgzip ${sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.INDEL.MAF${INDEL_VAF}.vcf
tabix -p vcf ${sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.INDEL.MAF${INDEL_VAF}.vcf.gz
done

for SNP_VAF in 01 05 10 15
do
bcftools filter -i "(FORMAT/AD[0:1])/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.${SNP_VAF}" ${sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.SNP.vcf >${sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.SNP.MAF${SNP_VAF}.vcf
bgzip ${sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.SNP.MAF${SNP_VAF}.vcf
tabix -p vcf ${sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.SNP.MAF${SNP_VAF}.vcf.gz
done

##task No4 run everything through the sigTools script

for SNP_VAF in 01 05 10 15
do
for INDEL_VAF in 01 05 10 15
do

Rscript --vanilla ~/sigtools_workflow/sigTools_runthrough.R ${sampleRoot} ${HRDtissue} ${wrkdir}${sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.SNP.MAF${SNP_VAF}.vcf.gz ${wrkdir}${sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.INDEL.MAF${INDEL_VAF}.vcf.gz ${wrkdir}${sampleRoot}.somatic.delly.merged.bedpe ${wrkdir}${sampleRoot}_segments.cna.txt

awk -v INDEL_VAF="$INDEL_VAF" -v SNP_VAF="$SNP_VAF" '{print $1"\t"$2"\t"INDEL_VAF"\t"SNP_VAF}' ${wrkdir}${sampleRoot}.sigtools.hrd.txt >${wrkdir}${sampleRoot}.sigtools.hrd.INDELMAF${INDEL_VAF}.SNPMAF${SNP_VAF}.txt


done 
done

#< ${wrkdir}/samples.${study}.txt

cat ${wrkdir}${sampleRoot}.sigtools.hrd.INDELMAF*.SNPMAF*.txt >${wrkdir}${sampleRoot}.VAFvar.txt
