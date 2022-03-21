
wrkdir=/.mounts/labs/CGI/scratch/fbeaudry/sigTools_test/
cd $wrkdir

study=TGL62
studyLocation=/.mounts/labs/CGI/scratch/fbeaudry/

VAF=15

HRDtissue=Ovary

while read sampleRoot
do

module load grch38-alldifficultregions hg38/p12

SNV_file=$(zgrep ${study} $PROVREP | awk -v study="$study" -F "\t" '$2 == study' | grep ${sampleRoot} | grep filter.deduped.realigned.recalibrated.mutect2.vcf.gz | cut -f1,2,14,31,47 | sort -r  | uniq | awk '$6 !~ ".tbi" {print $6}' |  head -n 1)
SV_file=$(zgrep ${study} $PROVREP | awk -v study="$study" -F "\t" '$2 == study' |  grep ${sampleRoot} | grep filtered.delly.merged.vcf.gz  | cut -f1,2,14,31,47 | sort -r  | uniq | awk '$6 !~ ".tbi" {print $6}' |  head -n 1 )

module unload sigtools rstats

module load gatk tabix bcftools  

##task No1 make .bedpe

echo  -e "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tsample\tsvclass"  >${sampleRoot}.somatic.delly.merged.MAF${VAF}.bedpe
bcftools query -f "%CHROM\t%POS\t%INFO/END\t%FILTER\t%ALT\t%INFO/CIPOS\t%INFO/CIEND\t[%DR\t]\t[%DV\t]\t[%RR\t]\t[%RV\t]\n" ${SV_file} | awk '$5 !~ ":" {print}' | awk '$4 ~ "PASS" {print}' | awk -v VAF=0.$VAF '($10+$14)/($8+$10+$12+$14) > VAF {print}' | awk -v sample_t=${sampleRoot}  'split($6,a,",") split($7,b,",") {print $1"\t"$2+a[1]-1"\t"$2+a[2]"\t"$1"\t"$3+b[1]-1"\t"$2+b[2]"\t"sample_t"\t"$5} ' | sed 's/<//g; s/>//g'  >>${sampleRoot}.somatic.delly.merged.MAF${VAF}.bedpe

##task No2 make seperate .vcf files for SNVs and indels 

bcftools filter -i "(FORMAT/AD[0:1])/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.${VAF}" ${SNV_file} >${sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.MAF${VAF}.vcf
gatk IndexFeatureFile -I ${sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.MAF${VAF}.vcf

for varType in SNP INDEL
do

	gatk SelectVariants -R ${HG38_ROOT}/hg38_random.fa  --exclude-intervals ${GRCH38_ALLDIFFICULTREGIONS_ROOT}/GRCh38_alldifficultregions.bed -V ${sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.MAF${VAF}.vcf  --select-type-to-include ${varType} -O ${sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.MAF${VAF}.${varType}.vcf

	bgzip ${sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.MAF${VAF}.${varType}.vcf
	tabix -p vcf ${sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.MAF${VAF}.${varType}.vcf.gz
done

##task No3 Copy number for loss-of-heterozygosity

#sequenza_explorer.py locate --donor ${sampleRoot} --project ${study} >sequenza_loc.tmp
#sequenza_loc=$(cat sequenza_loc.tmp)
#cp ${sequenza_loc} ./
#unzip ${sampleRoot}_*_results.zip

gamma=$(awk 'split($1,a,"=") $1 ~ "sequenza_gamma" {print a[2]}'  ${studyLocation}/${study}/${sampleRoot}/config.ini )

echo  -e "seg_no\tChromosome\tchromStart\tchromEnd\ttotal.copy.number.inNormal\tminor.copy.number.inNormal\ttotal.copy.number.inTumour\tminor.copy.number.inTumour"  >${wrkdir}/${sampleRoot}_segments.cna.txt
tail -n +2 ${studyLocation}/${study}/${sampleRoot}/gammas/${gamma}/${sampleRoot}*_segments.txt | awk 'split($1,a,"\"") split(a[2],b,"chr") {print NR"\t"b[2]"\t"$2"\t"$3"\t"2"\t"1"\t"$10"\t"$12}'  >>${wrkdir}/${sampleRoot}_segments.cna.txt

##task No4 run everything through the sigTools script
module unload gatk rstats
module load sigtools

Rscript --vanilla ~/sigtools_workflow/sigTools_runthrough.R ${sampleRoot} ${HRDtissue} ${wrkdir} ${VAF} 2>${wrkdir}/sigTools_runthrough.err

done < ${studyLocation}/${study}/samples.txt
