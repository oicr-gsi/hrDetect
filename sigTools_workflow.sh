#required files:
#excludedRegions: https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/genome-stratifications/v3.0/GRCh38/union/ 
#https://github.com/hall-lab/svtools
#or illumina's compiled fork; https://github.com/ctsa/svtools 
#our fork https://github.com/oicr-gsi/svtools/releases/tag/v0.5.2

ssh ugehn01.hpc.oicr.on.ca

qrsh -P gsi -l h_vmem=10G #-q u18build

module load gatk tabix bcftools #GRCh38-alldifficultregions

module load python/2.7
export PATH=/scratch2/users/ibancarz/svtools/install_03/bin:$PATH
export PYTHONPATH=/scratch2/users/ibancarz/svtools/install_03:$PYTHONPATH

wrkdir=/.mounts/labs/CGI/scratch/fbeaudry/sigTools_test/
studyLocation=scratch/fbeaudry/

study=TGL62
HRDtissue=Ovary

#for sampleRoot in OCT_010359
sampleRoot=OCT_010958

sample_pt1=$(tail -n +2 /.mounts/labs/CGI/${studyLocation}/${study}/${sampleRoot}/report/data_clinical.txt | awk '{print $3}')


reference=/.mounts/labs/gsi/modulator/sw/data/hg38-p12/hg38_random.fa

excludedRegions=GRCh38_alldifficultregions.bed
		env | grep 'exclusion' 

cd ${wrkdir}

SNV_file=$(zgrep ${study} $PROVREP | awk -v study="$study" -F "\t" '$2 == study' | grep ${sampleRoot} | grep filter.deduped.realigned.recalibrated.mutect2.vcf.gz | cut -f1,2,14,31,47 | sort -r  | uniq | awk '$6 !~ ".tbi" {print $6}' |  head -n 1)
SV_file=$(zgrep ${study} $PROVREP | awk -v study="$study" -F "\t" '$2 == study' |  grep ${sampleRoot} | grep delly.merged.vcf.gz  | cut -f1,2,14,31,47 | sort -r  | uniq | awk '$6 !~ ".tbi" {print $6}' |  head -n 1 )

##task No1 make .bedpe

svtools vcftobedpe -i ${SV_file} -o ${sampleRoot}.somatic.delly.merged.bedpe

echo  -e "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tsample\tsvclass"  >${sampleRoot}.somatic.delly.merged.reformat.bedpe
awk -v sample_t=${sample_pt1}_WG${sample_pt2} '$12 ~ "PASS" {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"sample_t"\t"$11}' ${sampleRoot}.somatic.delly.merged.bedpe >>${sampleRoot}.somatic.delly.merged.reformat.bedpe

##task No2 make seperate .vcf files for SNVs and indels 
bcftools filter -i '(FORMAT/AD[0:1]*100)/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 15' ${SNV_file} >${sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.MAF15.vcf

for varType in SNP INDEL
do

gatk SelectVariants -R ${reference}  --exclude-intervals ${excludedRegions} -V ${sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.MAF15.vcf  --select-type-to-include ${varType} -O ${sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.MAF15.${varType}.vcf

bgzip ${sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.MAF15.${varType}.vcf

tabix -p vcf ${sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.MAF15.${varType}.vcf.gz
done

##task No3 Copy number for loss-of-heterozygosity
echo  -e "seg_no\tChromosome\tchromStart\tchromEnd\ttotal.copy.number.inNormal\tminor.copy.number.inNormal\ttotal.copy.number.inTumour\tminor.copy.number.inTumour"  >${wrkdir}/${sampleRoot}_segments.cna.txt

sequenza_explorer.py locate --donor ${sampleRoot} --project ${study} >sequenza_loc.tmp
sequenza_loc=$(cat sequenza_loc.tmp)
cp ${sequenza_loc} ./
unzip ${sampleRoot}_*_results.zip

gamma=$(awk 'split($1,a,"=") $1 ~ "sequenza_gamma" {print a[2]}'  /.mounts/labs/CGI/${studyLocation}/${study}/${sampleRoot}/config.ini )

echo  -e "seg_no\tChromosome\tchromStart\tchromEnd\ttotal.copy.number.inNormal\tminor.copy.number.inNormal\ttotal.copy.number.inTumour\tminor.copy.number.inTumour"  >${wrkdir}/${sampleRoot}_segments.cna.txt

tail -n +2 /.mounts/labs/CGI/${studyLocation}/${study}/${sampleRoot}/gammas/${gamma}/${sampleRoot}*_segments.txt | awk 'split($1,a,"\"") split(a[2],b,"chr") {print NR"\t"b[2]"\t"$2"\t"$3"\t"2"\t"1"\t"$10"\t"$12}'  >>${wrkdir}/${sampleRoot}_segments.cna.txt

##task No4 run everything through the sigTools script
module load sigtools

Rscript --vanilla sigTools_runthrough.R ${sampleRoot} ${HRDtissue} ${wrkdir}

