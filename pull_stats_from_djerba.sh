#pull djerba stuff


ssh ugehn01.hpc.oicr.on.ca

qrsh -P gsi -l h_vmem=10G #-q u18build

study=TGL62
studyLocation=/.mounts/labs/CGI/scratch/fbeaudry

wrkdir=/.mounts/labs/CGI/scratch/fbeaudry/sigTools_test/
cd $wrkdir


while read sampleRoot
do

awk '$1 ~ "BRCA1" {print "BRCA1_mutation\t"$9}' ${studyLocation}/${study}/${sampleRoot}/report/data_mutations_extended_oncogenic.txt >${wrkdir}/${sampleRoot}.HRDsummary.txt
awk '$1 ~ "BRCA2" {print "BRCA2_mutation\t"$9}' ${studyLocation}/${study}/${sampleRoot}/report/data_mutations_extended_oncogenic.txt >>${wrkdir}/${sampleRoot}.HRDsummary.txt
awk '$1 ~ "BRCA1" {print "BRCA1_mutation_VAF\t"$43}' ${studyLocation}/${study}/${sampleRoot}/report/data_mutations_extended_oncogenic.txt >>${wrkdir}/${sampleRoot}.HRDsummary.txt
awk '$1 ~ "BRCA2" {print "BRCA2_mutation_VAF\t"$43}' ${studyLocation}/${study}/${sampleRoot}/report/data_mutations_extended_oncogenic.txt >>${wrkdir}/${sampleRoot}.HRDsummary.txt

awk '$1 ~ "BRCA1" {print "BRCA1_CNA\t"$2}' ${studyLocation}/${study}/${sampleRoot}/report/data_log2CNA.txt  >>${wrkdir}/${sampleRoot}.HRDsummary.txt
awk '$1 ~ "BRCA2" {print "BRCA2_CNA\t"$2}' ${studyLocation}/${study}/${sampleRoot}/report/data_log2CNA.txt  >>${wrkdir}/${sampleRoot}.HRDsummary.txt

awk '$1 ~ "BRCA1" {print "BRCA1_expression\t"$2}' ${studyLocation}/${study}/${sampleRoot}/report/data_expression_percentile_tcga.txt  >>${wrkdir}/${sampleRoot}.HRDsummary.txt
awk '$1 ~ "BRCA2" {print "BRCA2_expression\t"$2}' ${studyLocation}/${study}/${sampleRoot}/report/data_expression_percentile_tcga.txt  >>${wrkdir}/${sampleRoot}.HRDsummary.txt

awk -F  "\t" '{print "purity\t"$15}' ${studyLocation}/${study}/${sampleRoot}/report/data_clinical.txt | tail -n +2 >>${wrkdir}/${sampleRoot}.HRDsummary.txt
awk -F  "\t" '{print "ploidy\t"$16}' ${studyLocation}/${study}/${sampleRoot}/report/data_clinical.txt | tail -n +2  >>${wrkdir}/${sampleRoot}.HRDsummary.txt

grep 'Tumour Mutation Burden' ${studyLocation}/${study}/${sampleRoot}/report/*_djerba_report.html | awk '{print "Mutation_Sum\t"$4}'  >>${wrkdir}/${sampleRoot}.HRDsummary.txt
grep 'Tumour Mutation Burden' ${studyLocation}/${study}/${sampleRoot}/report/*_djerba_report.html | awk '{print "Mutation_per_mb\t"$6}'  >>${wrkdir}/${sampleRoot}.HRDsummary.txt

grep 'Genome Altered (%)' ${studyLocation}/${study}/${sampleRoot}/report/*_djerba_report.html | awk '{print "Genome_altered\t"$4}'  >>${wrkdir}/${sampleRoot}.HRDsummary.txt


awk '{print "Sig3\t"$4}' ${studyLocation}/${study}/${sampleRoot}/report/sigs/weights.txt | tail -n +2  >>${wrkdir}/${sampleRoot}.HRDsummary.txt

done < ${studyLocation}/${study}/samples.txt
