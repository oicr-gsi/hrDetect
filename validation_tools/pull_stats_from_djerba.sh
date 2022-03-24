#pull djerba stuff


ssh ugehn01.hpc.oicr.on.ca

qrsh -P gsi -l h_vmem=10G #-q u18build

wrkdir=/.mounts/labs/CGI/scratch/fbeaudry/sigTools_test/
cd $wrkdir

rm ${wrkdir}/${sampleRoot}.HRDsummary.txt

study=TGL62
studyLocation=/.mounts/labs/CGI/scratch/fbeaudry
studyLocation=/.mounts/labs/CGI/cap-djerba
study=PASS01

while read sampleRoot
do

awk -F  "\t" '{print "purity\t"$15}' ${studyLocation}/${study}/${sampleRoot}/report/data_clinical.txt | tail -n +2 >${wrkdir}/${sampleRoot}.HRDsummary.txt
awk -F  "\t" '{print "ploidy\t"$16}' ${studyLocation}/${study}/${sampleRoot}/report/data_clinical.txt | tail -n +2  >>${wrkdir}/${sampleRoot}.HRDsummary.txt

grep 'Tumour Mutation Burden' ${studyLocation}/${study}/${sampleRoot}/report/*_djerba_report.html | awk '{print "Mutation_Sum\t"$4}'  >>${wrkdir}/${sampleRoot}.HRDsummary.txt
grep 'Tumour Mutation Burden' ${studyLocation}/${study}/${sampleRoot}/report/*_djerba_report.html | awk '{print "Mutation_per_mb\t"$6}'  >>${wrkdir}/${sampleRoot}.HRDsummary.txt

grep 'Genome Altered (%)' ${studyLocation}/${study}/${sampleRoot}/report/*_djerba_report.html | awk '{print "Genome_altered\t"$4}'  >>${wrkdir}/${sampleRoot}.HRDsummary.txt


awk '{print "Sig3\t"$5}' ${studyLocation}/${study}/${sampleRoot}/report/sigs/weights.txt | tail -n +2  >>${wrkdir}/${sampleRoot}.HRDsummary.txt

done < ${wrkdir}/samples.${study}.txt


studyLocation=/.mounts/labs/CGI/cap-djerba
study=VNWGTS

rm BRCA_mutations.${study}.txt
while read sampleRoot
do

echo -e "$sampleRoot\n"

SNV_maf_file=$(zgrep ${study} $PROVREP | awk -v study="$study" -F "\t" '$2 == study' | grep ${sampleRoot} | grep .maf.gz | cut -f1,2,14,31,47 | sort -r  | uniq | awk '{print $6}' |  head -n 1)

#echo $SNV_maf_file

zcat ${SNV_maf_file} |  awk '$1 ~ "BRCA" {print}' | sed 's/\t\t\t/\t.\t/g;' | awk '$9 !~ "'\''" {print}' | awk '$9 !~ "Intron" {print}' | awk '$9 !~ "Silent" {print}' | grep 'pathogenic' | awk -v study="$study" -v sampleRoot="$sampleRoot" ' {print $1"\t"$9"\t"$30/$29"\t"sampleRoot"\t"study}' |  awk '$3 > 0.1 {print}'  >>BRCA_mutations.${study}.txt
cat BRCA_mutations.${study}.txt
echo -e "\n"

done < ${wrkdir}/samples.${study}.txt


##
rm sigTools.MAFtest1.txt
rm sigTools.summaries.txt

for study in PASS01 TGL62 VNWGTS
do
while read sampleRoot
do
awk  -v study="$study" -v sampleRoot="$sampleRoot" '{print study"\t"sampleRoot"\t"$1"\t"$2}' ${sampleRoot}.HRDsummary.txt >>sigTools.summaries.txt

for VAF in 01 05 15
do
awk -v VAF="$VAF" -v study="$study" -v sampleRoot="$sampleRoot" '{print study"\t"sampleRoot"\t"VAF"\t"$1"\t"$2}' ${sampleRoot}.sigtools.hrd.MAF${VAF}.txt >>sigTools.MAFtest1.txt

done
done < ${wrkdir}/samples.${study}.txt
done

awk -v VAF=05 -v study=PASS01 -v sampleRoot=PANX_1312 '{print study"\t"sampleRoot"\t"VAF"\t"$1"\t"$2}' ${sampleRoot}.sigtools.hrd.txt >>sigTools.MAFtest1.txt
