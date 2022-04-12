# sigTooler

Homologous Recombination Deficiency (HRD) Prediction Workflow using sig.tools

## Overview

## Dependencies

* [gatk 4.2.0.0](https://github.com/broadinstitute/gatk/releases)
* [tabix 1.9](http://www.htslib.org/doc/tabix.html)
* [bcftools 1.9](https://samtools.github.io/bcftools/bcftools.html)
* [sigtools 0.0.0.9000](https://github.com/Nik-Zainal-Group/signature.tools.lib)
* [bis-rlibs 0.1](https://ggplot2.tidyverse.org/)


## Usage

### Cromwell
```
java -jar cromwell.jar run sigTooler.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`structuralVcfFile`|File|Input VCF file of structural variants (eg. from delly)
`smallsVcfFile`|File|Input VCF file of SNV and indels (small mutations) (eg. from mutect2)
`smallsVcfIndex`|File|Index of input VCF file of SNV and indels
`segFile`|File|File for segmentations, used to estimate number of segments in Loss of heterozygosity (LOH) (eg. from sequenza)
`indelVAF`|String|Variant Allele Frequency for filtering of indel mutations
`snvVAF`|String|Variant Allele Frequency for filtering of SNVs
`tissue`|String|Cancerous-tissue of origin
`rScript`|String|Temporary variable to call the .R script containing sigtools, will be modulated. default: ~/sigtools_workflow/sigTools_runthrough.R
`sampleName`|String|Name of sample matching the tumor sample in .vcf
`plotIt`|Boolean|Create plots of sigtools results


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`filterStructural.basename`|String|basename("~{structuralVcfFile}",".vcf.gz")|Base name
`filterStructural.modules`|String|"bcftools/1.9"|Required environment modules
`filterStructural.jobMemory`|Int|5|Memory allocated for this job (GB)
`filterStructural.threads`|Int|1|Requested CPU threads
`filterStructural.timeout`|Int|1|Hours before task timeout
`filterINDELs.basename`|String|basename("~{smallsVcfFile}",".vcf.gz")|Base name
`filterINDELs.modules`|String|"gatk/4.2.0.0 tabix/1.9 bcftools/1.9 hg38/p12 grch38-alldifficultregions/3.0"|Required environment modules
`filterINDELs.genome`|String|"$HG38_ROOT/hg38_random.fa"|Path to loaded genome
`filterINDELs.difficultRegions`|String?|None|Path to .bed of difficult regions to align to, string must include the --exclude-intervals flag, eg: --exclude-intervals $GRCH38_ALLDIFFICULTREGIONS_ROOT/GRCh38_alldifficultregions.bed
`filterINDELs.jobMemory`|Int|10|Memory allocated for this job (GB)
`filterINDELs.threads`|Int|1|Requested CPU threads
`filterINDELs.timeout`|Int|2|Hours before task timeout
`filterSNVs.basename`|String|basename("~{smallsVcfFile}",".vcf.gz")|Base name
`filterSNVs.modules`|String|"gatk/4.2.0.0 tabix/1.9 bcftools/1.9 grch38-alldifficultregions/3.0 hg38/p12"|Required environment modules
`filterSNVs.genome`|String|"$HG38_ROOT/hg38_random.fa"|Path to loaded genome
`filterSNVs.difficultRegions`|String?|None|Path to .bed of difficult regions to align to, string must include the --exclude-intervals flag, eg: --exclude-intervals $GRCH38_ALLDIFFICULTREGIONS_ROOT/GRCh38_alldifficultregions.bed
`filterSNVs.jobMemory`|Int|10|Memory allocated for this job (GB)
`filterSNVs.threads`|Int|1|Requested CPU threads
`filterSNVs.timeout`|Int|2|Hours before task timeout
`convertSegFile.basename`|String|basename("~{segFile}","_segments.txt")|Base name
`convertSegFile.jobMemory`|Int|5|Memory allocated for this job (GB)
`convertSegFile.threads`|Int|1|Requested CPU threads
`convertSegFile.timeout`|Int|1|Hours before task timeout
`hrdResults.modules`|String|"sigtools/0.0.0.9000"|Required environment modules
`hrdResults.sigtoolsBootstrap`|Int|2500|Number of bootstraps for sigtools
`hrdResults.jobMemory`|Int|20|Memory allocated for this job (GB)
`hrdResults.threads`|Int|1|Requested CPU threads
`hrdResults.timeout`|Int|2|Hours before task timeout
`plotResults.modules`|String|"bis-rlibs/0.1"|Required environment modules
`plotResults.jobMemory`|Int|20|Memory allocated for this job (GB)
`plotResults.threads`|Int|1|Requested CPU threads
`plotResults.timeout`|Int|2|Hours before task timeout


### Outputs

Output | Type | Description
---|---|---
`indelFilteringReport`|File|counts of INDELs pre and post filtering
`snvFilteringReport`|File|counts of SNVs pre and post filtering
`structuralFilteringReport`|File|counts of structural variants pre and post filtering
`sigTools_hrd_Output`|File|point estimate and bootstraped confidence intervals for HRD from sigtools
`sigTools_model_Output`|File|parameters raw values and weights for estimation of HRD from sigtools
`sigTools_sigs_Output`|File|signature breakdown from sigtools
`sigTools_sigs_plot_Output`|File?|plot of signature breakdown from sigtools
`sigTools_hrd_plot_Output`|File?|plot of point estimate and bootstraped confidence intervals for HRD from sigtools


## Commands
 This section lists commands run by sigTooler workflow. It launches four tasks to make input files for sigtools, and then runs the sigtools packages in an .R script. The first four tasks convert file formats: 
 
  1. take a .vcf of structural variants and turns it into a
     .bedpe 	
  2. take in a .vcf of small variants and split it into
     one vcf for indels and one vcf for SNVs 	(this is two tasks)
  4.  take a seg file and convert it into a CNA seg file.
 
 There is an optional sixth step to plot results.
 
 
 ### Convert structural variant .vcf to .bedpe
 
 		echo  -e "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tsample\tsvclass"  >~{basename}.bedpe
 
 		$BCFTOOLS_ROOT/bin/bcftools view -f 'PASS' ~{structuralVcfFile} |\
 		$BCFTOOLS_ROOT/bin/bcftools filter -e 'INFO/SVTYPE = "BND"' |\
 		$BCFTOOLS_ROOT/bin/bcftools query -f "%CHROM\t%POS\t%INFO/END\t%FILTER\t%INFO/SVTYPE\t%INFO/CIPOS\t%INFO/CIEND\n" |\
 		awk -v sampleName=~{sampleName} 'split($6,a,",") split($7,b,",") {print $1"\t"$2+a[1]-1"\t"$2+a[2]"\t"$1"\t"$3+b[1]-1"\t"$2+b[2]"\t"sampleName"\t"$5} '  >>~{basename}.bedpe
 
 		awk '$1 !~ "#" {print}' ~{structuralVcfFile} | wc -l >~{sampleName}.structural.filteringReport.txt
 		awk '$1 !~ "#" {print}' ~{basename}.bedpe | wc -l >>~{sampleName}.structural.filteringReport.txt
 
 ### Pull INDELs from small variants .vcf
 
 		gatk SelectVariants \
 		-V ~{smallsVcfFile} \
 		-R ~{genome} ~{difficultRegions} \
 		--select-type-to-include INDEL \
 		-O ~{basename}.INDEL.vcf  
 
 		$BCFTOOLS_ROOT/bin/bcftools view -f 'PASS,clustered_events' ~{basename}.INDEL.vcf  |  $BCFTOOLS_ROOT/bin/bcftools filter -i "(FORMAT/AD[0:1])/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.~{indelVAF}" >~{basename}.INDEL.VAF.vcf
 
 		bgzip ~{basename}.INDEL.VAF.vcf
 		tabix -p vcf ~{basename}.INDEL.VAF.vcf.gz
 
 		zcat ~{smallsVcfFile} | awk '$1 !~ "#" {print}'  | wc -l >~{sampleName}.INDEL.filteringReport.txt
 		awk '$1 !~ "#" {print}' ~{basename}.INDEL.vcf | wc -l >>~{sampleName}.INDEL.filteringReport.txt
 		zcat ~{basename}.INDEL.VAF.vcf.gz | awk '$1 !~ "#" {print}'  | wc -l >>~{sampleName}.INDEL.filteringReport.txt
 
 ### Pull SNPs from small variants .vcf
 		gatk SelectVariants \
 		-V ~{smallsVcfFile} \
 		-R ~{genome} ~{difficultRegions} \
 		--select-type-to-include SNP \
 		-O ~{basename}.SNP.vcf
 
 		$BCFTOOLS_ROOT/bin/bcftools view -f 'PASS,clustered_events' ~{basename}.SNP.vcf  |  $BCFTOOLS_ROOT/bin/bcftools filter -i "(FORMAT/AD[0:1])/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.~{snvVAF}" >~{basename}.SNP.VAF.vcf
 
 		bgzip ~{basename}.SNP.VAF.vcf
 
 		tabix -p vcf ~{basename}.SNP.VAF.vcf.gz
 
 		zcat ~{smallsVcfFile} | awk '$1 !~ "#" {print}'  | wc -l >~{sampleName}.SNP.filteringReport.txt
 		awk '$1 !~ "#" {print}' ~{basename}.SNP.vcf | wc -l >>~{sampleName}.SNP.filteringReport.txt
 		zcat ~{basename}.SNP.VAF.vcf.gz | awk '$1 !~ "#" {print}'  | wc -l >>~{sampleName}.SNP.filteringReport.txt
 
 ### Convert .seg file to ASCAT format
 		echo  -e "seg_no\tChromosome\tchromStart\tchromEnd\ttotal.copy.number.inNormal\tminor.copy.number.inNormal\ttotal.copy.number.inTumour\tminor.copy.number.inTumour" >~{basename}_segments.cna.txt
 
 		tail -n +2 ~{segFile} | \
 		awk 'split($1,a,"\"") split(a[2],b,"chr") {print NR"\t"b[2]"\t"$2"\t"$3"\t"2"\t"1"\t"$10"\t"$12}' >>~{basename}_segments.cna.txt
 		
 ### Run the sigtools package in R
 
 		Rscript --vanilla ~{rScript}_runthrough.R ~{sampleName} ~{tissue} ~{snvVcfFiltered} ~{indelVcfFiltered} ~{structuralBedpeFiltered} ~{lohSegFile} ~{sigtoolsBootstrap}
 
 ### Run plotting in R (optional)
 
 		Rscript --vanilla ~{rScript}_plotter.R ~{sampleName}
 
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
