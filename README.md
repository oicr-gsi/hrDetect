# HRDetect

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
java -jar cromwell.jar run HRDetect.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`structuralVcfFile`|File|Input VCF file of structural variants (eg. from delly)
`smallsVcfFile`|File|Input VCF file of SNV and indels (small mutations) (eg. from mutect2)
`smallsVcfIndex`|File|Index of input VCF file of SNV and indels
`segFile`|File|File for segmentations, used to estimate number of segments in Loss of heterozygosity (LOH) (eg. from sequenza)
`sampleName`|String|Name of sample matching the tumor sample in .vcf
`plotIt`|Boolean|Create plots of sigtools results
`filterINDELs.VAF`|String|minimum variant allele frequency to retain variant
`filterSNVs.VAF`|String|minimum variant allele frequency to retain variant
`hrdResults.tissue`|String|Cancerous-tissue of origin
`hrdResults.sigtoolrScript`|String|Temporary variable to call the .R script containing sigtools, will be modulated. default: ~/sigtools_workflow/sigTools_runthrough.R
`plotResults.plotrScript`|String|Temporary variable to call the .R script containing sigtools, will be modulated. default: ~/sigtools_workflow/sigTools_plotter.R


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`filterStructural.basename`|String|basename("~{structuralVcfFile}",".vcf.gz")|Base name
`filterStructural.modules`|String|"bcftools/1.9"|Required environment modules
`filterStructural.structuralQUALfilter`|String|"'PASS'"|filter for filter calls to keep, eg. PASS
`filterStructural.structuralTYPEfilter`|String|"BND"|filter for tye of structural calls to remove, eg. BND
`filterStructural.jobMemory`|Int|5|Memory allocated for this job (GB)
`filterStructural.threads`|Int|1|Requested CPU threads
`filterStructural.timeout`|Int|1|Hours before task timeout
`filterINDELs.basename`|String|basename("~{smallsVcfFile}",".vcf.gz")|Base name
`filterINDELs.modules`|String|"gatk/4.2.0.0 tabix/1.9 bcftools/1.9 hg38/p12 grch38-alldifficultregions/3.0"|Required environment modules
`filterINDELs.genome`|String|"$HG38_ROOT/hg38_random.fa"|Path to loaded genome .fa
`filterINDELs.difficultRegions`|String?|None|Path to .bed of difficult regions to align to, string must include the --exclude-intervals flag, eg: --exclude-intervals $GRCH38_ALLDIFFICULTREGIONS_ROOT/GRCh38_alldifficultregions.bed
`filterINDELs.QUALfilter`|String|"'PASS,clustered_events'"|filter for filter calls to remove, eg. weak_evidence | strand_bias 
`filterINDELs.jobMemory`|Int|10|Memory allocated for this job (GB)
`filterINDELs.threads`|Int|1|Requested CPU threads
`filterINDELs.timeout`|Int|2|Hours before task timeout
`filterSNVs.basename`|String|basename("~{smallsVcfFile}",".vcf.gz")|Base name
`filterSNVs.modules`|String|"gatk/4.2.0.0 tabix/1.9 bcftools/1.9 hg38/p12 grch38-alldifficultregions/3.0"|Required environment modules
`filterSNVs.genome`|String|"$HG38_ROOT/hg38_random.fa"|Path to loaded genome .fa
`filterSNVs.difficultRegions`|String?|None|Path to .bed of difficult regions to align to, string must include the --exclude-intervals flag, eg: --exclude-intervals $GRCH38_ALLDIFFICULTREGIONS_ROOT/GRCh38_alldifficultregions.bed
`filterSNVs.QUALfilter`|String|"'PASS,clustered_events'"|filter for filter calls to remove, eg. weak_evidence | strand_bias 
`filterSNVs.jobMemory`|Int|10|Memory allocated for this job (GB)
`filterSNVs.threads`|Int|1|Requested CPU threads
`filterSNVs.timeout`|Int|2|Hours before task timeout
`convertSegFile.basename`|String|basename("~{segFile}","_segments.txt")|Base name
`convertSegFile.jobMemory`|Int|5|Memory allocated for this job (GB)
`convertSegFile.threads`|Int|1|Requested CPU threads
`convertSegFile.timeout`|Int|1|Hours before task timeout
`hrdResults.modules`|String|"sigtools/0.0.0.9000"|Required environment modules
`hrdResults.genomeVersion`|String|"hg38"|version of genome, eg hg38
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
  1. write header for file
  
 		echo  -e "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tsample\tsvclass"  >~{basename}.bedpe
 
 2. filter structural bed file for Quality and SVtype, then reformat
 
 		$BCFTOOLS_ROOT/bin/bcftools view -f ~{structuralQUALfilter} ~{structuralVcfFile} |\
 		$BCFTOOLS_ROOT/bin/bcftools filter -e 'INFO/SVTYPE = "~{structuralTYPEfilter}"' |\
 		$BCFTOOLS_ROOT/bin/bcftools query -f "%CHROM\t%POS\t%INFO/END\t%FILTER\t%INFO/SVTYPE\t%INFO/CIPOS\t%INFO/CIEND\n" |\
 		awk -v sampleName=~{sampleName} 'split($6,a,",") split($7,b,",") {print $1"\t"$2+a[1]-1"\t"$2+a[2]"\t"$1"\t"$3+b[1]-1"\t"$2+b[2]"\t"sampleName"\t"$5}' >>~{basename}.bedpe
 
 3. Make filtering report for troubleshooting
 
 		awk '$1 !~ "#" {print}' ~{structuralVcfFile} | wc -l >~{sampleName}.structural.filteringReport.txt
 		awk '$1 !~ "#" {print}' ~{basename}.bedpe | wc -l >>~{sampleName}.structural.filteringReport.txt
 
 	>>>
 	
  ### Filter small variants .vcf
 This step is run twice, for each of SNP and INDEL, under aliases.
 4. filter for type of small variant, and remove troublesome regions of alignment
 
 		gatk SelectVariants \
 		-V ~{smallsVcfFile} \
 		-R ~{genome} ~{difficultRegions} \
 		--select-type-to-include ~{smallType} \
 		-O ~{basename}.~{smallType}.vcf  
 5. filter structural bed file for Quality and VAF, then zip
 
 		$BCFTOOLS_ROOT/bin/bcftools view ~{basename}.~{smallType}.vcf  | \
 		$BCFTOOLS_ROOT/bin/bcftools filter -i "(FORMAT/AD[0:1])/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.~{VAF}" | \
 		$BCFTOOLS_ROOT/bin/bcftools filter -e "~{QUALfilter}" >~{basename}.VAF.vcf
 
 		bgzip ~{basename}.VAF.vcf
 		
 		tabix -p vcf ~{basename}.VAF.vcf.gz
 6. Make filtering report for troubleshooting
 
 		zcat ~{smallsVcfFile} | awk '$1 !~ "#" {print}'  | wc -l >~{sampleName}.filteringReport.txt
 		
 		awk '$1 !~ "#" {print}' ~{basename}.~{smallType}.vcf | wc -l >>~{sampleName}.filteringReport.txt
 		
 		zcat ~{basename}.VAF.vcf.gz | awk '$1 !~ "#" {print}'  | wc -l >>~{sampleName}.filteringReport.txt
 
 	>>>
  ### Convert .seg file to ASCAT format
 
  1. write header for file
 
 		echo  -e "seg_no\tChromosome\tchromStart\tchromEnd\ttotal.copy.number.inNormal\tminor.copy.number.inNormal\ttotal.copy.number.inTumour\tminor.copy.number.inTumour" >~{basename}_segments.cna.txt
 7. reformat for ASCAT
 
 		tail -n +2 ~{segFile} | \
 		awk 'split($1,a,"\"") split(a[2],b,"chr") {print NR"\t"b[2]"\t"$2"\t"$3"\t"2"\t"1"\t"$10"\t"$12}' >>~{basename}_segments.cna.txt
 	>>>
  ### Run the sigtools package in R
 
 		Rscript --vanilla ~{sigtoolrScript} ~{sampleName} ~{tissue} ~{snvVcfFiltered} ~{indelVcfFiltered} ~{structuralBedpeFiltered} ~{lohSegFile} ~{sigtoolsBootstrap} ~{genomeVersion}
 
  ### Run plotting in R (optional)
 
 
 		Rscript --vanilla ~{plotrScript} ~{sampleName} ~{sigTools_hrd_input} ~{sigTools_sigs_input}
 
 
 
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
