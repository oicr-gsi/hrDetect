# sigTooler

Homologous Recombination Deficiency (HRD) Prediction Workflow using sig.tools

## Overview

## Dependencies

* [gatk 4.2.0.0](https://github.com/broadinstitute/gatk/releases)
* [tabix 1.9](http://www.htslib.org/doc/tabix.html)
* [bcftools 1.9](https://samtools.github.io/bcftools/bcftools.html)
* [sigtools 0.0.0.9000](https://github.com/Nik-Zainal-Group/signature.tools.lib)
* [grch38-alldifficultregions 3.0](https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/genome-stratifications/v3.0/GRCh38/union/)
* [hg38 p12](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.38/)


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
`indelVAF`|Int|Variant Allele Frequency for filtering of indel mutations
`snvVAF`|Int|Variant Allele Frequency for filtering of SNVs
`tissue`|String|Cancerous-tissue of origin
`rScript`|String|Temporary variable to call the .R script containing sigtools, will be modulated. default: ~/sigtools_workflow/sigTools_runthrough.R
`sampleName`|String|Name of sample matching the tumor sample in .vcf


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`filterStructural.basename`|String|basename("~{structuralVcfFile}",".vcf.gz")|Base name
`filterStructural.modules`|String|"bcftools/1.9"|Required environment modules
`filterStructural.structuralVAF`|Int|0|VAF for structural variants
`filterStructural.jobMemory`|Int|5|Memory allocated for this job (GB)
`filterStructural.threads`|Int|1|Requested CPU threads
`filterStructural.timeout`|Int|1|Hours before task timeout
`filterINDELs.basename`|String|basename("~{smallsVcfFile}",".vcf.gz")|Base name
`filterINDELs.modules`|String|"gatk/4.2.0.0 tabix/1.9 bcftools/1.9 grch38-alldifficultregions/3.0 hg38/p12"|Required environment modules
`filterINDELs.genome`|String|"$HG38_ROOT/hg38_random.fa"|Path to loaded genome
`filterINDELs.difficultRegions`|String|"$GRCH38_ALLDIFFICULTREGIONS_ROOT/GRCh38_alldifficultregions.bed"|Path to loaded difficult regions to align to
`filterINDELs.jobMemory`|Int|10|Memory allocated for this job (GB)
`filterINDELs.threads`|Int|1|Requested CPU threads
`filterINDELs.timeout`|Int|2|Hours before task timeout
`filterSNVs.basename`|String|basename("~{smallsVcfFile}",".vcf.gz")|Base name
`filterSNVs.modules`|String|"gatk/4.2.0.0 tabix/1.9 bcftools/1.9 grch38-alldifficultregions/3.0 hg38/p12"|Required environment modules
`filterSNVs.genome`|String|"$HG38_ROOT/hg38_random.fa"|Path to loaded genome
`filterSNVs.difficultRegions`|String|"$GRCH38_ALLDIFFICULTREGIONS_ROOT/GRCh38_alldifficultregions.bed"|Path to loaded difficult regions to align to
`filterSNVs.jobMemory`|Int|10|Memory allocated for this job (GB)
`filterSNVs.threads`|Int|1|Requested CPU threads
`filterSNVs.timeout`|Int|2|Hours before task timeout
`convertSegFile.basename`|String|basename("~{segFile}","_segments.txt")|Base name
`convertSegFile.jobMemory`|Int|5|Memory allocated for this job (GB)
`convertSegFile.threads`|Int|1|Requested CPU threads
`convertSegFile.timeout`|Int|1|Hours before task timeout
`hrdResults.modules`|String|"sigtools/0.0.0.9000"|Required environment modules
`hrdResults.jobMemory`|Int|20|Memory allocated for this job (GB)
`hrdResults.threads`|Int|1|Requested CPU threads
`hrdResults.timeout`|Int|2|Hours before task timeout


### Outputs

Output | Type | Description
---|---|---
`sigToolsOutput`|File|point estimate and bootstraped confidence intervals for HRD from sigtools


## Commands
 This section lists commands run by sigtools_workflow. It launches four tasks to make input files for sigtools, and then runs the sigtools packages in an .R script.  

### Convert structural variant .vcf to .bedpe 

 
        echo  -e "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tsample\tsvclass"  >~{basename}.bedpe
 
        $BCFTOOLS_ROOT/bin/bcftools query -f "%CHROM\t%POS\t%INFO/END\t%FILTER\t%ALT\t%INFO/CIPOS\t%INFO/CIEND\t[%DR\t]\t[%DV\t]\t[%RR\t]\t[%RV\t]\n" ~{structuralVcfFile} | \
        awk '$5 !~ ":" {print}' | \
        awk '$4 ~ "PASS" {print}' | \
        awk -v VAF=0.~{structuralVAF} '($10+$14)/($8+$10+$12+$14) > VAF {print}' | \
        awk -v sampleName=~{sampleName} 'split($6,a,",") split($7,b,",") {print $1"\t"$2+a[1]-1"\t"$2+a[2]"\t"$1"\t"$3+b[1]-1"\t"$2+b[2]"\t"sampleName"\t"$5} ' | \
        sed 's/<//g; s/>//g' >>~{basename}.bedpe
 
### Pull INDELs from small variants .vcf
 
        gatk SelectVariants \
        -V ~{smallsVcfFile} \
        -R ~{genome}  \
        --exclude-intervals ~{difficultRegions} \
        --select-type-to-include INDEL \
        -O ~{basename}.INDEL.vcf
 
        $BCFTOOLS_ROOT/bin/bcftools filter -i "(FORMAT/AD[0:1])/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.~{indelVAF}" ~{basename}.INDEL.vcf >~{basename}.INDEL.VAF.vcf
 
        bgzip ~{basename}.INDEL.VAF.vcf
 
        tabix -p vcf ~{basename}.INDEL.VAF.vcf.gz

### Pull SNPs from small variants .vcf
 
        gatk SelectVariants \
        -V ~{smallsVcfFile} \
        -R ~{genome}  \
        --exclude-intervals ~{difficultRegions} \
        --select-type-to-include SNP \
        -O ~{basename}.SNP.vcf
 
        $BCFTOOLS_ROOT/bin/bcftools filter -i "(FORMAT/AD[0:1])/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.~{snvVAF}" ~{basename}.SNP.vcf >~{basename}.SNP.VAF.vcf
 
        bgzip ~{basename}.SNP.VAF.vcf
 
        tabix -p vcf ~{basename}.SNP.VAF.vcf.gz
 
    
### Convert .seg file to ASCAT format
 
        set -euo pipefail
 
        echo  -e "seg_no\tChromosome\tchromStart\tchromEnd\ttotal.copy.number.inNormal\tminor.copy.number.inNormal\ttotal.copy.number.inTumour\tminor.copy.number.inTumour" >~{basename}_segments.cna.txt
 
        tail -n +2 ~{segFile} | \
        awk 'split($1,a,"\"") split(a[2],b,"chr") {print NR"\t"b[2]"\t"$2"\t"$3"\t"2"\t"1"\t"$10"\t"$12}' >>~{basename}_segments.cna.txt

### Run the sigtools package in R
        Rscript --vanilla ~{rScript} ~{sampleName} ~{tissue} ~{snvVcfFiltered} ~{indelVcfFiltered} ~{structuralBedpeFiltered} ~{lohSegFile}

 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
