# hrDetect

Homologous Recombination Deficiency (HRD) Prediction Workflow using sig.tools

## Overview

## Dependencies

* [tabix 1.9](http://www.htslib.org/doc/tabix.html)
* [bcftools 1.9](https://samtools.github.io/bcftools/bcftools.html)
* [sigtools 0.0.0.9000](https://github.com/Nik-Zainal-Group/signature.tools.lib)


## Usage

### Cromwell
```
java -jar cromwell.jar run hrDetect.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`outputFileNamePrefix`|String|Name of sample matching the tumor sample in .vcf
`structuralVcfFile`|File|Input VCF file of structural variants (eg. from delly)
`smallsVcfFile`|File|Input VCF file of SNV and indels (small mutations) (eg. from mutect2)
`smallsVcfIndex`|File|Index file for smallsVcfFile
`segFile`|File|File for segmentations, used to estimate number of segments in Loss of heterozygosity (LOH) (eg. from sequenza)
`reference`|String|Reference genome version


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`filterStructural.modules`|String|"bcftools/1.9"|Required environment modules
`filterStructural.structuralQUALfilter`|String|"PASS"|filter for filter calls to keep, eg. PASS
`filterStructural.jobMemory`|Int|5|Memory allocated for this job (GB)
`filterStructural.threads`|Int|1|Requested CPU threads
`filterStructural.timeout`|Int|1|Hours before task timeout
`filterINDELs.VAF`|Float|0.01|minimum variant allele frequency to retain variant
`filterINDELs.QUALfilter`|String|"FILTER~'haplotype' | FILTER~'clustered_events' | FILTER~'multiallelic' | FILTER~'slippage' |FILTER~'weak_evidence' | FILTER~'strand_bias' | FILTER~'position' | FILTER~'normal_artifact' |  FILTER~'map_qual' | FILTER~'germline' | FILTER~'fragment' | FILTER~'contamination' | FILTER~'base_qual'"|filter for filter calls to remove, eg. FILTER~'weak_evidence' | FILTER~'strand_bias' 
`filterINDELs.jobMemory`|Int|10|Memory allocated for this job (GB)
`filterINDELs.threads`|Int|1|Requested CPU threads
`filterINDELs.timeout`|Int|2|Hours before task timeout
`filterSNVs.VAF`|Float|0.01|minimum variant allele frequency to retain variant
`filterSNVs.QUALfilter`|String|"FILTER~'haplotype' | FILTER~'clustered_events' | FILTER~'multiallelic' | FILTER~'slippage' |FILTER~'weak_evidence' | FILTER~'strand_bias' | FILTER~'position' | FILTER~'normal_artifact' |  FILTER~'map_qual' | FILTER~'germline' | FILTER~'fragment' | FILTER~'contamination' | FILTER~'base_qual'"|filter for filter calls to remove, eg. FILTER~'weak_evidence' | FILTER~'strand_bias' 
`filterSNVs.jobMemory`|Int|10|Memory allocated for this job (GB)
`filterSNVs.threads`|Int|1|Requested CPU threads
`filterSNVs.timeout`|Int|2|Hours before task timeout
`hrdResults.modules`|String|"sigtools/2.4.1 sigtools-data/1.0 hrdetect-rscript/1.5.8"|Required environment modules
`hrdResults.sigtoolrScript`|String|"$HRDETECT_RSCRIPT_ROOT/scripts/sigTools_runthrough.R"|.R script containing sigtools
`hrdResults.SVrefSigs`|String|"$SIGTOOLS_DATA_ROOT/RefSigv0_Rearr.tsv"|reference signatures for SVs
`hrdResults.SNVrefSigs`|String|"$SIGTOOLS_DATA_ROOT/COSMIC_v1_SBS_GRCh38.txt"|reference signatures for SNVs
`hrdResults.sigtoolsBootstrap`|Int|200|Number of bootstraps for sigtools
`hrdResults.indelCutoff`|Int|50|minimum number of indels to run analysis
`hrdResults.jobMemory`|Int|50|Memory allocated for this job (GB)
`hrdResults.threads`|Int|1|Requested CPU threads
`hrdResults.timeout`|Int|15|Hours before task timeout


### Outputs

Output | Type | Description | Labels
---|---|---|---
`hrd_signatures`|File|JSON file of hrdetect signatures|vidarr_label: hrd_signatures 
`SBS_exposures`|File|JSON of single basepair substitution signatures|
`SV_exposures`|File|JSON of structural variant signatures|vidarr_label: SV_exposures
`ID_catalog`|File|JSON cataloguing indels|vidarr_label: ID_catalog


This section lists command(s) run by hrDetect workflow
  
* Running hrDetect

### Generate list if PASS calls
  
```
	set -euo pipefail


	$BCFTOOLS_ROOT/bin/bcftools view -f '~{structuralQUALfilter}' ~{structuralVcfFile} >> ~{outputFileNamePrefix}.structural.PASS.vcf

	awk '$1 !~ "#" {print}' ~{structuralVcfFile} | wc -l >~{outputFileNamePrefix}.structural.filteringReport.txt
	awk '$1 !~ "#" {print}' ~{outputFileNamePrefix}.structural.PASS.vcf | wc -l >>~{outputFileNamePrefix}.structural.filteringReport.txt

```

### Normalize and Filter calls

```
	set -euo pipefail

	$BCFTOOLS_ROOT/bin/bcftools norm --multiallelics - --fasta-ref ~{genome} ~{difficultRegions} ~{smallsVcfFile} | \
	$BCFTOOLS_ROOT/bin/bcftools filter -i "TYPE='~{smallType}'" | \
	$BCFTOOLS_ROOT/bin/bcftools filter -e "~{QUALfilter}" | \
	$BCFTOOLS_ROOT/bin/bcftools filter -i "(FORMAT/AD[0:1])/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= ~{VAF}" >~{outputFileNamePrefix}.~{smallType}.VAF.vcf

	bgzip ~{outputFileNamePrefix}.~{smallType}.VAF.vcf
	tabix -p vcf ~{outputFileNamePrefix}.~{smallType}.VAF.vcf.gz

	zcat ~{smallsVcfFile} | awk '$1 !~ "#" {print}'  | wc -l >~{outputFileNamePrefix}.~{smallType}.filteringReport.txt
	zcat ~{outputFileNamePrefix}.~{smallType}.VAF.vcf.gz | awk '$1 !~ "#" {print}'  | wc -l >>~{outputFileNamePrefix}.~{smallType}.filteringReport.txt
  
```

### Generate report

```
	set -euo pipefail

	Rscript ~{sigtoolrScript} \
		--sampleName ~{outputFileNamePrefix} \
		--snvFile ~{snvVcfFiltered} \
		--indelFile  ~{indelVcfFiltered} \
		--SVFile ~{SV_vcf_location} \
		--LOHFile ~{lohSegFile} \
		--bootstraps ~{sigtoolsBootstrap} \
		--genomeVersion ~{genomeVersion} \
		--indelCutoff ~{indelCutoff} \
		--SVrefSigs ~{SVrefSigs} \
		--SNVrefSigs ~{SNVrefSigs}
```
  
## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
