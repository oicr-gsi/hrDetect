version 1.0

workflow sigTooler {

	input {
		File structuralVcfFile

    	File smallsVcfFile
    	File smallsVcfIndex

    	File segFile

    	Int variantAlleleFrequency
    	String tissue
	}

	parameter_meta {
    	structuralVcfFile: "Input VCF file of structural variants (eg. from delly)"

    	smallsVcfFile: "Input VCF file of SNV and indels (small mutations) (eg. from mutect2)"
    	smallsVcfIndex: "Index of input VCF file of SNV and indels"

    	segFile: "File for segmentations, used to estimate number of segments in Loss of heterozygosity (LOH) (eg. from sequenza)"

    	variantAlleleFrequency: "Variant Allele Frequency for filtering of small mutations"

    	tissue: "Cancerous-tissue of origin"
    	
  	}

	call filterStructural {
		input: 
			vcfFile = structuralVcfFile,
			outputFileNamePrefix = outputFileNamePrefix
	}

	call filterINDELs {
		input: 
			vcfFile = smallsVcfFile,
			outputFileNamePrefix = outputFileNamePrefix
	}

	call filterSNPs {
		input: 
			vcfFile = smallsVcfFile,
			outputFileNamePrefix = outputFileNamePrefix
	}

	call convertSegFile {
		input: 
			segFile = segFile,
			outputFileNamePrefix = outputFileNamePrefix
	}

	call hrdResults {
		input:	snvVcfFiltered = filterSNVs.snvVcfOutput,
				structuralVcfFiltered = filterStructural.structuralVcfOutput,
				outputFileNamePrefix = outputFileNamePrefix

	}

  	meta {
    	author: "Felix Beaudry"
    	email: "fbeaudry@oicr.on.ca"
    	description: "Homolog Recombination Deficiency Prediction Workflow using sig.tools"
    	dependencies: 
    	[
      		{
        		name: "gatk/4.2.0.0",
        		url: "https://github.com/broadinstitute/gatk/releases"
      		},
      		{
      			name: "tabix/1.9"
      			url: ""
      		},
      		{
      			name: "bcftools"
      			url: ""
      		}

    	]
    	output_meta: {
   			sigToolsOutput : "point estimate and bootstraped confidence intervals for HRD from sigtools"
    	}
	}

	output {
		File sigToolsOutput = "~{basename}.sigtools.hrd.txt"
  	}
  	
}

task filterStructural {
 	input {
	    File vcfFile 
	    String basename = basename("~{vcfFile}", ".vcf.gz")
	    String modules = "bcftools"
	    Int jobMemory = 32
	    Int threads = 4
	    Int timeout = 16
	 }

	parameter_meta {
	    vcfFile: "Vcf input file"
	    basename: "Base name"
	    modules: "Required environment modules"
	    jobMemory: "Memory allocated for this job (GB)"
	    threads: "Requested CPU threads"
	    timeout: "Hours before task timeout"
	}

	command <<<
	    set -euo pipefail

		echo  -e "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tsample\tsvclass"  >~{basename}.filtered.delly.merged.bedpe

		structural_VAF=0

		bcftools query -f "%CHROM\t%POS\t%INFO/END\t%FILTER\t%ALT\t%INFO/CIPOS\t%INFO/CIEND\t[%DR\t]\t[%DV\t]\t[%RR\t]\t[%RV\t]\n" ~{vcfFile} | awk '$5 !~ ":" {print}' | awk '$4 ~ "PASS" {print}' | awk -v VAF=0.~{structural_VAF} '($10+$14)/($8+$10+$12+$14) > VAF {print}' | awk -v sample_t=~{basename} 'split($6,a,",") split($7,b,",") {print $1"\t"$2+a[1]-1"\t"$2+a[2]"\t"$1"\t"$3+b[1]-1"\t"$2+b[2]"\t"sample_t"\t"$5} ' | sed 's/<//g; s/>//g'  >>~{basename}.filtered.delly.merged.bedpe

	>>> 

	runtime {
	    modules: "~{modules}"
	    memory:  "~{jobMemory} GB"
	    cpu:     "~{threads}"
	    timeout: "~{timeout}"
	}

	output {
    	File structuralbedpe = "~{basename}.filtered.delly.merged.bedpe"
	}
	
	meta {
	    	output_meta: {
      		structuralbedpe: "filtered structural .bedpe"
    	}
  	}
}

task filterSNPs {
 	input {
	    File vcfFile 
	    String basename = basename("~{vcfFile}", ".vcf.gz")
	    String modules = "gatk", "tabix", "bcftools", "grch38-alldifficultregions" "hg38/p12"
	    Int jobMemory = 32
	    Int threads = 4
	    Int timeout = 16
	 }

	parameter_meta {
	    vcfFile: "Vcf input file"
	    basename: "Base name"
	    modules: "Required environment modules"
	    jobMemory: "Memory allocated for this job (GB)"
	    threads: "Requested CPU threads"
	    timeout: "Hours before task timeout"
	}

	command <<<
	    set -euo pipefail

	gatk SelectVariants -R ~{HG38_ROOT}/hg38_random.fa  --exclude-intervals ~{GRCH38_ALLDIFFICULTREGIONS_ROOT}/GRCh38_alldifficultregions.bed -V ~{SNV_file}  --select-type-to-include SNP -O ~{sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.SNP.vcf

	bcftools filter -i "(FORMAT/AD[0:1])/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.~{SNP_VAF}" ~{sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.SNP.vcf >~{sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.SNP.MAF~{SNP_VAF}.vcf

	bgzip ~{sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.SNP.MAF*.vcf

	tabix -p vcf ~{sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.SNP.MAF*.vcf.gz

	>>> 

	runtime {
	    modules: "~{modules}"
	    memory:  "~{jobMemory} GB"
	    cpu:     "~{threads}"
	    timeout: "~{timeout}"
	}

	output {
		File indelVcfOutput = "~{basename}.filter.deduped.realigned.recalibrated.mutect2.indels.vcf.gz"
    	File snvVcfOutput = "~{basename}.filter.deduped.realigned.recalibrated.mutect2.SNV.vcf.gz"

	}
	
	meta {
	    	output_meta: {
      		indelVcfOutput: "filtered indel .vcf"
      		snvVcfOutput: "filtered SNV .vcf"

    	}
  	}
}

task filterINDELs {
 	input {
	    File vcfFile 
	    String basename = basename("~{vcfFile}", ".vcf.gz")
	    String modules = "gatk", "tabix", "bcftools", "grch38-alldifficultregions" "hg38/p12"
	    Int jobMemory = 32
	    Int threads = 4
	    Int timeout = 16
	 }

	parameter_meta {
	    vcfFile: "Vcf input file"
	    basename: "Base name"
	    modules: "Required environment modules"
	    jobMemory: "Memory allocated for this job (GB)"
	    threads: "Requested CPU threads"
	    timeout: "Hours before task timeout"
	}

	command <<<
	    set -euo pipefail

	varType="INDEL"

	gatk SelectVariants -R ${HG38_ROOT}/hg38_random.fa  --exclude-intervals ${GRCH38_ALLDIFFICULTREGIONS_ROOT}/GRCh38_alldifficultregions.bed -V ${SNV_file}  --select-type-to-include ${varType} -O ${sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.${varType}.vcf

	bcftools filter -i "(FORMAT/AD[0:1])/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.${indel_VAF}" ${sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.INDEL.vcf >${sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.INDEL.MAF${indel_VAF}.vcf

	bgzip ${sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.${varType}.MAF*.vcf

	tabix -p vcf ${sampleRoot}.filter.deduped.realigned.recalibrated.mutect2.${varType}.MAF*.vcf.gz

	>>> 

	runtime {
	    modules: "~{modules}"
	    memory:  "~{jobMemory} GB"
	    cpu:     "~{threads}"
	    timeout: "~{timeout}"
	}

	output {
		File indelVcfOutput = "~{basename}.filter.deduped.realigned.recalibrated.mutect2.indels.vcf.gz"
    	File snvVcfOutput = "~{basename}.filter.deduped.realigned.recalibrated.mutect2.SNV.vcf.gz"

	}
	
	meta {
	    	output_meta: {
      		indelVcfOutput: "filtered indel .vcf"
      		snvVcfOutput: "filtered SNV .vcf"

    	}
  	}
}

task convertSegFile {
 	input {
	    File segFile 
	    String basename = basename("~{segFile}", "_segments.txt")
	    Int jobMemory = 32
	    Int threads = 4
	    Int timeout = 16
	 }

	parameter_meta {
	    segFile: "segment input file from sequenza"
	    basename: "Base name"
	    jobMemory: "Memory allocated for this job (GB)"
	    threads: "Requested CPU threads"
	    timeout: "Hours before task timeout"
	}

	command <<<

	    set -euo pipefail

		gamma=$(awk 'split($1,a,"=") $1 ~ "sequenza_gamma" {print a[2]}'  ${studyLocation}/${study}/${sampleRoot}/config.ini )

		echo  -e "seg_no\tChromosome\tchromStart\tchromEnd\ttotal.copy.number.inNormal\tminor.copy.number.inNormal\ttotal.copy.number.inTumour\tminor.copy.number.inTumour"  >${wrkdir}/${sampleRoot}_segments.cna.txt
		
		tail -n +2 ${studyLocation}/${study}/${sampleRoot}/gammas/${gamma}/${sampleRoot}*_segments.txt | awk 'split($1,a,"\"") split(a[2],b,"chr") {print NR"\t"b[2]"\t"$2"\t"$3"\t"2"\t"1"\t"$10"\t"$12}'  >>${wrkdir}/${sampleRoot}_segments.cna.txt
	>>> 

	runtime {
	    modules: "~{modules}"
	    memory:  "~{jobMemory} GB"
	    cpu:     "~{threads}"
	    timeout: "~{timeout}"
	}

	output {
		File segmentsOutput = "~{basename}_segments.cna.txt"

	}
	
	meta {
	    	output_meta: {
      		segmentsOutput: "reformatted segmentation file"

    	}
  	}
}

task hrdResults {
	input {
		file basename
		String = ~{tissue}
		String modules = "sigtools"
		Int jobMemory = 32
		Int threads = 4
		Int timeout = 16
	}

	parameter_meta {
		snvVcfFiltered: "filtered SNV Vcf input file"
		structuralVcfFiltered: "filtered structural variant Vcf input file"
		modules: "Required environment modules"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail

		Rscript --vanilla ~/sigtools_workflow/sigTools_runthrough.R ~{sampleRoot} ~{HRDtissue} ~{snvFile_loc} ~{indel_vcf_file} ~{SV_bedpe_file} ~{LOH_seg_file}

	>>> 

	runtime {
	    modules: "~{modules}"
	    memory:  "~{jobMemory} GB"
	    cpu:     "~{threads}"
	    timeout: "~{timeout}"
	}

	output {
			File sigToolsOutput = "~{basename}.sigtools.hrd.txt"
	}

	meta {
    	output_meta: {
      		sigToolsOutput : "point estimate and bootstraped confidence intervals for HRD from sigtools"
	    }
	}
}





