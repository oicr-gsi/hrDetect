workflow sigTooler {

	input {
		File structuralVcfFile
    	File structuralVcfIndex

    	File smallsVcfFile
    	File smallsVcfIndex

    	File segFile

    	String VAF
    	String tissue
	}

	call filterStructural {
		input: vcfFile = structuralVcfFile
	}

	call filterSmalls {
		input: vcfFile = smallsVcfFile
	}

	call convertSegFile {
		input: segFile = segFile
	}

	call hrdResults {
		input:	snvVcfFiltered = filterSNVs.snvVcfOutput
				structuralVcfFiltered = filterStructural.structuralVcfOutput
	}

	parameter_meta {
    	smallsVcfFile: "Input VCF file"
    	smallsVcfIndex: "Input VCF index file"
    	targetBed: "Target bed file"
  	}

  	meta {
    	author: "Felix Beaudry"
    	email: "fbeaudry@oicr.on.ca"
    	description: "Homolog Recombination Deficiency Prediction Workflow using sig.tools"
    	dependencies: 
    	[
    		{
      			name: "svtools"
      			url: "https://github.com/ctsa/svtools"
      		},
      		{
        		name: "gatk/4.2.0.0",
        		url: "https://github.com/broadinstitute/gatk/releases"
      		},
      		{
      			name: "tabix/1.9"
      			url: ""
      		},
      		{
        		name: "rstats/4.1.2",
        		url: "https://cran.r-project.org/mirrors.html"
      		}
      		{
        		name: "svtools/0.5.1",
        		url: "https://github.com/hall-lab/svtools/archive/refs/tags/v0.5.1.tar.gz"
      		}
    	]
    	output_meta: {
   			
    	}
	}

	output {
		

  }
}

task filterStructural {
 	input {
	    File vcfFile 
	    String basename = basename("~{vcfFile}", ".vcf.gz")
	    String modules = "svtools"
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

		vcfToBedpe -i ~{structuralVcfFile} -o ~{basename}.somatic.delly.merged.bedpe

		echo  -e "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tsample\tsvclass"  >~{basename}.somatic.delly.merged.reformat.bedpe
		awk -v basename=~{basename} '$12 ~ "PASS" {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"basename"\t"$11}' ~{basename}.somatic.delly.merged.bedpe >>~{basename}.somatic.delly.merged.reformat.bedpe

	>>> 

	runtime {
	    modules: "~{modules}"
	    memory:  "~{jobMemory} GB"
	    cpu:     "~{threads}"
	    timeout: "~{timeout}"
	}

	output {
    	File structuralbedpe = "~{basename}.somatic.delly.merged.reformat.bedpe"
	}
	
	meta {
	    	output_meta: {
      		structuralbedpe: "filtered structural .bedpe"
    	}
  	}
}

task filterSmalls {
 	input {
	    File vcfFile 
	    File targetBed
	    String basename = basename("~{vcfFile}", ".vcf.gz")
	    String modules = "svtools"
	   	String referenceFasta
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

		module load gatk tabix

		gatk SelectVariants -R ~{referenceFasta} -V ~{vcfFile} --exclude-intervals ~{targetBed} -O ~{basename}.filter.deduped.realigned.recalibrated.mutect2.indels.vcf --select-type-to-include INDEL
		
		bgzip ~{basename}.filter.deduped.realigned.recalibrated.mutect2.indels.vcf
		
		tabix -p vcf ~{basename}.filter.deduped.realigned.recalibrated.mutect2.indels.vcf.gz

		gatk SelectVariants -R ~{referenceFasta} -V ~{vcfFile} --exclude-intervals ~{targetBed} -O ~{basename}.filter.deduped.realigned.recalibrated.mutect2.indels.vcf --select-type-to-include SNP
		
		bgzip ~{basename}.filter.deduped.realigned.recalibrated.mutect2.SNV.vcf
		
		tabix -p vcf ~{basename}.filter.deduped.realigned.recalibrated.mutect2.SNV.vcf.gz

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

		echo  -e "seg_no\tChromosome\tchromStart\tchromEnd\ttotal.copy.number.inNormal\tminor.copy.number.inNormal\ttotal.copy.number.inTumour\tminor.copy.number.inTumour"  >~{basename}_segments.cna.txt

		tail -n +2 segFile | awk 'split($1,a,"\"") split(a[2],b,"chr") {print NR"\t"b[2]"\t"$2"\t"$3"\t"2"\t"1"\t"$10"\t"$12}' >>~{basename}_segments.cna.txt

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
		String modules = "R"
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

		Rscript --vanilla sigTools_runthrough.R ~{basename} ~{tissue}

	>>> 

	runtime {
	    modules: "~{modules}"
	    memory:  "~{jobMemory} GB"
	    cpu:     "~{threads}"
	    timeout: "~{timeout}"
	}

	output {
	
	}

	meta {
    	output_meta: {

	    }
	}
}





