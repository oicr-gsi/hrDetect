version 1.0

workflow sigTooler {

	input {
		File structuralVcfFile
    	File smallsVcfFile
    	File smallsVcfIndex
    	File segFile
    	Int structuralVAF
    	Int indelVAF
    	Int snvVAF
    	String tissue
    	String rScript
    	String sampleName
	}

	parameter_meta {
    	structuralVcfFile: "Input VCF file of structural variants (eg. from delly)"
    	smallsVcfFile: "Input VCF file of SNV and indels (small mutations) (eg. from mutect2)"
    	smallsVcfIndex: "Index of input VCF file of SNV and indels"
    	segFile: "File for segmentations, used to estimate number of segments in Loss of heterozygosity (LOH) (eg. from sequenza)"
    	structuralVAF: "Variant Allele Frequency for filtering of structural variants. default: 0"
    	indelVAF: "Variant Allele Frequency for filtering of indel mutations"
    	snvVAF: "Variant Allele Frequency for filtering of SNVs"
    	tissue: "Cancerous-tissue of origin"
    	rScript: "Temporary variable to call the .R script containing sigtools, will be modulated. default: ~/sigtools_workflow/sigTools_runthrough.R"
    	sampleName: "Name of sample matching the tumor sample in .vcf"
  	}

	call filterStructural {
		input: 
			structuralVcfFile = structuralVcfFile,
			structuralVAF = structuralVAF,
			sampleName = sampleName
	}

	call filterINDELs {
		input: 
			smallsVcfFile = smallsVcfFile,
			indelVAF = indelVAF
	}

	call filterSNVs {
		input: 
			smallsVcfFile = smallsVcfFile,
			snvVAF = snvVAF
	}

	call convertSegFile {
		input: 
			segFile = segFile
	}

	call hrdResults {
		input:
			structuralBedpeFiltered = filterStructural.structuralbedpe,
			indelVcfFiltered = filterINDELs.indelVcfOutput,
			snvVcfFiltered = filterSNVs.snvVcfOutput,
			lohSegFile = convertSegFile.segmentsOutput,
			tissue = tissue,
    		rScript = rScript,
    		sampleName = sampleName
	}

  	meta {
    	author: "Felix Beaudry"
    	email: "fbeaudry@oicr.on.ca"
    	description: "Homologous Recombination Deficiency (HRD) Prediction Workflow using sig.tools"
    	dependencies: 
    	[
      		{
        		name: "gatk/4.2.0.0",
        		url: "https://github.com/broadinstitute/gatk/releases"
      		},
      		{
      			name: "tabix/1.9",
      			url: "http://www.htslib.org/doc/tabix.html"
      		},
      		{
      			name: "bcftools/1.9",
      			url: "https://samtools.github.io/bcftools/bcftools.html"
      		},
      		{
      			name: "sigtools/0.0.0.9000",
      			url: "https://github.com/Nik-Zainal-Group/signature.tools.lib"
      		},
      		{
      			name: "grch38-alldifficultregions/3.0",
      			url: "https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/genome-stratifications/v3.0/GRCh38/union/"
      		},
      		{
      			name: "hg38/p12",
      			url: "https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.38/"
      		}
    	]
    	output_meta: {
   			sigToolsOutput : "point estimate and bootstraped confidence intervals for HRD from sigtools"
    	}
	}

	output {
		File sigToolsOutput = "~{sampleName}.sigtools.hrd.txt"
  	}
}

task filterStructural {
	input {
		File structuralVcfFile 
		String basename = basename("~{structuralVcfFile}", ".vcf.gz")
		String modules = "bcftools/1.9"
		String sampleName
		Int structuralVAF = 0
		Int jobMemory = 5
		Int threads = 1
		Int timeout = 1
	}

	parameter_meta {
		structuralVcfFile: "Vcf input file"
		basename: "Base name"
		modules: "Required environment modules"
		sampleName: "Name of sample matching the tumor sample in .vcf"
		structuralVAF: "VAF for structural variants"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail

		echo  -e "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tsample\tsvclass"  >~{basename}.bedpe

		bcftools query -f "%CHROM\t%POS\t%INFO/END\t%FILTER\t%ALT\t%INFO/CIPOS\t%INFO/CIEND\t[%DR\t]\t[%DV\t]\t[%RR\t]\t[%RV\t]\n" ~{structuralVcfFile} | \
		awk '$5 !~ ":" {print}' | \
		awk '$4 ~ "PASS" {print}' | \
		awk -v VAF=0.~{structuralVAF} '($10+$14)/($8+$10+$12+$14) > VAF {print}' | \
		awk -v sampleName=~{sampleName} 'split($6,a,",") split($7,b,",") {print $1"\t"$2+a[1]-1"\t"$2+a[2]"\t"$1"\t"$3+b[1]-1"\t"$2+b[2]"\t"sampleName"\t"$5} ' | \
		sed 's/<//g; s/>//g' >>~{basename}.bedpe
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

task filterINDELs {
	input {
		File smallsVcfFile
		String basename = basename("~{smallsVcfFile}", ".vcf.gz")
		String modules = "gatk/4.2.0.0 tabix/1.9 bcftools/1.9 grch38-alldifficultregions/3.0 hg38/p12"
		String genome = "$HG38_ROOT/hg38_random.fa"
		String difficultRegions = "$GRCH38_ALLDIFFICULTREGIONS_ROOT}/GRCh38_alldifficultregions.bed"
		Int indelVAF
		Int jobMemory = 10
		Int threads = 1
		Int timeout = 2
	}

	parameter_meta {
		smallsVcfFile: "Vcf input file"
		basename: "Base name"
		modules: "Required environment modules"
		genome: "Path to loaded genome"
		difficultRegions: "Path to loaded difficult regions to align to"
		indelVAF: "VAF for indels"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail

		gatk SelectVariants \
		-V ~{smallsVcfFile} \
		-R ~{genome}  \
		--exclude-intervals ~{difficultRegions} \
		--select-type-to-include INDEL \
		-O ~{basename}.INDEL.vcf

		bcftools filter -i "(FORMAT/AD[0:1])/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.~{indelVAF}" ~{basename}.INDEL.vcf >~{basename}.INDEL.VAF.vcf

		bgzip ~{basename}.INDEL.VAF.vcf

		tabix -p vcf ~{basename}.INDEL.VAF.vcf.gz
	>>> 

	runtime {
		modules: "~{modules}"
		memory: "~{jobMemory} GB"
		cpu: "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File indelVcfOutput = "~{basename}.INDEL.VAF.vcf.gz"
	}
	
	meta {
		output_meta: {
			indelVcfOutput: "filtered indel .vcf"
		}
	}
}

task filterSNVs {
	input {
		File smallsVcfFile
		String basename = basename("~{smallsVcfFile}", ".vcf.gz")
		String modules = "gatk/4.2.0.0 tabix/1.9 bcftools/1.9 grch38-alldifficultregions/3.0 hg38/p12"
		String genome = "$HG38_ROOT/hg38_random.fa"
		String difficultRegions = "$GRCH38_ALLDIFFICULTREGIONS_ROOT}/GRCh38_alldifficultregions.bed"
		Int snvVAF
		Int jobMemory = 10
		Int threads = 1
		Int timeout = 2
	}

	parameter_meta {
		smallsVcfFile: "Vcf input file"
		snvVAF: "VAF for SNV filtering"
		basename: "Base name"
		modules: "Required environment modules"
		genome: "Path to loaded genome"
		difficultRegions: "Path to loaded difficult regions to align to"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail

		gatk SelectVariants \
		-V ~{smallsVcfFile} \
		-R ~{genome}  \
		--exclude-intervals ~{difficultRegions} \
		--select-type-to-include SNP \
		-O ~{basename}.SNP.vcf

		bcftools filter -i "(FORMAT/AD[0:1])/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.~{snvVAF}" ~{basename}.SNP.vcf >~{basename}.SNP.VAF.vcf

		bgzip ~{basename}.SNP.VAF.vcf

		tabix -p vcf ~{basename}.SNP.VAF.vcf.gz

	>>> 

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File snvVcfOutput = "~{basename}.SNP.VAF.vcf.gz"
	}

	meta {
		output_meta: {
			snvVcfOutput: "filtered SNV .vcf"
		}
	}
}

task convertSegFile {
	input {
		File segFile 
		String basename = basename("~{segFile}", "_segments.txt")
		Int jobMemory = 5
		Int threads = 1
		Int timeout = 1
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

		echo  -e "seg_no\tChromosome\tchromStart\tchromEnd\ttotal.copy.number.inNormal\tminor.copy.number.inNormal\ttotal.copy.number.inTumour\tminor.copy.number.inTumour" >~{basename}_segments.cna.txt

		tail -n +2 ~{segFile} | \
		awk 'split($1,a,"\"") split(a[2],b,"chr") {print NR"\t"b[2]"\t"$2"\t"$3"\t"2"\t"1"\t"$10"\t"$12}' >>~{basename}_segments.cna.txt
	>>> 

	runtime {
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
		File structuralBedpeFiltered
		File indelVcfFiltered
		File snvVcfFiltered
		File lohSegFile
		String tissue
		String rScript
		String sampleName
		String modules = "sigtools/0.0.0.9000"
		Int jobMemory = 20
		Int threads = 1
		Int timeout = 2
	}

	parameter_meta {
		structuralBedpeFiltered: "filtered structural variant .bedpe"
		indelVcfFiltered: "filtered INDEL .vcf"
		snvVcfFiltered: "filtered SNV .vcf"
		lohSegFile: "reformatted segmentation file"
		tissue: "Cancerous-tissue of origin"
		rScript: "Temporary variable to call the .R script containing sigtools, will be modulated. default: ~/sigtools_workflow/sigTools_runthrough.R"
		sampleName: "Name of sample matching the tumor sample in .vcf"		
		modules: "Required environment modules"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail

		Rscript --vanilla ~{rScript} ~{sampleName} ~{tissue} ~{snvVcfFiltered} ~{indelVcfFiltered} ~{structuralBedpeFiltered} ~{lohSegFile}
	>>> 

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File sigToolsOutput = "~{sampleName}.sigtools.hrd.txt"
	}

	meta {
		output_meta: {
			sigToolsOutput : "point estimate and bootstraped confidence intervals for HRD from sigtools"
		}
	}
}