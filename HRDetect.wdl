version 1.0

workflow HRDetect {
	input {
		File structuralVcfFile
		File smallsVcfFile
		File smallsVcfIndex
		File segFile
		String sampleName
	}

	parameter_meta {
		structuralVcfFile: "Input VCF file of structural variants (eg. from delly)"
		smallsVcfFile: "Input VCF file of SNV and indels (small mutations) (eg. from mutect2)"
		smallsVcfIndex: "Index of input VCF file of SNV and indels"
		segFile: "File for segmentations, used to estimate number of segments in Loss of heterozygosity (LOH) (eg. from sequenza)"
		sampleName: "Name of sample matching the tumor sample in .vcf"
	}

	call filterStructural {
		input: 
			structuralVcfFile = structuralVcfFile,
			sampleName = sampleName
	}

	call filterSMALLs as filterINDELs {
		input: 
			smallsVcfFile = smallsVcfFile,
			smallsVcfIndex = smallsVcfIndex,
			sampleName = sampleName,
			smallType = "indel"
	}

	call filterSMALLs as filterSNVs {
		input: 
			smallsVcfFile = smallsVcfFile,
			smallsVcfIndex = smallsVcfIndex,
			sampleName = sampleName,
			smallType = "snp"
	}

	call convertSegFile {
		input: 
			segFile = segFile
	}

	call hrdResults {
		input:
			structuralBedpeFiltered = filterStructural.structuralbedpe,
			indelVcfFiltered = filterINDELs.smallsVcfOutput,
			indelVcfIndexFiltered = filterINDELs.smallsVcfIndexOutput,
			snvVcfFiltered = filterSNVs.smallsVcfOutput,
			snvVcfIndexFiltered = filterSNVs.smallsVcfIndexOutput,
			lohSegFile = convertSegFile.segmentsOutput,
			sampleName = sampleName
	}

	meta {
		author: "Felix Beaudry"
		email: "fbeaudry@oicr.on.ca"
		description: "Homologous Recombination Deficiency (HRD) Prediction Workflow using sig.tools"
		dependencies: 
		[
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
			}
		]
		output_meta: {
			JSONout : "sigtools and CHORD results in JSON",
			indelFilteringReport: "counts of INDELs pre and post filtering",
			snvFilteringReport: "counts of SNVs pre and post filtering",
			structuralFilteringReport: "counts of structural variants pre and post filtering"
		}
	}
	output {
		File indelFilteringReport = "~{sampleName}.INDEL.filteringReport.txt"
		File snvFilteringReport = "~{sampleName}.SNP.filteringReport.txt"
		File structuralFilteringReport = "~{sampleName}.structural.filteringReport.txt"
		File JSONout = "~{sampleName}.signatures.json"
	}
}

task filterStructural {
	input {
		File structuralVcfFile 
		String basename = basename("~{structuralVcfFile}", ".vcf.gz")
		String modules = "bcftools/1.9"
		String sampleName
		String structuralQUALfilter = "PASS"
		String structuralTYPEfilter = "BND"
		Int jobMemory = 5
		Int threads = 1
		Int timeout = 1
	}

	parameter_meta {
		structuralVcfFile: "Vcf input file"
		basename: "Base name"
		modules: "Required environment modules"
		sampleName: "Name of sample matching the tumor sample in .vcf"
		structuralQUALfilter: "filter for filter calls to keep, eg. PASS"
		structuralTYPEfilter: "filter for tye of structural calls to remove, eg. BND"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail

		echo  -e "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tsample\tsvclass"  >~{basename}.bedpe

		$BCFTOOLS_ROOT/bin/bcftools view -f '~{structuralQUALfilter}' ~{structuralVcfFile} |\
		$BCFTOOLS_ROOT/bin/bcftools filter -e 'INFO/SVTYPE = "~{structuralTYPEfilter}"' |\
		$BCFTOOLS_ROOT/bin/bcftools query -f "%CHROM\t%POS\t%INFO/END\t%FILTER\t%INFO/SVTYPE\t%INFO/CIPOS\t%INFO/CIEND\n" |\
		awk -v sampleName=~{sampleName} 'split($6,a,",") split($7,b,",") {print $1"\t"$2+a[1]-1"\t"$2+a[2]"\t"$1"\t"$3+b[1]-1"\t"$3+b[2]"\t"sampleName"\t"$5}' >>~{basename}.bedpe

		awk '$1 !~ "#" {print}' ~{structuralVcfFile} | wc -l >~{sampleName}.structural.filteringReport.txt
		awk '$1 !~ "#" {print}' ~{basename}.bedpe | wc -l >>~{sampleName}.structural.filteringReport.txt

	>>>

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File structuralbedpe = "~{basename}.bedpe"
		File structuralFilteringReport = "~{sampleName}.structural.filteringReport.txt"
	}

	meta {
		output_meta: {
			structuralbedpe: "filtered structural .bedpe",
			structuralFilteringReport: "counts of variants pre and post filtering"
		}
	}
}

task filterSMALLs {
	input {
		String sampleName
		String smallType
		File smallsVcfFile
		File smallsVcfIndex
		String basename = basename("~{smallsVcfFile}", ".vcf.gz")
		String modules = "tabix/1.9 bcftools/1.9 hg38/p12 hg38-dac-exclusion/1.0"
		String genome = "$HG38_ROOT/hg38_random.fa"
		String? difficultRegions = "--regions-file $HG38_DAC_EXCLUSION_ROOT/hg38-dac-exclusion.v2.bed"
		Float VAF = 0.01
		String QUALfilter = "FILTER~'haplotype' | FILTER~'clustered_events' | FILTER~'multiallelic' | FILTER~'slippage' |FILTER~'weak_evidence' | FILTER~'strand_bias' | FILTER~'position' | FILTER~'normal_artifact' |  FILTER~'map_qual' | FILTER~'germline' | FILTER~'fragment' | FILTER~'contamination' | FILTER~'base_qual'"
		Int jobMemory = 10
		Int threads = 1
		Int timeout = 2
	}

	parameter_meta {
		smallsVcfFile: "Vcf input file"
		smallsVcfIndex: "Vcf input index file"
		basename: "Base name"
		modules: "Required environment modules"
		sampleName: "Name of sample matching the tumor sample in .vcf"
		genome: "Path to loaded genome .fa"
		difficultRegions: "Path to .bed excluding difficult regions, string must include the flag --regions-file "
		VAF: "minimum variant allele frequency to retain variant"
		smallType: "type of variant to keep: snp or indel"
		jobMemory: "Memory allocated for this job (GB)"
		QUALfilter: "filter for filter calls to remove, eg. FILTER~'weak_evidence' | FILTER~'strand_bias' "
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail

		$BCFTOOLS_ROOT/bin/bcftools norm --multiallelics - --fasta-ref ~{genome} ~{difficultRegions} ~{smallsVcfFile} |\
		$BCFTOOLS_ROOT/bin/bcftools filter -i "TYPE='~{smallType}'" |\
		$BCFTOOLS_ROOT/bin/bcftools filter -e "~{QUALfilter}" |\
		$BCFTOOLS_ROOT/bin/bcftools filter -i "(FORMAT/AD[0:1])/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= ~{VAF}" >~{basename}.VAF.vcf

		bgzip ~{basename}.VAF.vcf
		tabix -p vcf ~{basename}.VAF.vcf.gz

		zcat ~{smallsVcfFile} | awk '$1 !~ "#" {print}'  | wc -l >~{sampleName}.filteringReport.txt
		zcat ~{basename}.VAF.vcf.gz | awk '$1 !~ "#" {print}'  | wc -l >>~{sampleName}.filteringReport.txt

	>>> 

	runtime {
		modules: "~{modules}"
		memory: "~{jobMemory} GB"
		cpu: "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File smallsVcfOutput = "~{basename}.VAF.vcf.gz"
		File smallsVcfIndexOutput = "~{basename}.VAF.vcf.gz.tbi"
		File smallsFilteringReport = "~{sampleName}.filteringReport.txt"
	}

	meta {
		output_meta: {
			smallsVcfOutput: "filtered .vcf",
			smallsVcfIndexOutput: "filtered .vcf.tbi indexed",
			smallsFilteringReport: "counts of variants pre and post filtering"
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
		File indelVcfIndexFiltered
		File snvVcfFiltered
		File snvVcfIndexFiltered
		File lohSegFile
		String oncotree
		String sampleName
		String modules = "hrdetect-scripts/1.3"
		String sigtoolrScript = "$HRDETECT_SCRIPTS_ROOT/bin/sigTools_runthrough.R"
		String SVrefSigs = "$SIGTOOLS_ROOT/lib/R/library/signature.tools.lib/data/RefSigv1_Rearr.tsv"
		String genomeVersion = "hg38"
		Int sigtoolsBootstrap = 2500
		Int indelCutoff = 10
		Int jobMemory = 30
		Int threads = 1
		Int timeout = 5
	}

	parameter_meta {
		structuralBedpeFiltered: "filtered structural variant .bedpe"
		SVrefSigs: "reference signatures for SVs"
		indelVcfFiltered: "filtered INDEL .vcf"
		snvVcfFiltered: "filtered SNV .vcf"
		indelVcfIndexFiltered: "filtered INDEL .vcf.tbi (indexed)"
		snvVcfIndexFiltered: "filtered SNV .vcf.tbi (indexed)"
		lohSegFile: "reformatted segmentation file"
		oncotree: "oncotree code of cancer"
		sigtoolrScript: ".R script containing sigtools"
		sampleName: "Name of sample matching the tumor sample in .vcf"		
		modules: "Required environment modules"
		genomeVersion: "version of genome, eg hg38"
		sigtoolsBootstrap: "Number of bootstraps for sigtools"
		indelCutoff: "minimum number of indels to run analysis"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail

		Rscript ~{sigtoolrScript} -s ~{sampleName} \
			-S ~{snvVcfFiltered} -I  ~{indelVcfFiltered} \
			-V ~{structuralBedpeFiltered} -L ~{lohSegFile} \
			-b ~{sigtoolsBootstrap} -g ~{genomeVersion} -i ~{indelCutoff} \
			-r ~{SVrefSigs}

	>>> 

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File JSONout = "~{sampleName}.signatures.json"
	}

	meta {
		output_meta: {
			JSONout : "JSON file of sigtools and CHORD signatures"
		}
	}
}