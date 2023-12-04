version 1.0

struct GenomeResources {
	String filterSMALLsModules
	String genome
	String difficultRegions
}

workflow hrDetect {
	input {
		String outputFileNamePrefix
		File structuralVcfFile
		File smallsVcfFile
		File smallsVcfIndex
		File segFile
		String reference
	}

	parameter_meta {
		structuralVcfFile: "Input VCF file of structural variants (eg. from delly)"
		smallsVcfFile: "Input VCF file of SNV and indels (small mutations) (eg. from mutect2)"
		smallsVcfIndex: "Index file for smallsVcfFile"
		segFile: "File for segmentations, used to estimate number of segments in Loss of heterozygosity (LOH) (eg. from sequenza)"
		outputFileNamePrefix: "Name of sample matching the tumor sample in .vcf"
		reference: "Reference genome version"
	}

	Map[String,GenomeResources] resources = {
		"hg38": {
			"filterSMALLsModules": "tabix/1.9 bcftools/1.9 hg38/p12 hg38-dac-exclusion/1.0",
			"genome": "$HG38_ROOT/hg38_random.fa",
			"difficultRegions": "--regions-file $HG38_DAC_EXCLUSION_ROOT/hg38-dac-exclusion.v2.bed"
		}
	}

	call filterStructural {
		input: 
			outputFileNamePrefix = outputFileNamePrefix,
			structuralVcfFile = structuralVcfFile
	}

	call filterSMALLs as filterINDELs {
		input: 
			outputFileNamePrefix = outputFileNamePrefix,
			smallsVcfFile = smallsVcfFile,
			smallsVcfIndex = smallsVcfIndex,
			smallType = "indel",
			modules = resources[reference].filterSMALLsModules,
			genome = resources[reference].genome,
			difficultRegions = resources[reference].difficultRegions
	}

	call filterSMALLs as filterSNVs {
		input: 
			outputFileNamePrefix = outputFileNamePrefix,
			smallsVcfFile = smallsVcfFile,
			smallsVcfIndex = smallsVcfIndex,
			smallType = "snp",
			modules = resources[reference].filterSMALLsModules,
			genome = resources[reference].genome,
			difficultRegions = resources[reference].difficultRegions
	}

	call hrdResults {
		input:
			outputFileNamePrefix = outputFileNamePrefix,
			SV_vcf_location = filterStructural.SV_vcf_location,
			indelVcfFiltered = filterINDELs.smallsVcfOutput,
			indelVcfIndexFiltered = filterINDELs.smallsVcfIndexOutput,
			snvVcfFiltered = filterSNVs.smallsVcfOutput,
			snvVcfIndexFiltered = filterSNVs.smallsVcfIndexOutput,
			lohSegFile = segFile,
			genomeVersion = reference
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
			hrd_signatures : "JSON file of hrdetect signatures",
			SV_exposures: "JSON of structural variant signatures",
			SV_catalog: "JSON cataloguing structural variants",
			ID_catalog: "JSON cataloguing indels"
		}
	}
	output {
		File hrd_signatures = hrdResults.hrd_signatures
		File SBS_exposures = hrdResults.SBS_exposures
		File SV_exposures = hrdResults.SV_exposures
		File ID_catalog = hrdResults.ID_catalog
	}
}

task filterStructural {
	input {
		String outputFileNamePrefix
		File structuralVcfFile 
		String modules = "bcftools/1.9"
		String structuralQUALfilter = "PASS"
		Int jobMemory = 5
		Int threads = 1
		Int timeout = 1
	}

	parameter_meta {
		structuralVcfFile: "Vcf input file"
		modules: "Required environment modules"
		outputFileNamePrefix: "Name of sample matching the tumor sample in .vcf"
		structuralQUALfilter: "filter for filter calls to keep, eg. PASS"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail


		$BCFTOOLS_ROOT/bin/bcftools view -f '~{structuralQUALfilter}' ~{structuralVcfFile} >> ~{outputFileNamePrefix}.structural.PASS.vcf

		awk '$1 !~ "#" {print}' ~{structuralVcfFile} | wc -l >~{outputFileNamePrefix}.structural.filteringReport.txt
		awk '$1 !~ "#" {print}' ~{outputFileNamePrefix}.structural.PASS.vcf | wc -l >>~{outputFileNamePrefix}.structural.filteringReport.txt

	>>>

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File SV_vcf_location = "~{outputFileNamePrefix}.structural.PASS.vcf"
		File structuralFilteringReport = "~{outputFileNamePrefix}.structural.filteringReport.txt"
	}

	meta {
		output_meta: {
			SV_vcf_location: "filtered structural .vcf",
			structuralFilteringReport: "counts of variants pre and post filtering"
		}
	}
}


task filterSMALLs {
	input {
		String outputFileNamePrefix
		File smallsVcfIndex
		String smallType
		File smallsVcfFile
		String modules
		String genome
		String? difficultRegions
		Float VAF = 0.01
		String QUALfilter = "FILTER~'haplotype' | FILTER~'clustered_events' | FILTER~'multiallelic' | FILTER~'slippage' |FILTER~'weak_evidence' | FILTER~'strand_bias' | FILTER~'position' | FILTER~'normal_artifact' |  FILTER~'map_qual' | FILTER~'germline' | FILTER~'fragment' | FILTER~'contamination' | FILTER~'base_qual'"
		Int jobMemory = 10
		Int threads = 1
		Int timeout = 2
	}

	parameter_meta {
		smallsVcfFile: "Vcf input file"
		modules: "Required environment modules"
		outputFileNamePrefix: "Name of sample matching the tumor sample in .vcf"
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

		$BCFTOOLS_ROOT/bin/bcftools norm --multiallelics - --fasta-ref ~{genome} ~{difficultRegions} ~{smallsVcfFile} | \
		$BCFTOOLS_ROOT/bin/bcftools filter -i "TYPE='~{smallType}'" | \
		$BCFTOOLS_ROOT/bin/bcftools filter -e "~{QUALfilter}" | \
		$BCFTOOLS_ROOT/bin/bcftools filter -i "(FORMAT/AD[0:1])/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= ~{VAF}" >~{outputFileNamePrefix}.~{smallType}.VAF.vcf

		bgzip ~{outputFileNamePrefix}.~{smallType}.VAF.vcf
		tabix -p vcf ~{outputFileNamePrefix}.~{smallType}.VAF.vcf.gz

		zcat ~{smallsVcfFile} | awk '$1 !~ "#" {print}'  | wc -l >~{outputFileNamePrefix}.~{smallType}.filteringReport.txt
		zcat ~{outputFileNamePrefix}.~{smallType}.VAF.vcf.gz | awk '$1 !~ "#" {print}'  | wc -l >>~{outputFileNamePrefix}.~{smallType}.filteringReport.txt

	>>> 

	runtime {
		modules: "~{modules}"
		memory: "~{jobMemory} GB"
		cpu: "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File smallsVcfOutput = "~{outputFileNamePrefix}.~{smallType}.VAF.vcf.gz"
		File smallsVcfIndexOutput = "~{outputFileNamePrefix}.~{smallType}.VAF.vcf.gz.tbi"
		File smallsFilteringReport = "~{outputFileNamePrefix}.~{smallType}.filteringReport.txt"
	}

	meta {
		output_meta: {
			smallsVcfOutput: "filtered .vcf",
			smallsVcfIndexOutput: "filtered .vcf.tbi indexed",
			smallsFilteringReport: "counts of variants pre and post filtering"
		}
	}
}

task hrdResults {
	input {
		String outputFileNamePrefix
		File SV_vcf_location
		File indelVcfFiltered
		File indelVcfIndexFiltered
		File snvVcfFiltered
		File snvVcfIndexFiltered
		File lohSegFile
		String modules = "sigtools/2.4.1 sigtools-data/1.0 sigtools-rscript/1.5.8"
		String sigtoolrScript = "$SIGTOOLS_RSCRIPT_ROOT/scripts/sigTools_runthrough.R"
		String SVrefSigs = "$SIGTOOLS_DATA_ROOT/RefSigv0_Rearr.tsv"
		String SNVrefSigs = "$SIGTOOLS_DATA_ROOT/COSMIC_v1_SBS_GRCh38.txt"
		String genomeVersion
		Int sigtoolsBootstrap = 200
		Int indelCutoff = 50
		Int jobMemory = 50
		Int threads = 1
		Int timeout = 15
	}

	parameter_meta {
		SV_vcf_location: "structural variant vcf"
		SVrefSigs: "reference signatures for SVs"
		SNVrefSigs: "reference signatures for SNVs"
		indelVcfFiltered: "filtered INDEL .vcf"
		snvVcfFiltered: "filtered SNV .vcf"
		indelVcfIndexFiltered: "filtered INDEL .vcf.tbi (indexed)"
		snvVcfIndexFiltered: "filtered SNV .vcf.tbi (indexed)"
		lohSegFile: "reformatted segmentation file"
		sigtoolrScript: ".R script containing sigtools"
		outputFileNamePrefix: "Name of sample matching the tumor sample in .vcf"		
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

	>>> 

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File hrd_signatures = "~{outputFileNamePrefix}.signatures.json"
		File SBS_exposures = "~{outputFileNamePrefix}.exposures.SBS.json"
		File SV_exposures = "~{outputFileNamePrefix}.exposures.SV.json"
		File ID_catalog = "~{outputFileNamePrefix}.catalog.ID.json"
	}

	meta {
		output_meta: {
			hrd_signatures : "JSON file of hrdetect signatures",
			SBS_exposures: "JSON of single basepair substitution signatures",
			SV_exposures: "JSON of structural variant signatures",
			ID_catalog: "JSON cataloguing indels"
		}
	}
}
