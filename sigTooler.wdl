version 1.0

workflow sigTooler {
	input {
		File structuralVcfFile
		File smallsVcfFile
		File smallsVcfIndex
		File segFile
		String indelVAF
		String snvVAF
		String tissue
		String rScript
		String sampleName
		Boolean plotIt
	}

	parameter_meta {
		structuralVcfFile: "Input VCF file of structural variants (eg. from delly)"
		smallsVcfFile: "Input VCF file of SNV and indels (small mutations) (eg. from mutect2)"
		smallsVcfIndex: "Index of input VCF file of SNV and indels"
		segFile: "File for segmentations, used to estimate number of segments in Loss of heterozygosity (LOH) (eg. from sequenza)"
		indelVAF: "Variant Allele Frequency for filtering of indel mutations"
		snvVAF: "Variant Allele Frequency for filtering of SNVs"
		tissue: "Cancerous-tissue of origin"
		rScript: "Temporary variable to call the .R script containing sigtools, will be modulated. default: ~/sigtools_workflow/sigTools_runthrough.R"
		sampleName: "Name of sample matching the tumor sample in .vcf"
	}

	call filterStructural {
		input: 
			structuralVcfFile = structuralVcfFile,
			sampleName = sampleName
	}

	call filterINDELs {
		input: 
			smallsVcfFile = smallsVcfFile,
			smallsVcfIndex = smallsVcfIndex,
			indelVAF = indelVAF,
			sampleName = sampleName
	}

	call filterSNVs {
		input: 
			smallsVcfFile = smallsVcfFile,
			smallsVcfIndex = smallsVcfIndex,
			snvVAF = snvVAF,
			sampleName = sampleName
	}

	call convertSegFile {
		input: 
			segFile = segFile
	}

	call hrdResults {
		input:
			structuralBedpeFiltered = filterStructural.structuralbedpe,
			indelVcfFiltered = filterINDELs.indelVcfOutput,
			indelVcfIndexFiltered = filterINDELs.indelVcfIndexOutput,
			snvVcfFiltered = filterSNVs.snvVcfOutput,
			snvVcfIndexFiltered = filterSNVs.snvVcfIndexOutput,
			lohSegFile = convertSegFile.segmentsOutput,
			tissue = tissue,
			rScript = rScript,
			sampleName = sampleName
	}

	if(plotIt == true){
		call plotResults {
			input:
				sigTools_hrd_input = hrdResults.sigTools_hrd_Output ,
				sigTools_sigs_input = hrdResults.sigTools_sigs_Output,
				rScript = rScript,
				sampleName = sampleName
		}
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
				name: "bis-rlibs/0.1",
				url: "https://ggplot2.tidyverse.org/"
			}
		]
		output_meta: {
			sigTools_hrd_Output : "point estimate and bootstraped confidence intervals for HRD from sigtools",
			sigTools_model_Output : "parameters raw values and weights for estimation of HRD from sigtools",
			sigTools_sigs_Output : "signature breakdown from sigtools" ,
			sigTools_sigs_plot_Output: "plot of signature breakdown from sigtools",
			sigTools_hrd_plot_Output: "plot of point estimate and bootstraped confidence intervals for HRD from sigtools",
			indelFilteringReport: "counts of INDELs pre and post filtering",
			snvFilteringReport: "counts of SNVs pre and post filtering",
			structuralFilteringReport: "counts of structural variants pre and post filtering"
		}
	}
	output {
		File indelFilteringReport = "~{sampleName}.INDEL.filteringReport.txt"
		File snvFilteringReport = "~{sampleName}.SNP.filteringReport.txt"
		File structuralFilteringReport = "~{sampleName}.structural.filteringReport.txt"
		File sigTools_hrd_Output = "~{sampleName}.sigtools.hrd.txt"
		File sigTools_model_Output = "~{sampleName}.sigtools.model.txt"
		File sigTools_sigs_Output = "~{sampleName}.sigtools.sigs.txt"
		File? sigTools_sigs_plot_Output = "~{sampleName}.sigtools.sigs.png"
		File? sigTools_hrd_plot_Output = "~{sampleName}.sigtools.hrd.png"
	}
}

task filterStructural {
	input {
		File structuralVcfFile 
		String basename = basename("~{structuralVcfFile}", ".vcf.gz")
		String modules = "bcftools/1.9"
		String sampleName
		Int jobMemory = 5
		Int threads = 1
		Int timeout = 1
	}

	parameter_meta {
		structuralVcfFile: "Vcf input file"
		basename: "Base name"
		modules: "Required environment modules"
		sampleName: "Name of sample matching the tumor sample in .vcf"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail

		echo  -e "chrom1\tstart1\tend1\tchrom2\tstart2\tend2\tsample\tsvclass"  >~{basename}.bedpe

		$BCFTOOLS_ROOT/bin/bcftools view -f 'PASS' ~{structuralVcfFile} |\
		$BCFTOOLS_ROOT/bin/bcftools filter -e 'INFO/SVTYPE = "BND"' |\
		$BCFTOOLS_ROOT/bin/bcftools query -f "%CHROM\t%POS\t%INFO/END\t%FILTER\t%INFO/SVTYPE\t%INFO/CIPOS\t%INFO/CIEND\n" |\
		awk -v sampleName=~{sampleName} 'split($6,a,",") split($7,b,",") {print $1"\t"$2+a[1]-1"\t"$2+a[2]"\t"$1"\t"$3+b[1]-1"\t"$2+b[2]"\t"sampleName"\t"$5} '  >>~{basename}.bedpe

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

task filterINDELs {
	input {
		File smallsVcfFile
		File smallsVcfIndex
		String basename = basename("~{smallsVcfFile}", ".vcf.gz")
		String modules = "gatk/4.2.0.0 tabix/1.9 bcftools/1.9 hg38/p12 grch38-alldifficultregions/3.0"
		String sampleName
		String genome = "$HG38_ROOT/hg38_random.fa"
		String? difficultRegions
		String indelVAF
		Int jobMemory = 10
		Int threads = 1
		Int timeout = 2
	}

	parameter_meta {
		smallsVcfFile: "Vcf input file"
		basename: "Base name"
		modules: "Required environment modules"
		sampleName: "Name of sample matching the tumor sample in .vcf"
		genome: "Path to loaded genome"
		difficultRegions: "Path to .bed of difficult regions to align to, string must include the --exclude-intervals flag, eg: --exclude-intervals $GRCH38_ALLDIFFICULTREGIONS_ROOT/GRCh38_alldifficultregions.bed"
		indelVAF: "VAF for indels"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail

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

	>>> 

	runtime {
		modules: "~{modules}"
		memory: "~{jobMemory} GB"
		cpu: "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File indelVcfOutput = "~{basename}.INDEL.VAF.vcf.gz"
		File indelVcfIndexOutput = "~{basename}.INDEL.VAF.vcf.gz.tbi"
		File indelFilteringReport = "~{sampleName}.INDEL.filteringReport.txt"
	}

	meta {
		output_meta: {
			indelVcfOutput: "filtered INDEL .vcf",
			indelVcfIndexOutput: "filtered INDEL .vcf.tbi indexed",
			indelFilteringReport: "counts of variants pre and post filtering"
		}
	}
}

task filterSNVs {
	input {
		File smallsVcfFile
		File smallsVcfIndex
		String basename = basename("~{smallsVcfFile}", ".vcf.gz")
		String modules = "gatk/4.2.0.0 tabix/1.9 bcftools/1.9 grch38-alldifficultregions/3.0 hg38/p12"
		String sampleName
		String genome = "$HG38_ROOT/hg38_random.fa"
		String? difficultRegions 
		String snvVAF
		Int jobMemory = 10
		Int threads = 1
		Int timeout = 2
	}

	parameter_meta {
		smallsVcfFile: "Vcf input file"
		snvVAF: "VAF for SNV filtering"
		basename: "Base name"
		modules: "Required environment modules"
		sampleName: "Name of sample matching the tumor sample in .vcf"
		genome: "Path to loaded genome"
		difficultRegions: "Path to .bed of difficult regions to align to, string must include the --exclude-intervals flag, eg: --exclude-intervals $GRCH38_ALLDIFFICULTREGIONS_ROOT/GRCh38_alldifficultregions.bed"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail

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

	>>> 

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File snvVcfOutput = "~{basename}.SNP.VAF.vcf.gz"
		File snvVcfIndexOutput = "~{basename}.SNP.VAF.vcf.gz.tbi"
		File snvFilteringReport = "~{sampleName}.SNP.filteringReport.txt"
	}

	meta {
		output_meta: {
			snvVcfOutput: "filtered SNV .vcf",
			snvVcfIndexOutput: "filtered SNV .vcf.tbi (indexed)",
			snvFilteringReport: "counts of variants pre and post filtering"
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
		String tissue
		String rScript
		String sampleName
		String modules = "sigtools/0.0.0.9000"
		Int sigtoolsBootstrap = 2500
		Int jobMemory = 20
		Int threads = 1
		Int timeout = 2
	}

	parameter_meta {
		structuralBedpeFiltered: "filtered structural variant .bedpe"
		indelVcfFiltered: "filtered INDEL .vcf"
		snvVcfFiltered: "filtered SNV .vcf"
		indelVcfIndexFiltered: "filtered INDEL .vcf.tbi (indexed)"
		snvVcfIndexFiltered: "filtered SNV .vcf.tbi (indexed)"
		lohSegFile: "reformatted segmentation file"
		tissue: "Cancerous-tissue of origin"
		rScript: "Temporary variable to call the .R script containing sigtools, will be modulated. default: ~/sigtools_workflow/sigTools_runthrough.R"
		sampleName: "Name of sample matching the tumor sample in .vcf"		
		modules: "Required environment modules"
		sigtoolsBootstrap: "Number of bootstraps for sigtools"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail

		Rscript --vanilla ~{rScript}_runthrough.R ~{sampleName} ~{tissue} ~{snvVcfFiltered} ~{indelVcfFiltered} ~{structuralBedpeFiltered} ~{lohSegFile} ~{sigtoolsBootstrap}

	>>> 

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File sigTools_hrd_Output = "~{sampleName}.sigtools.hrd.txt"
		File sigTools_model_Output = "~{sampleName}.sigtools.model.txt"
		File sigTools_sigs_Output = "~{sampleName}.sigtools.sigs.txt"
	}

	meta {
		output_meta: {
			sigTools_hrd_Output : "point estimate and bootstraped confidence intervals for HRD from sigtools",
			sigTools_model_Output : "parameters raw values and weights for estimation of HRD from sigtools",
			sigTools_sigs_Output : "signature breakdown from sigtools" 
		}
	}
}

task plotResults {
	input {
		File sigTools_hrd_input 
		File sigTools_sigs_input 
		String rScript
		String sampleName
		String modules = "bis-rlibs/0.1"
		Int jobMemory = 20
		Int threads = 1
		Int timeout = 2
	}

	parameter_meta {
		rScript: "Temporary variable to call the .R script containing sigtools, will be modulated. default: ~/sigtools_workflow/sigTools_runthrough.R"
		sampleName: "Name of sample matching the tumor sample in .vcf"		
		modules: "Required environment modules"
		jobMemory: "Memory allocated for this job (GB)"
		threads: "Requested CPU threads"
		timeout: "Hours before task timeout"
	}

	command <<<
		set -euo pipefail

		Rscript --vanilla ~{rScript}_plotter.R ~{sampleName}

	>>> 

	runtime {
		modules: "~{modules}"
		memory:  "~{jobMemory} GB"
		cpu:     "~{threads}"
		timeout: "~{timeout}"
	}

	output {
		File sigTools_sigs_plot_Output = "~{sampleName}.sigtools.sigs.png"
		File sigTools_hrd_plot_Output = "~{sampleName}.sigtools.hrd.png"
	}

	meta {
		output_meta: {
			sigTools_sigs_plot_Output: "plot of signature breakdown from sigtools",
			sigTools_hrd_plot_Output: "plot of point estimate and bootstraped confidence intervals for HRD from sigtools"
		}
	}
}
