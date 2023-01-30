nextflow.enable.dsl=2

process TRIMMOMATIC_trimmomatic {

	label "process_high"
	
	input:
	tuple val(sample), path(read1), path(read2),val(genotype),val(readlen),val(group)
	
	output:
	tuple val(sample), path ("${sample}.trimmed_1P.fastq.gz"), path ("${sample}.trimmed_2P.fastq.gz"),val(genotype),val(group)
	tuple val(sample),val(genotype),path("${sample}.trimmed.summary"),val(group)
	
	script:
	"""
	trimmomatic PE -threads ${task.cpus} -summary ${sample}.trimmed.summary ${params.trimmo_params} ${read1} ${read2} \
	-baseout ${sample}.trimmed.fastq.gz \
	LEADING:${params.trimmo_leading_qual_min} \
	TRAILING:${params.trimmo_trailing_qual_min} \
	SLIDINGWINDOW:${params.trimmo_sliding_window_size}:${params.trimmo_sliding_window_qual} \
	MINLEN:${params.trimmo_min_length}
	"""

}
