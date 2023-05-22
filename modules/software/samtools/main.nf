nextflow.enable.dsl=2

process SAMTOOLS_filter {

	tag "${sample}"
	label("process_medium")
	publishDir("${params.outdir}/02_bam/filtered/")

	input:
	path(refseq)
	tuple val(sample),path(samgz),val(genotype),val(group)
	
	
	output:
	tuple val(sample),path("${sample}.${params.aln_filter_mode}.bam"),val(genotype),val(group)
	
	script:
	if (params.aln_filter_mode == 'paired_rmdup')
	{
		"""
		bgzip -cd ${samgz} \
		| samtools view --threads=${task.cpus} -bT ${refseq} -f 0x02 - \
		| samtools fixmate -c -r - - \
		| samtools sort -m ${params.samtools_max_memory} - \
		| samtools rmdup - - \
		| samtools view --threads=${task.cpus} -b -q 1 -F 4 -F 256 > ${sample}.${params.aln_filter_mode}.bam
		"""
	}
	else if (params.aln_filter_mode == 'rmdup')
	{
		"""
		bgzip -cd ${samgz} \
		| samtools view --threads=${task.cpus} -bT ${refseq} - \
		| samtools fixmate -c -r - - \
		| samtools sort -m ${params.samtools_max_memory} - \
		| samtools rmdup - - \
		| samtools view --threads=${task.cpus} -b -q 1 -F 4 -F 256 > ${sample}.${params.aln_filter_mode}.bam
		"""
	}
	else if (params.aln_filter_mode == 'paired')
	{
		"""
		bgzip -cd ${samgz} \
		| samtools view --threads=${task.cpus} -bT ${refseq} -f 0x02 - \
		| samtools fixmate -c -r - - \
		| samtools sort -m ${params.samtools_max_memory} - \
		| samtools view --threads=${task.cpus} -b -q 1 -F 4 -F 256 > ${sample}.${params.aln_filter_mode}.bam
		"""
	}
	
}


process SAMTOOLS_mpileup {
	tag "${sample}"
	label("process_medium")
	
	input:
	path(refseq)
	tuple val(sample),path(bam),val(genotype),val(group)
	
	
	output:
	tuple val(group),path("${group}.mpileup.gz"),val(gt)


	
	script:
	gt = genotype.first()
	"""
	samtools merge --threads=${task.cpus} - ${bam} | samtools mpileup -a ${params.mpileup_parameters} -f ${refseq} - | bgzip -c > ${group}.mpileup.gz
	"""
}

process SAMTOOLS_mpileupRecall {

	tag "${sample}"
	label("process_medium")
	
	input:
	path(positions)
	path(refseq)
	tuple val(sample),path(bam),val(genotype),val(group)
	
	
	output:
	tuple val(group),path("${group}.mpileup.cns.gz"),val(gt)


	
	script:
	gt = genotype.first()
	"""
	samtools merge --threads=${task.cpus} - ${bam} | samtools mpileup ${params.mpileuprecall_parameters} -f ${refseq} -l ${positions} - | bgzip -c > ${group}.mpileup.cns.gz
	"""

}
	
process SAMTOOLS_countaln {

	label("process_medium")
	
	input:
	tuple val(sample), path(bam), val(genotype),val(group)
	
	output:
	path("${sample}.aln.filtered")
	
	script:
	"""
	ALN=`samtools flagstat --threads=${task.cpus} ${bam} | head -1 | cut -f1 -d ' '`
	COVERAGE=`samtools depth -aa ${bam}| awk '{c++;s+=\$3}END{print s/c}'`

	echo -e "${sample}\t${group}\t${genotype}\t\$ALN\t\$COVERAGE" > ${sample}.aln.filtered
	"""
}
	
