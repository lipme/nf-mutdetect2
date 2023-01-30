nextflow.enable.dsl=2

process BWA_index {
	tag ("${refseq.baseName}")
	storeDir "store/"

	input:
	path(refseq)
	  
	output:
	path ("${refseq.baseName}-index")
	  
	script:
	"""
	mkdir ${refseq.baseName}-index
	bwa index -p ${refseq.baseName}-index/genome.fa -a bwtsw $refseq
	"""
}

process BWA_mem {

	tag("${sample}")
	label("process_high")
	
	input:
	path(index)
	tuple val(sample), path (read1), path (read2),val(type),val(group)
	
	output:
	tuple val(sample), path ("${sample}.sam.gz"),val(type),val(group)
	
	script:
	"""
	bwa mem -R "@RG\\tID:$sample\\tLB:$sample\\tSM:$group\\tPL:ILLUMINA" -t $task.cpus -M ${index}/genome.fa ${read1} ${read2} | gzip -c > ${sample}.sam.gz
	"""

}
