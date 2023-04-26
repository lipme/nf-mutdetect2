nextflow.enable.dsl=2

process BCFTOOLS_snpeff_exclude {

	errorStrategy { task.exitStatus in 0..1 ? 'ignore' : 'terminate' }

	input:
	tuple val(sample), path(snpeff_vcf), path(tbi)
	
	output:
	tuple val(sample), path("${sample}.noWT.ann_filter.vcf.gz"),path("${sample}.noWT.ann_filter.vcf.gz.tbi")
	path("${sample}.noWT.ann_filter.count")
	
	script:
	"""
	BEFORE=`bgzip -cd ${snpeff_vcf} | (grep -v '^#' || true)  | wc -l`
	bcftools  filter -i 'ANN[*] !~ "${params.snpeff_exclude_ann}"'  ${snpeff_vcf} | bgzip -c > ${sample}.noWT.ann_filter.vcf.gz
	tabix ${sample}.noWT.ann_filter.vcf.gz
	AFTER=`bgzip -cd ${sample}.noWT.ann_filter.vcf.gz | (grep -v '^#' || true) | wc -l`
	echo -e "${sample}\t\$BEFORE\t\$AFTER" > ${sample}.noWT.ann_filter.count
	"""

}

process BCFTOOLS_getDeletions {
	
	publishDir("${params.outdir}/03_vcf/del")
	
	input:
	tuple val(sample), path(indel_vcf),path(tbi)
	
	output:
	tuple val(sample), path("${sample}.del.vcf.gz"),path("${sample}.del.vcf.gz.tbi")
	path("${sample}.del.count")
	
	script:
	"""
	bcftools view --types indels ${indel_vcf} | bcftools norm -m - |   bcftools filter --include 'strlen(REF)>strlen(ALT)' |   bcftools view -H | bgzip -c > ${sample}.del.vcf.gz
	tabix ${sample}.del.vcf.gz
	COUNT=`bgzip -cd  ${sample}.del.vcf.gz | (grep -c -v '^#' || true)`
	echo -e "${sample}\t\$COUNT" >  ${sample}.del.count
	"""
	
}

process BCFTOOLS_merge {

	label("process_medium")
 
 	input:
	path(cns_vcf)
	val (prefix)
	val (suffix)
	
	output:
	tuple val(prefix), path("${prefix}.${suffix}.raw.vcf.gz"), path("${prefix}.${suffix}.raw.vcf.gz.tbi")
	
	script:
	"""
	find . -name '*.cns.vcf.gz' > input.files
	bcftools merge --threads ${task.cpus} --info-rules ADP:avg,HOM:sum,HET:sum,NC:sum,WT:sum --file-list input.files | bgzip -c > ${prefix}.${suffix}.raw.vcf.gz
	tabix ${prefix}.${suffix}.raw.vcf.gz
	"""
	
	
}

process BCFTOOLS_norm {

	label("process_medium")
	
	input:
	tuple val(prefix), path(vcf), path(tbi)
	path(reference)
	val (prefix)
	
	output:
	tuple val(prefix), path("${prefix}.norm.vcf.gz")
	
	script:
	"""
	bcftools norm --threads ${task.cpus} -m-any --check-ref x --fasta-ref ${reference} --output ${prefix}.norm.vcf.gz --output-type z ${vcf}  2> norm.log
	"""
	
}

process BCFTOOLS_stats {

	label("process_medium")
	
	publishDir("${params.outdir}/06_summaries/")
	
	input:
	tuple val(sample), path(vcf),path(tbi)
	
	output:
	path("${vcf.baseName}.stats.tsv")
	path("${vcf.baseName}.stats.perSample.tsv")
	
	script:
	"""
	bcftools stats --threads ${task.cpus} --samples -  ${vcf} > ${vcf.baseName}.stats.tsv
	grep -B2 ^PS ${vcf.baseName}.stats.tsv | cut -f3- > ${vcf.baseName}.stats.perSample.tsv
	"""
}
