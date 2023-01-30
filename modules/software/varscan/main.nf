nextflow.enable.dsl=2


process VARSCAN_mpileup2cns {

	tag "${sample}"
	publishDir("${params.outdir}/03_vcf/cns/");
	
	input:
	tuple val(sample),path(mpileup),val(genotype)

	output:
	path("${sample}.cns.vcf.gz*")
	
	script:
	"""
	bgzip -cd ${mpileup} | varscan mpileup2cns --output-vcf 1 \
	 --min-coverage ${params.varscan_min_coverage} \
	 --min-reads2 ${params.varscan_min_reads2} \
	 --min-avg-qual ${params.varscan_avg_qual} \
	 --min-var-freq ${params.varscan_var_freq} \
	 --min-var-freq-for-hom ${params.varscan_var_freq_for_hom} \
	 --p-value ${params.varscan_pvalue}  | sed "s/Sample1/${sample}/" | bgzip -c > ${sample}.cns.vcf.gz
	tabix ${sample}.cns.vcf.gz
	"""
	
}

process VARSCAN_mpileup2snp {

	tag "${sample}"
	
	input:
	tuple val(sample),path(mpileup),val(genotype)
	
	output:
	tuple val(sample),path("${sample}.snp.vcf.gz"),val(genotype)

	script:
	"""
	bgzip -cd ${mpileup} | varscan mpileup2snp --output-vcf 1 \
	 --min-coverage ${params.varscan_min_coverage} \
	 --min-reads2 ${params.varscan_min_reads2} \
	 --min-avg-qual ${params.varscan_avg_qual} \
	 --min-var-freq ${params.varscan_var_freq} \
	 --min-var-freq-for-hom ${params.varscan_var_freq_for_hom} \
	 --p-value ${params.varscan_pvalue}  | sed "s/Sample1/${sample}/" | bgzip -c > ${sample}.snp.vcf.gz
	"""
}

process VARSCAN_mpileup2indel {

	tag "${sample}"
	
	input:
	tuple val(sample),path(mpileup),val(genotype)
	
	output:
	tuple val(sample),path("${sample}.indel.vcf.gz"),val(genotype)

	script:
	"""
	bgzip -cd ${mpileup} | varscan mpileup2indel --output-vcf 1 \
	 --min-coverage ${params.varscan_min_coverage} \
	 --min-reads2 ${params.varscan_min_reads2} \
	 --min-avg-qual ${params.varscan_avg_qual} \
	 --min-var-freq ${params.varscan_var_freq} \
	 --min-var-freq-for-hom ${params.varscan_var_freq_for_hom} \
	 --p-value ${params.varscan_pvalue}  | sed "s/Sample1/${sample}/" | bgzip -c > ${sample}.indel.vcf.gz
	"""
}
