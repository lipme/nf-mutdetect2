nextflow.enable.dsl=2


process VARSCAN_mpileup2cns {

	tag "${sample}"
	publishDir("${params.outdir}/03_vcf/cns/");
	
	input:
	tuple val(sample),path(mpileup),val(genotype)

	output:
	path("${sample}.cns.vcf.gz*")
	
	script:
	
	
	min_coverage = params.varscan_mutant_min_coverage
	min_reads2 = params.varscan_mutant_min_reads2
	avg_qual = params.varscan_mutant_avg_qual
	var_freq = params.varscan_mutant_var_freq
	var_freq_for_hom = params.varscan_mutant_var_freq_for_hom
	pvalue = params.varscan_mutant_pvalue
	
	if (genotype == params.wt_genotype) {
		min_coverage = params.varscan_wt_min_coverage
		min_reads2 = params.varscan_wt_min_reads2
		avg_qual = params.varscan_wt_avg_qual
		var_freq = params.varscan_wt_var_freq
		var_freq_for_hom = params.varscan_wt_var_freq_for_hom
		pvalue = params.varscan_wt_pvalue
	}
	
	
	"""
	bgzip -cd ${mpileup} | varscan mpileup2cns --output-vcf 1 \
	 --min-coverage ${min_coverage} \
	 --min-reads2 ${min_reads2} \
	 --min-avg-qual ${avg_qual} \
	 --min-var-freq ${var_freq} \
	 --min-var-freq-for-hom ${var_freq_for_hom} \
	 --p-value ${pvalue}  | sed "s/Sample1/${sample}/" | bgzip -c > ${sample}.cns.vcf.gz
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
	
	min_coverage = params.varscan_mutant_min_coverage
	min_reads2 = params.varscan_mutant_min_reads2
	avg_qual = params.varscan_mutant_avg_qual
	var_freq = params.varscan_mutant_var_freq
	var_freq_for_hom = params.varscan_mutant_var_freq_for_hom
	pvalue = params.varscan_mutant_pvalue
	
	if (genotype == params.wt_genotype) {
		min_coverage = params.varscan_wt_min_coverage
		min_reads2 = params.varscan_wt_min_reads2
		avg_qual = params.varscan_wt_avg_qual
		var_freq = params.varscan_wt_var_freq
		var_freq_for_hom = params.varscan_wt_var_freq_for_hom
		pvalue = params.varscan_wt_pvalue
	}
	
	"""
	bgzip -cd ${mpileup} | varscan mpileup2snp --output-vcf 1 \
	 --min-coverage ${min_coverage} \
	 --min-reads2 ${min_reads2} \
	 --min-avg-qual ${avg_qual} \
	 --min-var-freq ${var_freq} \
	 --min-var-freq-for-hom ${var_freq_for_hom} \
	 --p-value ${pvalue}  | sed "s/Sample1/${sample}/" | bgzip -c > ${sample}.snp.vcf.gz
	"""
}

process VARSCAN_mpileup2indel {

	tag "${sample}"
	
	input:
	tuple val(sample),path(mpileup),val(genotype)
	
	output:
	tuple val(sample),path("${sample}.indel.vcf.gz"),val(genotype)

	script:
	
	min_coverage = params.varscan_mutant_min_coverage
	min_reads2 = params.varscan_mutant_min_reads2
	avg_qual = params.varscan_mutant_avg_qual
	var_freq = params.varscan_mutant_var_freq
	var_freq_for_hom = params.varscan_mutant_var_freq_for_hom
	pvalue = params.varscan_mutant_pvalue
	
	if (genotype == params.wt_genotype) {
		min_coverage = params.varscan_wt_min_coverage
		min_reads2 = params.varscan_wt_min_reads2
		avg_qual = params.varscan_wt_avg_qual
		var_freq = params.varscan_wt_var_freq
		var_freq_for_hom = params.varscan_wt_var_freq_for_hom
		pvalue = params.varscan_wt_pvalue
	}
	
	
	
	"""
	bgzip -cd ${mpileup} | varscan mpileup2indel --output-vcf 1 \
	 --min-coverage ${min_coverage} \
	 --min-reads2 ${min_reads2} \
	 --min-avg-qual ${avg_qual} \
	 --min-var-freq ${var_freq} \
	 --min-var-freq-for-hom ${var_freq_for_hom} \
	 --p-value ${pvalue}  | sed "s/Sample1/${sample}/" | bgzip -c > ${sample}.indel.vcf.gz
	"""
}
