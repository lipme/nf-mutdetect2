nextflow.enable.dsl=2

process VCFTOOLS_excludePostions {

	tag("$sample")
	
	
	input:
	path(positions)
	tuple val(sample), path(vcf), path(tbi)
	
	output:
	tuple val(sample), path("${sample}.noWT.vcf.gz"), path("${sample}.noWT.vcf.gz.tbi")
	path("${sample}.WTcommon.count")
	
	script:
	"""
	vcftools --gzvcf $vcf --exclude-positions $positions --recode --recode-INFO-all  --stdout| bgzip -c > ${sample}.noWT.vcf.gz
	tabix ${sample}.noWT.vcf.gz
	BEFORE=`bgzip -cd $vcf  | (grep -v '^#' || true) | wc -l`
	AFTER=`bgzip -cd ${sample}.noWT.vcf.gz  |  (grep -v '^#' || true) | wc -l`
	DIFF=\$((\$BEFORE-\$AFTER))
	echo -e "${sample}\t\$DIFF" > ${sample}.WTcommon.count
	"""
	
}

process VCFTOOLS_convertToTsv {

	tag("$sample")
	publishDir("${params.outdir}/05_tsv")
	
	input:
	tuple val(sample), path(vcf), path(tbi)
	
	output:
	path("${vcf.baseName}.tsv")
	
	script:
	"""
	vcftools --gzvcf ${vcf}  --get-INFO ANN --stdout  > ANN.tsv
	vcftools --gzvcf ${vcf}  --extract-FORMAT-info GT  --stdout | cut -f3- > GT.tsv
	paste ANN.tsv GT.tsv > ${vcf.baseName}.tsv
	"""

}

