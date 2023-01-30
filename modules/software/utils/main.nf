nextflow.enable.dsl=2


process UTILS_countSnps {

	input:
	tuple val(sample),path(vcffile),val(genotype)
	
	output:
	path("${sample}.snps.count")
	
	script:
	"""
	COUNT=`gzip -cd ${vcffile} | (grep -v ^# -c || true)` 
	echo -e "${sample}\t${genotype}\t\$COUNT" > ${sample}.snps.count
	"""
	
}

process UTILS_getSnpPositions {
	
	publishDir("${params.outdir}/04_positions/");
	input:
	path(vcffiles)
	val(genotype)
	
	output:
	path("${genotype}.snpPositions.tsv")
	
	script:
	"""
	zcat $vcffiles | (grep -v ^# || true) | cut -f1,2 | sort -u | sort -k1,1 -k2n,2 > ${genotype}.snpPositions.tsv
	"""
}

process UTILS_genomeSize {

	input:
	path(refseq)
	
	output:
	stdout
	
	script:
	"""
	grep -v '^>' ${refseq} | tr -d "\\n" | wc -c
	"""
}

process UTILS_libSize {
	tag "$sample"

	errorStrategy { task.exitStatus in [0,141] ? 'ignore' : 'terminate' }
	
	input:
	
	tuple val(sample), path(read1), path(read2),val(genotype), val(readlen), val(group),val(genomesize)

	output:
	path("${sample}.libsize")

	script:

	"""
	zcat $read1 | echo \$(((`wc -l`/4)*2)) > ${sample}.readcount.tmp
	echo "${sample}\t${group}\t${genotype}" > ${sample}.name
	READCOUNT=`cat ${sample}.readcount.tmp`
	BASES=\$((${readlen} * \$READCOUNT))
	echo \$BASES > ${sample}.bases.tmp
	COVERAGE=\$((\$BASES / ${genomesize}))
	echo \$COVERAGE > ${sample}.coverage.tmp
	paste ${sample}.name ${sample}.readcount.tmp ${sample}.bases.tmp ${sample}.coverage.tmp> ${sample}.libsize
	"""
}

process UTILS_rawaln {

	input:
	tuple val(sample), path(samgz), val(genotype),val(group)
	
	output:
	path("${sample}.aln.raw")
	"""
	ALN=`gzip -cd    ${samgz} | (grep -c -v ^@ || true)`
	echo -e "${sample}\t${group}\t${genotype}\t\$ALN" > ${sample}.aln.raw
	"""
}

process UTILS_filterTransitions {

	publishDir("${params.outdir}/03_vcf/transitions");

	errorStrategy { task.exitStatus in 0..1 ? 'ignore' : 'terminate' }
	
	input:
	tuple val(sample), path(snpeff_fil_vcf),path (tbi)
	
	output:
	tuple val(sample), path("${sample}.snpeff.ann_filter.transitions.vcf.gz"),path("${sample}.snpeff.ann_filter.transitions.vcf.gz.tbi")
	path("${sample}.snpeff.ann_filter.transitions.count")
	
	
	script:
	"""
	bgzip -cd  ${snpeff_fil_vcf} | awk -F '\t' '((\$4 == "A" && \$5 == "G") || (\$4 == "G" && \$5 == "A") || (\$4 == "C" && \$5 == "T") || (\$4 == "T" && \$5 == "C") || \$0 ~ /^#/ )' | bgzip -c > ${sample}.snpeff.ann_filter.transitions.vcf.gz
	tabix ${sample}.snpeff.ann_filter.transitions.vcf.gz
	AFTER=`bgzip -cd ${sample}.snpeff.ann_filter.transitions.vcf.gz | (grep -v '^#' || true) | wc -l`
	echo -e "${sample}\t\$AFTER" > ${sample}.snpeff.ann_filter.transitions.count
	"""
	
}

process UTILS_summary{

	publishDir("${params.outdir}/06_summaries")
	input:
	path(libsize_counts)
	path(trimmomatic_counts)
	path(rawaln_count)
	path(filteredaln_count)
	
	output:
	path("libsize_aln.tsv")
	
	
	script:
	"""
	sort  ${trimmomatic_counts} | (grep -v ^# || true)|  cut -f4-> trimmomatic_counts.tmp
	TRIMMOMATIC_HEADER=`grep -h -m 1 ^# ${trimmomatic_counts} | uniq | sed 's/\t/\tTRIMMOMATIC /g' | cut -f4-`;
	echo "#sample\tgroup\tgenotype\treadCount\tbases\tcoverage\t\$TRIMMOMATIC_HEADER\trawAln\tfilteredAln\tfilteredAlnCoverage" >  libsize_aln.tsv
	sort ${libsize_counts} > size.tmp
	sort ${rawaln_count}  | cut -f 4  > aln.tsv.raw
	sort ${filteredaln_count} | cut -f 4- > aln.tsv.fil
	paste size.tmp trimmomatic_counts.tmp aln.tsv.raw  aln.tsv.fil >>  libsize_aln.tsv
	
	"""
}

process UTILS_reformatTrimmomaticSummary {

	input:
	tuple val(sample), val(genotype), path(summary),val(group)
	
	output:
	path("${sample}.trimmed.summary.tsv")
	
	script:
	"""
	HEADER=`sed -E 's/: .+//' ${summary} | tr '\n' '\t'`
	SUMMARY=`sed -E 's/.+: //' ${summary} | tr '\n' '\t'`
	echo -e "#sample\tgroup\tgenotype\t\$HEADER" > ${sample}.trimmed.summary.tsv
	echo -e "${sample}\t${group}\t${genotype}\t\$SUMMARY" >> ${sample}.trimmed.summary.tsv
	"""
}

process UTILS_concatvcf {
	
	input:
	tuple val(sample), path(snp_indel_vcf),val(genotype)
	
	output:
	tuple val(sample), path("${sample}.raw.vcf.gz"),val(genotype)
	
	script:
	"""
	(zgrep -h '^#' *.snp.vcf.gz || true) > ${sample}.raw.vcf
	(zgrep -h -v '^#' ${snp_indel_vcf} || true) | sort -k1,1 -k2n,2 >> ${sample}.raw.vcf
	gzip ${sample}.raw.vcf
	"""
	
}


process UTILS_publishRawSnpEff {

	publishDir("${params.outdir}/03_vcf/merge")
	
	input:
	tuple val(sample),path(vcffile),path(tbi)
	
	output:
	path(vcffile)
	path(tbi)
	
	script:
	"""
	touch ${vcffile}
	touch ${tbi}
	"""
}

process UTILS_publishNoWTSnpEff{

	publishDir("${params.outdir}/03_vcf/merge/noWT/")
	
	input:
	tuple val(sample),path(vcffile),path(tbi)
	
	output:
	path(vcffile)
	
	script:
	"""
	touch ${vcffile}
	touch ${tbi}
	"""
}


process UTILS_publishANNFilterSnpEff{

	publishDir("${params.outdir}/03_vcf/merge/noWT/ANN_filter")
	
	input:
	tuple val(sample),path(vcffile),path(tbi)
	
	output:
	path(vcffile)
	path(tbi)
	
	script:
	"""
	touch ${vcffile}
	touch ${tbi}
	"""
}


