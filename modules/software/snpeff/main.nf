nextflow.enable.dsl=2

process SNPEFF_builddb {

	storeDir "store/"
	
	input:
	path(refseq)
	path(refannot)

	output:
	path("${refannot.baseName}-snpeff-db")

	script:
	"""
	mkdir -p ${refannot.baseName}-snpeff-db/data/genomes
	mkdir -p ${refannot.baseName}-snpeff-db/data/refannot
	cp $refseq ${refannot.baseName}-snpeff-db/data/genomes/refannot.fa
	cp $refannot ${refannot.baseName}-snpeff-db/data/refannot/genes.gff
	echo "refannot.genome : refannot" > ${refannot.baseName}-snpeff-db/snpeff.config
	snpEff build -gff3 -noCheckProtein -noCheckCds -c ${refannot.baseName}-snpeff-db/snpeff.config refannot
	"""
}

process SNPEFF_eff {

	input:
	path(snpeff_db)
	tuple val(sample), path(vcf)

	output:
	tuple val(sample), path("${sample}.raw.vcf.gz"),path("${sample}.raw.vcf.gz.tbi")

	script:
	snpeff_options = "-upDownStreamLen ${params.snpeff_up_downstream_len} -no-downstream  "

	"""
	snpEff eff ${snpeff_options}  -c ${snpeff_db}/snpeff.config refannot ${vcf} | bgzip -c > ${sample}.raw.vcf.gz
	tabix ${sample}.raw.vcf.gz
	"""
}
