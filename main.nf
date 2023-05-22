nextflow.enable.dsl=2

include { TRIMMOMATIC_trimmomatic } from "./modules/software/trimmomatic/main.nf"
include { BWA_index; BWA_mem } from "./modules/software/bwa/main.nf"
include { SAMTOOLS_filter; SAMTOOLS_mpileup; SAMTOOLS_countaln} from "./modules/software/samtools/main.nf"
include { VARSCAN_mpileup2snp;VARSCAN_mpileup2indel;VARSCAN_mpileup2cns;VARSCAN_mpileup2cns as VARSCAN_mpileup2cnsWT } from "./modules/software/varscan/main.nf"
include { VCFTOOLS_excludePostions;VCFTOOLS_convertToTsv as VCFTOOLS_convertToTsvRaw; VCFTOOLS_convertToTsv as VCFTOOLS_convertToTsvNoWT; VCFTOOLS_convertToTsv as VCFTOOLS_convertToTsvNoWTAnn } from "./modules/software/vcftools/main.nf"
include { SNPEFF_builddb;SNPEFF_eff; SNPEFF_eff as SNPEFF_effWT} from "./modules/software/snpeff/main.nf"
include { BCFTOOLS_stats as BCFTOOLS_stats_WT;BCFTOOLS_stats as BCFTOOLS_stats_noWT; BCFTOOLS_stats as BCFTOOLS_stats_noWT_ANNfilter; BCFTOOLS_stats as BCFTOOLS_stats_raw;BCFTOOLS_getDeletions;BCFTOOLS_snpeff_exclude; BCFTOOLS_merge as BCFTOOLS_mergeRawMutants; BCFTOOLS_merge as BCFTOOLS_mergeRawWT; BCFTOOLS_norm as BCFTOOLS_normRawMutants; BCFTOOLS_norm as BCFTOOLS_normRawWT} from "./modules/software/bcftools/main.nf"
include { UTILS_publishANNFilterSnpEff;UTILS_publishNoWTSnpEff;UTILS_publishRawSnpEff;UTILS_publishRawSnpEff as UTILS_publishRawSnpEffWT;UTILS_concatvcf;UTILS_countSnps as UTILS_countWTSnps;UTILS_countSnps as  UTILS_countMUTANTSnps;  UTILS_getSnpPositions as UTILS_getMUTANTSnpPositions;UTILS_getSnpPositions as UTILS_getWTSnpPositions; UTILS_genomeSize; UTILS_libSize; UTILS_rawaln; UTILS_summary;UTILS_reformatTrimmomaticSummary} from "./modules/software/utils/main.nf"
include { UTILS_mpileupRecall; UTILS_mpileupRecall as UTILS_mpileupRecallWT } from "./modules/software/utils/main.nf"
include { BCFTOOLS_fixref as BCFTOOLS_fixrefMutants; BCFTOOLS_fixref as BCFTOOLS_fixrefRawWT;BCFTOOLS_convertToTsvBySample } from "./modules/software/bcftools/main.nf"

workflow {

	Utils.assertFileExists('design',params.design)

	sample_ch = Channel
    .fromPath(params.design)
    .splitCsv(header: true)
    .map { row -> tuple (row.sample, file(row.fastq_1), file(row.fastq_2), row.genotype,row.readlen, row.group)}

	refseq_ch = Channel.fromPath(params.refseq).collect()
	refannot_ch = Channel.fromPath(params.refannot).collect()
	
	File cfgDir = new File(params.outdir + "/01_cfg")
	cfgDir.mkdir()
	Utils.dumpParams("${params.outdir}/01_cfg/run.parameters.cfg",params)
	
	genomesize = UTILS_genomeSize(refseq_ch).collect().flatten()
	snpeffdb_ch = SNPEFF_builddb(refseq_ch,refannot_ch).collect()
	bwa_index = BWA_index(refseq_ch).collect()

	libsize_ch = UTILS_libSize(sample_ch.combine(genomesize))

	if ( params.bamfiles ) {
		Utils.assertFileExists('bamfiles',params.bamfiles)
	
		bamfiles_ch = Channel
			.fromPath(params.bamfiles)
			.splitCsv(header: true)
			.map { row -> tuple (row.sample, file(row.bamfile))}
	
		bam_ch = sample_ch.concat(bamfiles_ch)
			.groupTuple()
			.map { t -> tuple(t.get(0), t.get(1).last(), t.get(3).first(),t.get(5).first()) }

		trimmomatic_summary_ch = Channel.empty();
		samgz_ch = Channel.empty();
		
	}
	else {
		// CLEANING
		(trimmed_ch,trimmomatic_summary_ch) = TRIMMOMATIC_trimmomatic(sample_ch)
		
		// MAPPING
		samgz_ch = BWA_mem(bwa_index,trimmed_ch)
		
		//MAPPING FILTER
		bam_ch = SAMTOOLS_filter(refseq_ch,samgz_ch)
	}
	
	// SNP CALLING BY LIB
	mpileup_ch = SAMTOOLS_mpileup(refseq_ch,bam_ch.groupTuple(by:[3]))
	snps_ch = VARSCAN_mpileup2snp(mpileup_ch)
	indels_ch = VARSCAN_mpileup2indel(mpileup_ch)
	vcf_ch = UTILS_concatvcf(snps_ch.concat(indels_ch).groupTuple(by:[0,2]))
	
	// SNP RECALL ON WT
	wtsample_ch = vcf_ch.filter {it.get(2) == params.wt_genotype}
	wt_snp_positions = UTILS_getWTSnpPositions(wtsample_ch.map{ t -> t.get(1)}.flatten().collect(),'WT', params.min_line_wt_position,999999).collect()
	//mpileup_wt_ch = SAMTOOLS_mpileupRecallWT(wt_snp_positions,refseq_ch,bam_ch.filter {it.get(2) == params.wt_genotype}.groupTuple(by:[3]))
	mpileup_wt_ch = UTILS_mpileupRecallWT(wt_snp_positions,mpileup_ch.filter {it.get(2) == params.wt_genotype})
	cns_wt_ch = VARSCAN_mpileup2cnsWT(mpileup_wt_ch)
	fixref_cns_wt_ch = BCFTOOLS_fixrefRawWT(cns_wt_ch,refseq_ch)
	
	merged_raw_wt_vcf = BCFTOOLS_mergeRawWT(fixref_cns_wt_ch.collect(), 'WT','raw')
	merged_wt_vcf = BCFTOOLS_normRawWT(merged_raw_wt_vcf,refseq_ch,'WT')
	
	
	vcfSnpeff_wt_ch = SNPEFF_effWT(snpeffdb_ch,merged_wt_vcf.collect())
	UTILS_publishRawSnpEffWT(vcfSnpeff_wt_ch)

	// SNP RECALL ON MUTANTS
	mutantsample_ch = vcf_ch.filter {it.get(2) != params.wt_genotype}
	mutant_snp_positions = UTILS_getMUTANTSnpPositions(mutantsample_ch.map{ t -> t.get(1)}.flatten().collect(),'MUTANT',params.min_line_mutant_position,params.max_line_mutant_position).collect()
	//mpileup_mutants_ch = SAMTOOLS_mpileupRecall(mutant_snp_positions,refseq_ch,bam_ch.filter {it.get(2) != params.wt_genotype}.groupTuple(by:[3]))
	mpileup_mutants_ch = UTILS_mpileupRecall(mutant_snp_positions,mpileup_ch.filter {it.get(2) != params.wt_genotype})
	
	
	cns_mutants_ch = VARSCAN_mpileup2cns(mpileup_mutants_ch)
	fixref_cns_mutants_ch = BCFTOOLS_fixrefMutants(cns_mutants_ch,refseq_ch)
	merged_raw_mutants_vcf = BCFTOOLS_mergeRawMutants(fixref_cns_mutants_ch.collect(),'MUTANTS','raw')
	merged_mutants_vcf = BCFTOOLS_normRawMutants(merged_raw_mutants_vcf,refseq_ch,'MUTANTS')
	vcfSnpeff_mutants_ch = SNPEFF_eff(snpeffdb_ch,merged_mutants_vcf.collect())
	UTILS_publishRawSnpEff(vcfSnpeff_mutants_ch)


	// SNP FILTERING ON MUTANTS
	// Remove positions found in WT
	(noWTvcf_ch,removedSnps_ch) = VCFTOOLS_excludePostions(wt_snp_positions,vcfSnpeff_mutants_ch)
	UTILS_publishNoWTSnpEff(noWTvcf_ch)
	
	BCFTOOLS_convertToTsvBySample(noWTvcf_ch.combine(mpileup_mutants_ch).map { t -> tuple(t.get(3), t.get(1), t.get(2)) })
	
	// Remove positions based on SnpEff ANN
	(noWT_ANNfilter_vcf_ch, snps_count_ch) = BCFTOOLS_snpeff_exclude(noWTvcf_ch)
	UTILS_publishANNFilterSnpEff(noWT_ANNfilter_vcf_ch)
	
	
	// Convert to matrix
	VCFTOOLS_convertToTsvRaw(vcfSnpeff_mutants_ch)
	VCFTOOLS_convertToTsvNoWT(noWTvcf_ch)
	VCFTOOLS_convertToTsvNoWTAnn(noWT_ANNfilter_vcf_ch)
	
	
	// STATS
	
	trimmomatic_count_ch = UTILS_reformatTrimmomaticSummary(trimmomatic_summary_ch)
	rawaln_count_ch = UTILS_rawaln(samgz_ch)
	filteredaln_count_ch = SAMTOOLS_countaln(bam_ch)
	(mapping,wtcalling,othercalling) = UTILS_summary(libsize_ch.collect(),trimmomatic_count_ch.collect(),rawaln_count_ch.collect(),filteredaln_count_ch.collect())
	
	BCFTOOLS_stats_WT(vcfSnpeff_wt_ch)
	BCFTOOLS_stats_noWT_ANNfilter(noWT_ANNfilter_vcf_ch)
	BCFTOOLS_stats_noWT(noWTvcf_ch)
	BCFTOOLS_stats_raw(vcfSnpeff_mutants_ch)
	
}


