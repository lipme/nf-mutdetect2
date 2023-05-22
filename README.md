# NF-MUTDETECT2

A Nextflow pipeline to detect EMS polymorphic positions in a mutant library (Illumina paired end sequencing).
Inspired by: https://github.com/ijpb-bioinformatics/mutdetect-pipeline

## DESCRIPTION

1. trim raw reads based on quality [trimmomatic]
2. map reads [bwa-mem]
3. filter alignments [samtools]
4. call SNP/INDELS for each sample (technical replicates are pooled) [samtools mpileup | varscan]
5. WT processing
5.1. get all WT polymorphic positions 
5.2. recall positions on all WT samples [varscsan mpileup2cns]
5.3. predict polymorphism effect on WT samples [snpeff]
6. MUTANT processing
6.1. get all MUTANT polymorphic positions 
6.2. recall positions on all MUTANT samples [varscsan mpileup2cns]
6.3. predict polymorphism effect on MUTANT samples [snpeff]
7. MUTANT filtering
7.1. Exclude WT polymorphisms (from 5.1)
7.2. Exclude MUTANTS polymorphisms found in more than `max_line_mutant_position`
7.2. Exclude non usefull effect positions
8. Produces summary files

## RUN

### PRE-REQUISITES

- NextFlow (> 21.04.0)  (Java 8 or later is required) : <https://www.nextflow.io/index.html#GetStarted>
- Singularity container app (> 3.0) : <https://sylabs.io/guides/3.0/user-guide/index.html>

### USAGE
```bash
nextflow run $DIR/nf-mutdetect2/main.nf \
	--design $DIR/nf-mutdetect2/test.csv \
	--outdir $DIR/nf-mutdetect2/test \
	--refseq $DIR/reference/AtCol0-TAIR10.genome.fasta \
	--refannot $DIR/reference/AtCol0-TAIR10.gff3 -profile genotoul,debug -resume -ansi-log false
```

## OUTPUT
```bash
├── 01_cfg
│   └── run.parameters.cfg
├── 02_bam
│   └── filtered #Filtered by sample alignments
├── 03_vcf
│   └── merge
│       ├── MUTANTS.raw.vcf.gz # Raw polymorphic sites in MUTANT samples in VCF format
│       ├── MUTANTS.raw.vcf.gz.tbi
│       ├── noWT
│       │   ├── MUTANTS.noWT.vcf.gz # Non WT polymorphic sites in MUTANT samples in VCF format
│       │   └── MUTANTS.noWT.vcf.gz.tbi
│       │   └── ANN_filter
│       │       ├── MUTANTS.noWT.ann_filter.vcf.gz # Non WT and ANN filtered polymorphic sites in MUTANT samples in VCF format
│       │       └── MUTANTS.noWT.ann_filter.vcf.gz.tbi
│       ├── WT.raw.vcf.gz
│       └── WT.raw.vcf.gz.tbi
├── 04_positions
│   ├── MUTANT.snpPositions.tsv # MUTANT polymorphic sites (compliant with min number of lines)
│   └── WT.snpPositions.tsv  # WT polymorphic sites (compliant with min number of lines)
├── 05_tsv
│   ├── bySample
│   │      └── XXXXXX.tsv.gz 
│   ├── MUTANTS.raw.vcf.tsv # Raw polymorphic sites in MUTANT samples table
│   ├── MUTANTS.noWT.vcf.tsv  # Non WT polymorphic sites in MUTANT samples table
│   └── MUTANTS.noWT.ann_filter.vcf.tsv # Non WT and ANN filtered polymorphic sites in MUTANT samples table
├── 06_summaries
│   ├── libsize_aln.tsv # Raw data and alignment statistics
│   ├── MUTANTS.raw.vcf.stats.tsv # Statistics about raw polymorphic sites in MUTANT samples 
│   ├── MUTANTS.noWT.vcf.stats.tsv # Statistics about WT polymorphic sites in MUTANT samples 
│   ├── MUTANTS.noWT.ann_filter.vcf.stats.tsv # Statistics about non WT and ANN filtered polymorphic sites in MUTANT samples
│   ├── MUTANTS.raw.vcf.stats.perSample.tsv # Per sample statistics
│   ├── MUTANTS.noWT.vcf.stats.perSample.tsv
│   ├── MUTANTS.noWT.ann_filter.vcf.stats.perSample.tsv
│   ├── WT.raw.vcf.stats.tsv #WT statistics
│   └── WT.raw.vcf.stats.perSample.tsv #WT per sample statistics
└── 08_pipeline_info
    ├── execution_report.html
    ├── execution_timeline.html
    ├── execution_trace.txt
    └── pipeline_dag.svg
```
## PARAMETERS

### MANDATORY

- design : a CSV file describing your samples

```csv
sample,group,fastq_1,fastq_2,genotype,readlen
SEQ1,SEQ1,/path/to/SEQ1_R1.fastq.gz,/path/to/SEQ1_R2.fastq.gz,EMS,150
SEQ2,SEQ2,/path/to/SEQ2_R1.fastq.gz,/path/to/SEQ2_R2.fastq.gz,WT,150
SEQ2-2,SEQ2,/path/to/SEQ2-2_R1.fastq.gz,/path/to/SEQ2_2_R2.fastq.gz,WT,150
```

NB: samples from same group will be merged during SNP calling step

-refseq : reference genome FASTA file

-refannot : reference annotation GFF3 file

### ADVANCED

- bamfiles : a CSV file describing path to precomputed clean alignment files (already mapped and filtered, ready to SNPcalling)
```csv
sample,bamfile
SEQ1,/path/to/SEQ1_rmdup.bam
SEQ2,/path/to/SEQ2_rmdup.bam
SEQ2-2,/path/to/SEQ2-2_rmdup.bam
```

- Other default values are set in `conf/params/default.config`:

```java
  //genotype in design file corresponding to wildtype (parental lineage) samples
  wt_genotype = 'WT'
  
  // TRIMMOMATIC parameters
  trimmo_leading_qual_min = 20 // LEADING: Cut bases off the start of a read, if below a threshold quality
  trimmo_trailing_qual_min = 20 // TRAILING: Cut bases off the end of a read, if below a threshold quality
  trimmo_sliding_window_size = 4  // SLIDINGWINDOW: Performs a sliding window trimming approach. It starts scanning at the 
  trimmo_sliding_window_qual = 20 // 5‟ end and clips the read once the average quality within the window falls below a threshold.
  trimmo_min_length = 50 //MINLEN: Drop the read if it is below a specified length
  trimmo_params = ' -validatePairs -phred33 ' // CUSTOM parameters
  
  // SAMTOOLS Filtering parameters
  aln_filter_mode = 'paired_rmdup' //values : paired_rmdup, rmdup, paired
  samtools_max_memory = '8G'
  // MPILEUP parameters
  mpileup_parameters = ' -B --max-depth 100 '
  mpileuprecall_parameters = '-B'
  
  // VARSCAN MUTANT SAMPLE PARAMETERS
  varscan_mutant_min_coverage = 3 // Minimum read depth at a position to make a call
  varscan_mutant_min_reads2 = 3 // Minimum supporting reads at a position to call variants 
  varscan_mutant_avg_qual = 15 // Minimum base quality at a position to count a read
  varscan_mutant_var_freq = 0.2 // Minimum variant allele frequency threshold
  varscan_mutant_var_freq_for_hom = 0.8 // 	Minimum frequency to call homozygote
  varscan_mutant_pvalue = 0.01 // Default p-value threshold for calling variants
  
  // VARSCAN WILDTYPE SAMPLE PARAMETERS
  varscan_wt_min_coverage = 5
  varscan_wt_min_reads2 = 4
  varscan_wt_avg_qual = 15
  varscan_wt_var_freq = 0.2
  varscan_wt_var_freq_for_hom = 0.8
  varscan_wt_pvalue = 0.01
  
  // POSITION FILTERING PARAMETERS
  min_line_wt_position = 3 // Minimum number of WT samples to support a polymorphic position to be kept
  min_line_mutant_position = 1  // Minimum number of MUTANT samples to support a polymorphic position  to be kept
  max_line_mutant_position = 4 // Maximum number of MUTANT samples to support a polymorphic position  to be kept
  
  snpeff_up_downstream_len = 1000 // SNPEff downstream length 
  snpeff_exclude_ann = 'LOW\\|intergenic_variant\\|upstream_gene_variant\\|intergenic_region' // SNPEff ANN field filtering string
```
## USE CASE

This pipeline has been developped and applied to the <i>A. thaliana</i> HEM library (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6157545/)
Results are available here: https://lipm-browsers.toulouse.inrae.fr/pub/ATHEM

## CONTAINERS

This pipeline uses a set of Singularity/Apptainer containers.
Pre-built images are available here https://lipm-browsers.toulouse.inra.fr/pub/singularity-repository/mutdetect2/ and automatically downloaded.
Definition files are provided in the `containers` directory

## CONTACT
sebastien.carrere@inrae.fr

## CITATION

 - Sebastien Carrere, Patrick Laufs, Jean-Marc Routaboul, Raphael Mercier & Laurent Noel, TO BE PUBLISHED

## CREDITS

- Sebastien Carrere (INRAE/CNRS - LIPME)
