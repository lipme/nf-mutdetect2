# NF-MUTDETECT2

```bash
nextflow run $DIR/nf-mutdetect2/main.nf \
	--design $DIR/nf-mutdetect2/test.csv \
	--outdir $DIR/nf-mutdetect2/test \
	--refseq $DIR/reference/AtCol0-TAIR10.genome.fasta \
	--refannot $DIR/reference/AtCol0-TAIR10.gff3 -profile genotoul,debug -resume -ansi-log false
```

- test.csv

```tsv
sample,group,fastq_1,fastq_2,genotype,readlen
SEQ1,SEQ1,/path/to/SEQ1_R1.fastq.gz,/path/to/SEQ1_R2.fastq.gz,EMS,150
SEQ2,SEQ2,/path/to/SEQ2_R1.fastq.gz,/path/to/SEQ2_R2.fastq.gz,WR,150
```

* group: merge technical replicate for SNP calling step
