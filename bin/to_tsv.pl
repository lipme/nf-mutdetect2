#!/usr/bin/perl


use strict;

MAIN:
{
	
	my $vcf = $ARGV[0];
	my $sample = $ARGV[1];
	
	die "Sample must be defined - $0 <VCF> <SAMPLE>" if ($sample eq '');
	die "VCF file not fount - $0 <VCF> <SAMPLE>" if (! -e $vcf);
	
	my @a_doc = ('ADP="Average per-sample depth of bases with Phred score = 30"',
'WT="Number of samples called reference (wild-type)"',
'HET="Number of samples called heterozygous-variant"',
'HOM="Number of samples called homozygous-variant"',
'NC="Number of samples not called"',
'GT="Genotype"',
'GQ="Genotype Quality"',
'SDP="Raw Read Depth as reported by SAMtools"',
'DP="Quality Read Depth of bases with Phred score = 30"',
'RD="Depth of reference-supporting bases (reads1)"',
'AD="Depth of variant-supporting bases (reads2)"',
'FREQ="Variant allele frequency"',
'PVAL="P-value from Fisher Exact Test"',
'RBQ="Average quality of reference-supporting bases (qual1)"',
'ABQ="Average quality of variant-supporting bases (qual2)"',
'RDF="Depth of reference-supporting bases on forward strand (reads1plus)"',
'RDR="Depth of reference-supporting bases on reverse strand (reads1minus)"',
'ADF="Depth of variant-supporting bases on forward strand (reads2plus)"',
'ADR="Depth of variant-supporting bases on reverse strand (reads2minus)"',
'VariationType="SNP/INS/DEL"',
'VariationCode="REF_ALT"',
'Putative_impact_Level="HIGH=3, MODERATE=2, LOW=1, MODIFIER=0}"',
'Allele="ALT"',
'Annotation="Effect Annotated using Sequence Ontology terms"',
'Putative_impact="A simple estimation of putative impact / deleteriousness : {HIGH, MODERATE, LOW, MODIFIER}"',
'Gene Name="Common gene name (HGNC). Optional: use closest gene when the variant is intergenic."',
'Gene ID="Gene ID"',
'Feature type="Which type of feature is in the next field (e.g. transcript, motif, miRNA, etc.)"',
'Feature ID="Depending on the annotation, this may be: Transcript ID (preferably using version number), Motif ID, miRNA, ChipSeq peak, Histone mark, etc. Note: Some features may not have ID (e.g. histone marks from custom Chip-Seq experiments may not have a unique ID)."',
'Transcript biotype="The bare minimum is at least a description on whether the transcript is {Coding, Noncoding}"',
'Rank / total="Exon or Intron rank / total number of exons or introns."',
'HGVS.c="Variant using HGVS notation (DNA level)"',
'HGVS.p="If variant is coding, this field describes the variant using HGVS notation (Protein level). Since transcript ID is already mentioned in feature ID, it may be omitted here."',
'cDNA_position / cDNA_len="Position in cDNA and trancripts cDNA length (one based)."',
'CDS_position / CDS_len="Position and number of coding bases (one based includes START and STOP codons)."',
'Protein_position / Protein_len="Position and number of AA (one based, including START, but not STOP)."',
'Distance to feature="All items in this field are options, so the field could be empty;  Up/Downstream: Distance to first / last codon; Intergenic: Distance to closest gene; Distance to closest Intron boundary in exon (+/- up/downstream). If same, use positive number.; Distance to closest exon boundary in Intron (+/- up/downstream); Distance to first base in MOTIF;   Distance to first base in miRNA; Distance to exon-intron boundary in splice_site or splice _region; ChipSeq peak: Distance to summit (or peak center); Histone mark / Histone state: Distance to summit (or peak center);"',
'Errors="Warnings or Information messages"');

	my $header = '%CHROM\t%POS\t%REF\t%ALT\t%ADP\t%WT\t%HET\t%HOM\t%NC\t[%SAMPLE\t%GT\t%GQ\t%SDP\t%DP\t%RD\t%AD\t%FREQ\t%PVAL\t%RBQ\t%ABQ\t%RDF\t%RDR\t%ADF\t%ADR]\t%ANN\n';
	my @a_ann_fields = ('VariationType','VariationCode', 'Putative_impact_Level','Allele','Annotation','Putative_impact','Gene Name','Gene ID','Feature type','Feature ID','Transcript biotype','Rank / total','HGVS.c','HGVS.p','cDNA_position / cDNA_len','CDS_position / CDS_len','Protein_position / Protein_len','Distance to feature','Errors, Warnings or Information messages');
	my $ann_header = join ("\t",@a_ann_fields); 
	
	
	my $cmd = "bcftools query -s $sample -f  '$header' $vcf | grep -P '\\t[01]/1\\t'";
	print STDERR "[CMD]\t$cmd\n";
	my @a_lines = `$cmd`;
	
	my $tsv_header = $header;
	$tsv_header =~ s/\\t/\t/g;
	$tsv_header =~ s/%//g;
	$tsv_header =~ s/ANN\\n//g;
	$tsv_header =~ s/\[//g;
	$tsv_header =~ s/\]//g;
	
	chomp @a_lines;
	print "##$sample\n##";
	print join("\n##", @a_doc),"\n";
	print "#$tsv_header$ann_header\n" ;
	
	foreach my $l (@a_lines)
	{
			#print STDERR "[INFO]\t$l\n";
			my @a_f = split ("\t", $l, -1);
			
			my $type = 'SNP';
			my $code = $a_f[2] . '_' . $a_f[3];
			if (length $a_f[2] > 1) {
				$type = 'DEL';
			} elsif (length $a_f[3] > 1 ) {
				$type = 'INS';
			}
			
			
			
			my $ann = pop(@a_f);
			
			my @a_ann = split (',', $ann);
			foreach my $ann_f (@a_ann)
			{
				my @a_ann_values = split ('\|', $ann_f, scalar @a_ann_fields);
				my $putative_impact_level = 0;
				
				if ($a_ann_values[2] eq 'HIGH') { 
					$putative_impact_level = 3;
				} elsif ($a_ann_values[2] eq 'MODERATE') { 
					$putative_impact_level = 2;
				} elsif ($a_ann_values[2] eq 'LOW') { 
					$putative_impact_level = 1;
				}
				print join("\t", @a_f), "\t" , $type, "\t",  $code,  "\t",$putative_impact_level , "\t", join ("\t", @a_ann_values) , "\n";
			}
	}

	exit 0;
}


	
=pod

CREATE TABLE VARIATION(
CHROM TEXT,
POS NUMERIC,
REF TEXT,
ALT TEXT,
ADP NUMERIC,
WT NUMERIC,
HET NUMERIC,
HOM NUMERIC,
NC NUMERIC,
SAMPLE TEXT,
GT TEXT,
GQ NUMERIC,
SDP NUMERIC,
DP NUMERIC,
RD NUMERIC,
AD NUMERIC,
FREQ NUMERIC,
PVAL NUMERIC,
RBQ NUMERIC,
ABQ NUMERIC,
RDF NUMERIC,
RDR NUMERIC,
ADF NUMERIC,
ADR NUMERIC,
VARIATIONTYPE TEXT,
VARIATIONCODE TEXT, 
PUTATIVE_IMPACT_LEVEL NUMERIC, 
ALLELE TEXT,
ANNOTATION TEXT,
PUTATIVE_IMPACT TEXT,
GENE_NAME TEXT,
GENE_ID TEXT,
FEATURE_TYPE TEXT,
FEATURE_ID TEXT,
TRANSCRIPT_BIOTYPE TEXT,
RANK_TOTAL TEXT,
HGVS_C TEXT,
HGVS_P TEXT,
CDNA_POSITION_CDNA_LEN TEXT,
CDS_POSITION_CDS_LEN TEXT,
PROTEIN_POSITION_PROTEIN_LEN TEXT,
DISTANCE_TO_FEATURE TEXT,
ERRORS TEXT);

CREATE INDEX SAMPLE_INDEX ON VARIATION (SAMPLE);
CREATE INDEX GENE_INDEX ON VARIATION (GENE_NAME);
CREATE INDEX EFFECT_INDEX ON VARIATION (ANNOTATION);
CREATE INDEX IMPACT_INDEX ON VARIATION (PUTATIVE_IMPACT);
CREATE INDEX IMPACT_INDEX ON VARIATION (PUTATIVE_IMPACT_LEVEL);
CREATE INDEX VARIATIONTYPE_INDEX ON VARIATION (VARIATIONTYPE);
CREATE INDEX VARIATIONCODE_INDEX ON VARIATION (VARIATIONCODE);

	## sqlite3 -header  -cmd ".mode tab"

1. Quels sont les gènes qui sont mutés dans au moins deux lignées parmi les XXX lignées listees

	"select * from VARIATION where SAMPLE IN ('$sample_list')  AND PUTATIVE_IMPACT_LEVEL >= $level"

2. what are the mutations affecting a gene in the line ES1M5S10185?
3. what are the mutations  in the line ES1M5S10185?
4. quels sont les mutants présentant un polymorphisme dans le gène AT1G01020


	"select * from VARIATION where GENE_NAME == '$gene'  AND PUTATIVE_IMPACT_LEVEL >= $level"

5. quels sont les mutants présentant des polymorphismes sur la région Chr1:5000-10000

	"select * from VARIATION where CHROM == '$chrom' AND POS >= $start AND POS <= $end "

6. quels sont les gènes mutés chez ES1M5S10185
	
	"select * from VARIATION where SAMPLE == '$sample'  AND PUTATIVE_IMPACT_LEVEL >= $level"

7. Quels sont les genes mutés en commun entre les mutants XXXX (liste de lignées)

	"select * from VARIATION where SAMPLE IN ('$sample_list')  AND PUTATIVE_IMPACT_LEVEL >= $level"


=cut
