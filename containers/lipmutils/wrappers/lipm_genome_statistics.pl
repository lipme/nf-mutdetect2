#!/usr/bin/perl

# $Id: lipm_genome_statistics.pl 1387 2014-09-12 07:35:35Z sallet $

=pod

=head1 NAME

lipm_genome_statistics.pl

=head1  SYNOPSIS

./lipm_genome_statistics.pl --in "annotation/PLHAL710p0*.gff3" --out /PLHAL710_STAT --scaffolds PLHAL710.fna


=head1 DESCRIPTION

Get statistics about annotation. Load the GFF3 file(s) and generate two files:
    * out.general_statistics.xls  global statistics about genes/mRNA/CDS/UTR/intron/ncRNA/intergenic
    * out.statistics_per_gene.xls individual gene statistics


=cut

BEGIN
{
    $ENV{LANG} = "C";
}


use strict;
use General;
use GeneralBioinfo;
use ParamParser;
use Carp;


MAIN:
{
	my $o_param = New ParamParser('GETOPTLONG', \&Usage, 'help', 'in=s', 'out=s', 'scaffolds=s', 'kingdom=s');

	$o_param->AssertDefined('in', 'out');
	$o_param->AssertFileExists('scaffolds');
	$o_param->SetUnlessDefined('kingdom', "eukaryote");
	$o_param->AssertAllowedValue('kingdom', 'eukaryote', 'prokaryote');

	my $seq     = $o_param->Get('scaffolds');
	my $ra_gff3 = $o_param->GetListOfFiles('in');
	my $out     = $o_param->Get('out');
	my $kingdom = $o_param->Get('kingdom');

	my %h_sequences = ();
	&LoadSeq($seq, \%h_sequences);
	# initialize intergenic sequences with all the sequences
	my %h_intergenicsequences = ();
	&LoadSeq($seq, \%h_intergenicsequences);

	print "Input file(s):\n";
	foreach my $gff (@$ra_gff3)
	{
		print "   $gff\n";
	}

	my %h_genetype = ();
	my %h_stat     = ();
	my $nb_all_genes = 0; # gene of protein + gene of ncrna

	foreach my $gff (@$ra_gff3)
	{
		my ($dir, $f) = $gff =~ /(.+)\/(.+)/; # /the/in/path/file.gff3

		my $prev_prefixid = '';
		my $prev_type     = '';
		my %h_prev_res    = ();

		my $f_in = &GetStreamIn($gff);
		while(my $l = <$f_in>)
		{
			my %h_res = ();
			&ParseGFF3Line(\$l, 2, \%h_res);
			next if ($l =~ /^###/ || ! defined($h_res{ID}));
			my $seqid = $h_res{seqid};
			my $id    = $h_res{ID};
			my ($type, $prefixid, $altcode) = &SplitIds($id);
			# Consider only the first variant for statistics
			if ( ! defined($altcode) || $altcode eq 'a' )
			{
				# EXON
				if ($h_res{type} eq 'exon' )
				{
					my $s = &ExtractSeq (\$h_sequences{$seqid}, $h_res{start}, $h_res{end});
					$h_stat{$prefixid}{exon_gccount} += &GetGC($s);
					$h_stat{$prefixid}{exon_len}     += &GetLen($s);
					$h_stat{$prefixid}{exon_nb}++;
					if ( $prefixid eq $prev_prefixid && $h_prev_res{type} eq 'exon' ) # EK: peut etre plutot verifier le parent???
					{
						# overlapping exons 
						if ( $h_res{start} < $h_prev_res{end} )
						{
							Carp::croak("$h_res{start} < $h_prev_res{end} \n");
						}
						# Extract intron
						if ( ($h_res{start} -1) >= ($h_prev_res{end}+1) )
						{
							my $s = &ExtractSeq (\$h_sequences{$seqid}, $h_prev_res{end}+1, $h_res{start}-1);
							$h_stat{$prefixid}{intron_gccount} += &GetGC($s);
							$h_stat{$prefixid}{intron_len}     += &GetLen($s);
						}
					}
				}
				
				# CDS
				if ( $h_res{type} eq 'CDS' )
				{
					my $s = &ExtractSeq (\$h_sequences{$seqid}, $h_res{start},$h_res{end});
					$h_stat{$prefixid}{CDS_gccount} += &GetGC($s);
					$h_stat{$prefixid}{CDS_len}     += &GetLen($s);
					$h_stat{$prefixid}{CDS_len_withN} += length($s);
				}
				# mRNA
				if ($h_res{type} eq 'mRNA' )
				{
					$h_genetype{$prefixid} = 'mRNA';
					my $s = &ExtractSeq (\$h_sequences{$seqid}, $h_res{start}, $h_res{end});
					$h_stat{$prefixid}{mRNA_gene_len} = &GetLen($s); # length of the gene (not the mature mrna)
				}
				# ncRNA
				if ($h_res{type} eq 'ncRNA' || $h_res{type} eq 'tRNA' || $h_res{type} eq 'rRNA')
				{
					$h_genetype{$prefixid} = 'ncRNA';
					my $s = &ExtractSeq (\$h_sequences{$seqid}, $h_res{start}, $h_res{end});
					$h_stat{$prefixid}{ncRNA_gccount}  += &GetGC($s);
					$h_stat{$prefixid}{ncRNA_gene_len} += &GetLen($s);
				}
				# five_prime_UTR
				if ($h_res{type} eq 'five_prime_UTR')
				{
					my $s = &ExtractSeq (\$h_sequences{$seqid}, $h_res{start}, $h_res{end});
					$h_stat{$prefixid}{five_UTR_len}     += &GetLen($s);
					$h_stat{$prefixid}{five_UTR_gccount} += &GetGC($s);
					if ( $prefixid eq $prev_prefixid && $h_prev_res{type} eq 'five_prime_UTR' ) # peut etre plutot verifier le parent???,
					{
						# Extract intron UTR
						if ( ($h_res{start} -1) >= ($h_prev_res{end}+1) )
						{
							my $s = &ExtractSeq (\$h_sequences{$seqid}, $h_prev_res{end}+1, $h_res{start}-1);
							$h_stat{$prefixid}{five_UTR_intron_gccount} += &GetGC($s);
							$h_stat{$prefixid}{five_UTR_intron_len}     += &GetLen($s);
						}
					}
				}
				# three_prime_UTR
				if ($h_res{type} eq 'three_prime_UTR')
				{
					my $s = &ExtractSeq (\$h_sequences{$seqid}, $h_res{start}, $h_res{end});
					$h_stat{$prefixid}{three_UTR_len}     += &GetLen($s);
					$h_stat{$prefixid}{three_UTR_gccount} += &GetGC($s);
					if ( $prefixid eq $prev_prefixid && $h_prev_res{type} eq 'three_prime_UTR' ) # peut etre plutot verifier le parent???,
					{
						# Extract intron UTR
						if ( ($h_res{start} -1) >= ($h_prev_res{end}+1) )
						{
							my $s = &ExtractSeq (\$h_sequences{$seqid}, $h_prev_res{end}+1, $h_res{start}-1);
							$h_stat{$prefixid}{three_UTR_intron_gccount} += &GetGC($s);
							$h_stat{$prefixid}{three_UTR_intron_len}     += &GetLen($s);
						}
					}
				}
			}
	
			if ( $h_res{type} eq 'gene' )
			{
				$nb_all_genes++; # included ncrna gene

				my $s = &ExtractSeq (\$h_sequences{$seqid}, $h_res{start}, $h_res{end});
				$h_stat{$prefixid}{gene_len} = &GetLen($s);

				# replace the gene subsequence by 'X' in the intergenic sequence
				my $gene_len = $h_res{end} - $h_res{start} + 1; # il faut laisser les n
				substr($h_intergenicsequences{$seqid}, $h_res{start}-1, $gene_len, 'x' x $gene_len);
			}

			$h_stat{$prefixid}{nb_variants}++ if ( $h_res{type} eq 'mRNA' || $h_res{type} eq 'ncRNA' || $h_res{type} eq 'tRNA' || $h_res{type} eq 'rRNA');
			
			$prev_prefixid = $prefixid;
			%h_prev_res    = %h_res;
		}
	}

	my $f_stat    = &GetStreamOut("$out.statistics_per_gene.xls");
	my $f_general = &GetStreamOut("$out.general_statistics.xls");

	if ($kingdom eq "eukaryote")
	{
		print $f_stat "gene_EGN_ID\ttype\tmRNA_length\tExons_nb\tIntron_len\tVariant_nb";
		print $f_stat "\tfive_UTR_len\tfive_UTR_intron_len\tthree_UTR_len\tthree_UTR_intron_len\n";
	}
	else
	{
		print $f_stat "gene_EGN_ID\ttype\tmRNA_length\tVariant_nb";
		print $f_stat "\tfive_UTR_len\tthree_UTR_len\n";
	}

	# global count
	my $g_intron_nb      			= 0; # global intron nb	
	my $g_intron_len                = 0;
	my $g_exon_nb                   = 0; # global exon nb
	my $g_exon_len                  = 0, 
	my $g_mrna_gene_len             = 0; # length of the gene coding for protein
	my $g_cds_len       			= 0; # global cds len
	my $min_cds = 1000000;
	my $max_cds = 0;
	my $g_ncrna_gene_len            = 0; # lg of the ncrna gene
	my $g_ncrna_mature_len          = 0; # lg of the matured ncrna
	my $g_ncrna_intron_nb           = 0; # global intron ncrna nb
	my $g_ncrna_intron_len          = 0; 
	my $min_ncrna = 1000000;
	my $max_ncrna = 0;
	my $g_ncrna_exon_nb             = 0;
	my $gene_length       			= 0; # gene (ncrna and mrna) length

	my $g_five_UTR_len              = 0;
	my $g_five_UTR_nb               = 0;
	my $g_three_UTR_len             = 0;
	my $g_three_UTR_nb              = 0;

	my $nb_mRNA_genes  				 = 0;
	my $nb_mRNA_genes_with_introns   = 0;
	my $nb_mRNA_genes_with_five_UTR  = 0;
	my $nb_mRNA_genes_with_three_UTR = 0;
	my $nb_ncRNA_genes               = 0;
	my $nb_ncRNA_genes_with_introns  = 0;

	my $cds_gc       = 0;
	my $intron_gc    = 0;
	my $exon_gc      = 0;
	my $ncrna_gc     = 0; # gc% of the ncrna mature
	my $five_utr_gc  = 0;
	my $three_utr_gc = 0;

	foreach my $id (sort keys(%h_stat))
	{
		if ($h_genetype{$id} eq 'mRNA')
		{
			if  ( defined($h_stat{$id}{five_UTR_len}) )
			{
				$nb_mRNA_genes_with_five_UTR++;
			}
			else
			{
				$h_stat{$id}{five_UTR_len} = 0 ;
			}
			if (defined $h_stat{$id}{three_UTR_len})
			{
				$nb_mRNA_genes_with_three_UTR++;
			}
			else
			{
				$h_stat{$id}{three_UTR_len} = 0; 
			}
			
			$h_stat{$id}{intron_len}           = 0  if ( ! defined($h_stat{$id}{intron_len}   ) );
			$h_stat{$id}{five_UTR_intron_len}  = 0  if ( ! defined($h_stat{$id}{five_UTR_intron_len} ) );		
			$h_stat{$id}{three_UTR_intron_len} = 0  if ( ! defined($h_stat{$id}{three_UTR_intron_len} ) );
			
			if ($kingdom eq "eukaryote")
			{
				print $f_stat "$id\tmRNA\t$h_stat{$id}{mRNA_gene_len}\t$h_stat{$id}{exon_nb}\t$h_stat{$id}{intron_len}\t$h_stat{$id}{nb_variants}";
				print $f_stat "\t$h_stat{$id}{five_UTR_len}\t$h_stat{$id}{five_UTR_intron_len}\t$h_stat{$id}{three_UTR_len}\t$h_stat{$id}{three_UTR_intron_len}\n";
			}
			else
			{
				print $f_stat "$id\tmRNA\t$h_stat{$id}{mRNA_gene_len}\t$h_stat{$id}{nb_variants}";
				print $f_stat "\t$h_stat{$id}{five_UTR_len}\t$h_stat{$id}{three_UTR_len}\n";
			}

			# Update global counts
			$g_intron_nb   += ($h_stat{$id}{exon_nb}-1);			
			$g_intron_len  += $h_stat{$id}{intron_len};
			$g_exon_nb     += $h_stat{$id}{exon_nb};
			$g_mrna_gene_len += $h_stat{$id}{mRNA_gene_len}; 
			$g_exon_len    += $h_stat{$id}{exon_len};
			$g_cds_len     += $h_stat{$id}{CDS_len};
			$min_cds       = $h_stat{$id}{CDS_len_withN} if ($h_stat{$id}{CDS_len_withN} > 0 && $h_stat{$id}{CDS_len_withN} < $min_cds);
			$max_cds       = $h_stat{$id}{CDS_len_withN} if ($h_stat{$id}{CDS_len_withN} > $max_cds);
			$cds_gc        += ($h_stat{$id}{CDS_gccount});
			$intron_gc     += ($h_stat{$id}{intron_gccount});
			$exon_gc       += ($h_stat{$id}{exon_gccount});
			$five_utr_gc   += $h_stat{$id}{five_UTR_gccount};
			$three_utr_gc  += $h_stat{$id}{three_UTR_gccount};
			
			$g_five_UTR_len  += $h_stat{$id}{five_UTR_len};
			$g_three_UTR_len += $h_stat{$id}{three_UTR_len};
		   
			$nb_mRNA_genes++;
			$nb_mRNA_genes_with_introns++ if ( $h_stat{$id}{exon_nb} > 1 );
		}
		elsif ($h_genetype{$id} eq 'ncRNA')
		{
			$h_stat{$id}{intron_len} = 0  if ( ! defined($h_stat{$id}{intron_len}  ) );
			$h_stat{$id}{exon_nb}    = 1  if ( ! defined($h_stat{$id}{exon_nb}     ) );

			# there are introns in ncrna
			if ( $h_stat{$id}{exon_nb} > 1 )
			{
				$ncrna_gc  += $h_stat{$id}{exon_gccount};
				$g_ncrna_intron_nb  += ($h_stat{$id}{exon_nb}-1); 			
				$g_ncrna_intron_len += $h_stat{$id}{intron_len};	
				$nb_ncRNA_genes_with_introns++;
			}
			else
			{
				$ncrna_gc  += ($h_stat{$id}{ncRNA_gccount});
			}
			
			if ($kingdom eq "eukaryote")
			{
				print $f_stat "$id\tncRNA\t$h_stat{$id}{ncRNA_gene_len}\t$h_stat{$id}{exon_nb}\t$h_stat{$id}{intron_len}\t$h_stat{$id}{nb_variants}";
				print $f_stat "\t-\t-\t-\t-\n";
			}
			else
			{
				print $f_stat "$id\tncRNA\t$h_stat{$id}{ncRNA_gene_len}\t$h_stat{$id}{nb_variants}";
				print $f_stat "\t-\t-\n";
				
			}
			$g_ncrna_gene_len += $h_stat{$id}{ncRNA_gene_len};
			$nb_ncRNA_genes++;
			$g_ncrna_exon_nb += $h_stat{$id}{exon_nb};
			$min_ncrna = $h_stat{$id}{ncRNA_gene_len} if ($h_stat{$id}{ncRNA_gene_len} < $min_ncrna);
			$max_ncrna = $h_stat{$id}{ncRNA_gene_len} if ($h_stat{$id}{ncRNA_gene_len} > $max_ncrna);
		}
		else
		{
			# unknow type!!!
		}
		$gene_length  += $h_stat{$id}{gene_len}; # gene (ncrna and mrna) length
	}
	$f_stat->close;

	my $seq_gc  = 0;
	my $seq_len = 0;
	foreach my $seqid (keys(%h_sequences))
	{
		$seq_gc   += &GetGC($h_sequences{$seqid});
		$seq_len  += &GetLen($h_sequences{$seqid});
	}


	printf $f_general "Number of nucleotides (without 'N')\t$seq_len\n";
	printf $f_general "\tPer cent GC\t%.2f\n", $seq_gc*100/$seq_len;
	printf $f_general "Total number of genes\t " . ($nb_mRNA_genes+$nb_ncRNA_genes) . "\n";
	printf $f_general "\tTotal nucleotides (bp)\t%d\n",  $gene_length;


    printf $f_general "\n** Protein coding genes\n";

	printf $f_general "Number of protein coding genes\t$nb_mRNA_genes\n";
	printf $f_general "\tMean gene length (bp)\t%.2f\n", $g_mrna_gene_len/$nb_mRNA_genes if ($nb_mRNA_genes > 0);
	printf $f_general "\tCoding nucleotides (bp)\t%d\n", $g_cds_len;
	if ($nb_mRNA_genes > 0)
	{
		printf $f_general "\tPer cent genes with introns\t%2.f\n",   $nb_mRNA_genes_with_introns  *100/$nb_mRNA_genes if ($kingdom eq "eukaryote");
		printf $f_general "\tPer cent genes with five UTR\t%2.f\n",  $nb_mRNA_genes_with_five_UTR *100/$nb_mRNA_genes;
		printf $f_general "\tPer cent genes with three UTR\t%2.f\n", $nb_mRNA_genes_with_three_UTR*100/$nb_mRNA_genes;
	}
	
	printf $f_general "Exons\n";
	if ($g_exon_nb > 0)
	{
		printf $f_general "\tMean number per gene\t%.2f\n",  $g_exon_nb/$nb_mRNA_genes if ($kingdom eq "eukaryote");
		printf $f_general "\tMean length (bp)\t%.2f\n",      $g_exon_len/$g_exon_nb;
		printf $f_general "\tGC per cent\t%.2f\n",           $exon_gc*100/$g_exon_len;
	}
	if ($kingdom eq "eukaryote")
	{
		printf $f_general "Introns\n";
		if ($g_intron_nb > 0)
		{
			printf $f_general "\tMean number per gene\t%.2f\n",  $g_intron_nb/$nb_mRNA_genes; 
			printf $f_general "\tMean length (bp)\t%.2f\n",      $g_intron_len/$g_intron_nb;
			printf $f_general "\tGC per cent\t%.2f\n",           $intron_gc*100/$g_intron_len;
		}
	}

	printf $f_general "CDS\n";
	if ($g_cds_len)
	{
		printf $f_general "\tMean length (bp)\t%.2f\n",      $g_cds_len/$nb_mRNA_genes;
		printf $f_general "\tMin length (bp)\t%.2f\n",      $min_cds;
		printf $f_general "\tMax length (bp)\t%.2f\n",      $max_cds;
		printf $f_general "\tGC per cent\t%.2f\n",           $cds_gc*100/$g_cds_len;
	}

	printf $f_general "five_prime_UTR\n";
	if ($nb_mRNA_genes_with_five_UTR > 0)
	{
		printf $f_general "\tMean length (bp)\t%.2f\n",      $g_five_UTR_len/$nb_mRNA_genes_with_five_UTR;
		printf $f_general "\tGC per cent\t%.2f\n",           $five_utr_gc*100/$g_five_UTR_len;
	}

	printf $f_general "three_prime_UTR\n";
	if ($nb_mRNA_genes_with_three_UTR > 0)
	{
		printf $f_general "\tMean length (bp)\t%.2f\n",      $g_three_UTR_len/$nb_mRNA_genes_with_three_UTR;
		printf $f_general "\tGC per cent\t%.2f\n",           $three_utr_gc*100/$g_three_UTR_len;
	}


	printf $f_general "\n** Non protein coding genes\n";
	printf $f_general "Number of non protein coding genes\t$nb_ncRNA_genes\n";
	if ($nb_ncRNA_genes > 0)
	{
		printf $f_general "\tMean ncRNA gene length (bp)\t%.2f\n", $g_ncrna_gene_len/$nb_ncRNA_genes;
		printf $f_general "\tMin length (bp) %i\n",  $min_ncrna;
		printf $f_general "\tMax length (bp) %i\n",  $max_ncrna;
		printf $f_general "\tGC per cent\t%.2f\n",           $ncrna_gc*100/$g_ncrna_gene_len;
		if ($kingdom eq "eukaryote")
		{
			printf $f_general "\tPer cent ncRNA genes with introns\t%2.f\n", $nb_ncRNA_genes_with_introns*100/$nb_ncRNA_genes;
			printf $f_general "\tMean exon number per ncRNA gene\t%.2f\n",  $g_ncrna_exon_nb/$nb_ncRNA_genes;
		}
	}


	print $f_general "\n** Intergenic (inter protein-coding genes)\n";
	my $ig_len = 0;
	my $ig_gc  = 0;
	foreach my $seqid (keys(%h_sequences))
	{
		$h_intergenicsequences{$seqid} =~ s/x//g; # x means covered by a gene. thus remove the x regions.
		$ig_gc  += &GetGC($h_intergenicsequences{$seqid});
		$ig_len += length($h_intergenicsequences{$seqid});
		
	}
	if ($nb_all_genes > 0)
	{
		printf $f_general "\tMean length\t%.2f\n",  $ig_len/$nb_all_genes;
	}
	printf $f_general "\tGC per cent\t%.2f\n",  $ig_gc*100/$ig_len;
	$f_general->close;
	
	exit 0;
}


=head2 procedure LoadSeq

 Title        : LoadSeq
 Usage        : 
 Prerequisite : none
 Function     : Fill $rh_seq. Key=> seq id, value=> sequence
 Args         : $file      Filename
                $rh_seq    Ref to a hashtable
 Access       : public
 Globals      : none 

=cut 
sub LoadSeq
{
	my($file, $rh_seq) = @_;
	
	my $f_seq = &GetStreamIn($file);
	my $id    = '';
	while(my $l=<$f_seq>)
	{
		chomp($l);
		if ( $l =~ /^>(\S+)/ )
		{
			$id = $1;
		}
		else
		{
			$$rh_seq{$id} .= $l;
		}
	}
	$f_seq->close;

	return;
}


=head2 function GetGC

 Title        : GetGC
 Usage        : my $nb = &GetGC($seq)
 Function     : return the number of 'G' or 'C' in the sequence
 Args         : $seq string
 Return       : the number of G or C in the string $seq
 Globals      : none

=cut
sub GetGC
{
	my($seq) = @_;
	
	my ($c) = $seq =~ tr/GCgc/GCgc/;
	
	return $c;
}

=head2 function GetLen

 Title        : GetLen
 Usage        : my $len = &GetLen($seq)
 Function     : return the length of the string without counting the 'N'
 Args         : $seq string
 Return       : the length of the string without counting the 'N'
 Globals      : none

=cut
sub GetLen
{
	my($seq) = @_;
	
	$seq =~ s/N//gi;
	
	return length($seq);
}


=head2 fonction SplitIds

 Title        : SplitIds
 Usage        : 
 Prerequisite : none
 Function     : Extract information from the ID of a feature: type, prefix, code alternatif
 Returns      : ($type, $prefixid, $altcode). Per example CDS, PLHAL710p01.1, "?"
 Args         : $id       string   ID of an element. Example CDS:PLHAL710p01.15.1
 Access       : public
 Globals      : none 

=cut 
sub SplitIds
{
	my($id) = @_;
	# exon:PLHAL710p07.24.1             --> PLHAL710p07.24 (enleve fin)
	# exon:GMI11495-Rm2011G.a.4.1       --> GMI11495-Rm2011G.a.4 (enleve fin)
	# exon:PLHAL710p07.24.utr0          --> PLHAL710p07.24 (enleve fin)
    # exon:Mt0001_00004.1               --> Mt0001_00004 (enleve la fin)
    # [gene|mRNA|ncRNA]:scaffold51_69.2 --> scaffold51_69.2 (tout)
	# mRNA:GMI11495-Rm2011G.a.4         --> GMI11495-Rm2011G.a.4 (tout)
	my @a_id     = split(/[:\.]/, $id);  
	my $type     = $a_id[0];              # exon 
	#my $prefixid = (scalar(@a_id)==4 || $type eq "gene" || $type eq "mRNA" || $type eq "ncRNA") ? "$a_id[1].$a_id[2]" : $a_id[1]; 
    #  scalar= 4 = exon:PLHAL710p07.24.utr0 ou 3 = exon:Mt0001_00004.1
	my ($altcode) = $id =~ /([a-z])$/;   # last letter of an alternative gene 

	my $prefixid;
	if ($type eq "gene" || $type eq "mRNA" || $type eq "ncRNA")
	{
		($prefixid) = $id =~ /.+:(.+)$/;
	}
	else
	{
		($prefixid) = $id =~ /.+:(.+)\.[^\.]+$/;
	}
	return ($type, $prefixid, $altcode);
}

sub Usage	
{
print STDERR<<END

Get statistics about annotation. Load the GFF3 file(s) and generate two files:
    * out.general_statistics.xls  global statistics about genes/mRNA/CDS/UTR/intron/ncRNA/intergenic
    * out.statistics_per_gene.xls individual gene statistics

------------------- USAGE --------------------
$0
[Mandatory]
    --help
    --in        file, file of files (.fof) or list of files   GFF3 file(s) of genes. Recognized elts: gene/mRNA/exon*CDS/ncRNA
    --out       filename                                      Root filename of output files. 
    --scaffolds fastafile                                     Multifasta of the corresponding genomic sequence(s)
    --kingdom   <eukaryote|prokaryote>                        Kingdom. Default=eukaryote

    --help                                                    This message

END
}
