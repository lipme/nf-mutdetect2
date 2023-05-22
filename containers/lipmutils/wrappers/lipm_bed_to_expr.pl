#!/usr/bin/perl

# $Id: lipm_bed_to_expr.pl 1096 2013-04-11 08:39:45Z carrere $

# Copyright INRA

# Jerome.Gouzy@toulouse.inra.fr
# Erika.Sallet@toulouse.inra.fr
# Sebastien.Carrere@toulouse.inra.fr
# Emmanuel.Courcelle@toulouse.inra.fr

# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".

# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.

# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.

# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.

BEGIN
{
    $ENV{LANG} = "C";
    our $VERSION = do {my @r = (q$Rev: 1096 $ =~ /\d+/g); $r[0]};
}

use strict;
use warnings;
use IO::File;
use ParamParser;
use Carp;
use General;
use GeneralBioinfo;
use RnaSeq;

my $DEFAULT_TYPES =
  "mRNA gene region repeat_region tRNA rRNA tmRNA insertion_sequence ncRNA transcript_region primary_transcript";
my $DEFAULTBEDSEARCH = " -s -f 1 ";    # one pair on the same strand
my $SCALINGFACTOR    = 1;
my $SCOREMAX         = &MAXINT;
my $MIN_EXON_LEN     = 18;

MAIN:
{
    my $o_param = New ParamParser(
                                  'GETOPTLONG',         \&Usage,
                                  "gff3=s",             "bedfile=s",
                                  "bedpe",              "outprefix=s",
                                  "option=s",           "num_of_mapped=i",
                                  "mode=s",             "scalingfactor=i",
                                  "scoremax=i",         "types=s",
                                  "tag=s",              "min_mean_coverage=i",
                                  "min_coverage=s",     "min_exon_len=i",
                                  "interval",           "verbose",
                                  'data_series_desc=s', 'data_series_name=s',
                                  'read_counts=s',      'covmax=i',
                                  'min_coverage_percent=i'
                                  );

    $o_param->AssertAllowedValue('mode', 'rpkm', 'depth', 'annotation', 'wig');
    $o_param->SetUnlessDefined('tag',    "SRCTAG");
    $o_param->SetUnlessDefined('types',  $DEFAULT_TYPES);

    #SC:20130411 Par defaut, on ne trimme pas (mail JG 20130411)
    #$o_param->SetUnlessDefined('covmax',$COVMAX);
    $o_param->SetUnlessDefined('option', $DEFAULTBEDSEARCH);

    my $mode = $o_param->Get('mode');
    if ($mode =~ /rpkm/)
    {
        #pairToBed and intersectBed in the path
        &Rpkm($o_param);
    }
    elsif ($mode eq 'depth' || $mode eq 'annotation')
    {
        &Depth($o_param, $mode);
    }
    elsif ($mode eq 'wig')
    {
        &Wig($o_param, $mode);
    }
}

sub LoadCov
{
    my ($o_param, $r_seqname, $ra_m, $ra_p, $ra_n, $ra_non_ig) = @_;

    $o_param->AssertDefined('bedfile');
    my $f_map     = $o_param->GetStreamIn('bedfile');
    my $countfile = $o_param->GetIfDefined('read_counts');
    my %h_counts  = ();
    if (defined($countfile))
    {
        my $f_counts = $o_param->GetStreamIn('read_counts');
        while (my $l = <$f_counts>)
        {
            chomp($l);
            my ($readid, $count) = split("\t", $l);
            $h_counts{$readid} = $count;
        }
        $f_counts->close;
    }
    my $count = 1;
    if ($o_param->IsDefined('bedpe'))
    {
        my $interval = ($o_param->IsDefined('interval')) ? 1 : 0;
        while (my $l = <$f_map>)
        {
            chomp($l);

#bedpe
#SMc     1819517 1819548 SMc     1819519 1819548 HWIEAS269_0008:6:100:10028:17684#CAGATC 62      +       +       reg=1819518-1819549 len=32
# in bed, start is 0-based and end is 1-based!
            my ($id, $min1, $max1, $t1, $min2, $max2, $readid, $t3, $strand1, $strand2) = split("\t", $l);

            # je me recale sur du 1-based
            $min1++;
            $min2++;
            if (defined($countfile))
            {
                next if (!defined($h_counts{$readid}));
                $count = $h_counts{$readid};
            }
            if (!defined($$r_seqname))
            {
                $$r_seqname = $id;
            }
            elsif ($id ne $$r_seqname)
            {
                Carp::croak("Fatal error: this program can't work on multi sequence files\n");
            }
            my $minmin = &Min($min1, $min2);
            my $maxmax = &Max($max1, $max2);
            &Add($ra_non_ig, $minmin, $maxmax, $count);
            if ($interval)
            {
                if ($strand1 eq $strand2)    # oriented
                {
                    my $ra_ = ($strand1 eq '+') ? $ra_p : $ra_m;
                    print STDERR "$strand1 $minmin,$maxmax\n" if ($o_param->IsDefined('verbose'));
                    &Add($ra_, $minmin, $maxmax, $count);
                }
                else
                {
                    &Add($ra_n, $minmin, $maxmax, $count);
                }
            }
            else
            {
                if ($strand1 eq $strand2)    # oriented
                {
                    my $ra_ = ($strand1 eq '+') ? $ra_p : $ra_m;
                    print STDERR "$strand1 $min1,$max1  $min2,$max2\n" if ($o_param->IsDefined('verbose'));
                    &Add($ra_, $min1, $max1, $count);
                    &Add($ra_, $min2, $max2, $count);
                }
                else
                {
                    &Add($ra_n, $min1, $max1, $count);
                    &Add($ra_n, $min2, $max2, $count);
                }
            }
        }
    }
    else
    {
        while (my $l = <$f_map>)
        {
            chomp($l);

            #bed >= 6
            #SMc     1819517 1819548 HWIEAS269_0008:6:100:10028:17684#CAGATC/1    62   +
            my ($id, $min, $max, $target, $score, $strand) = split("\t", $l);
            if (defined($countfile))
            {
                next if (!defined($h_counts{$target}));
                $count = $h_counts{$target};
            }

            # je me recale sur du 1-based
            $min++;
            if (!defined($$r_seqname))
            {
                $$r_seqname = $id;
            }
            elsif ($id ne $$r_seqname)
            {
                Carp::croak("Fatal error: this program can't work on multi sequence files\n");
            }
            if (defined($strand) && $strand =~ /[\+\-]/)
            {
                my $ra_ = ($strand eq '+') ? $ra_p : $ra_m;
                &Add($ra_, $min, $max, $count);
            }
            else
            {
                &Add($ra_n, $min, $max, $count);
            }
        }
    }
    $f_map->close;
}

# genere les fichiers wig pour gbrowse et apollo

sub Wig
{
    my ($o_param, $mode) = @_;

    my @a_m      = ();
    my @a_p      = ();
    my @a_n      = ();
    my @a_non_ig = ();
    my $seqname  = undef;
    &LoadCov($o_param, \$seqname, \@a_m, \@a_p, \@a_n, \@a_non_ig);

    $o_param->SetUnlessDefined('min_coverage', 1);
    my $min_coverage = $o_param->Get('min_coverage');
    my $min_exon_len = 1;
    my $name         = $o_param->Get('data_series_name');

    if (scalar(@a_p) > 0 && scalar(@a_m) > 0)
    {
        my $desc = $o_param->Get('data_series_desc');

     #&DumpRegionsWIG($seqname,$f_out,$name,$desc,$min_coverage,$min_exon_len,$o_param->Get('covmax'),'+-',\@a_p,\@a_m);
        foreach my $strand ('+', '-')
        {
            my $f_out = &GetStreamOut($o_param->Get('outprefix') . ".gbrowse.$strand.wig");
            my $ra_ = ($strand eq '+') ? \@a_p : \@a_m;
            &DumpRegionsWIG($seqname, $f_out, $name, $desc, $min_coverage, $min_exon_len, $o_param->Get('covmax'),
                            $strand, $ra_);
            $f_out->close;
        }
    }
    else
    {
        if (scalar(@a_p) > 0 && scalar(@a_m) == 0)
        {
            my $desc  = $o_param->Get('data_series_desc');
            my $f_out = &GetStreamOut($o_param->Get('outprefix') . ".gbrowse.+.wig");
            $desc .= " strand+";
            &DumpRegionsWIG($seqname, $f_out, $name, $desc, $min_coverage, $min_exon_len, $o_param->Get('covmax'),
                            '+', \@a_p);
            $f_out->close;
        }
        if (scalar(@a_m) > 0 && scalar(@a_p) == 0)
        {
            my $desc  = $o_param->Get('data_series_desc');
            my $f_out = &GetStreamOut($o_param->Get('outprefix') . ".gbrowse.-.wig");
            $desc .= " strand-";
            &DumpRegionsWIG($seqname, $f_out, $name, $desc, $min_coverage, $min_exon_len, $o_param->Get('covmax'),
                            '-', \@a_m);
            $f_out->close;
        }
    }
    if (scalar(@a_n) > 0)
    {
        my $desc  = $o_param->Get('data_series_desc');
        my $f_out = &GetStreamOut($o_param->Get('outprefix') . ".gbrowse.nostrand.wig");
        $desc .= " strand.";
        &DumpRegionsWIG($seqname, $f_out, $name, $desc, $min_coverage, $min_exon_len, $o_param->Get('covmax'),
                        '.', \@a_n);
        $f_out->close;
    }
    if (0)    # n'apporte rien par rapport a --interval
    {
        if (scalar(@a_non_ig) > 0 && $o_param->IsDefined('bedpe'))
        {
            my $desc = $o_param->Get('data_series_desc');
            $desc .= " transcribed_region";
            my $f_out = &GetStreamOut($o_param->Get('outprefix') . ".gbrowse.noig.wig");
            &DumpRegionsWIG($seqname, $f_out, $name, $desc, $min_coverage, $min_exon_len, $o_param->Get('covmax'),
                            '+', \@a_non_ig);
            $f_out->close;
        }
    }
}

sub Depth
{
    my ($o_param, $mode) = @_;

    my @a_m      = ();
    my @a_p      = ();
    my @a_n      = ();
    my @a_non_ig = ();
    my $seqname  = undef;
    &LoadCov($o_param, \$seqname, \@a_m, \@a_p, \@a_n, \@a_non_ig);

    $o_param->SetUnlessDefined('scalingfactor', $SCALINGFACTOR);
    $o_param->SetUnlessDefined('scoremax',      &MAXINT);
    my $scalef   = $o_param->Get('scalingfactor');
    my $maxvalue = $o_param->Get('scoremax');
    my $maxi     = &Max(&Max($#a_p, $#a_m), $#a_n);
    my $cpt      = 1;
    my $tag      = $o_param->Get('tag');

    my $so    = "transcript_region";
    my $f_out = &GetStreamOut($o_param->Get('outprefix') . ".position.gff3");
    $o_param->SetUnlessDefined('min_coverage', 1);
    my $min_coverage = $o_param->Get('min_coverage');

	if ($o_param->IsDefined('min_coverage_percent'))
	{	
		my @a_definedseqdata;
		for (my $i = 1; $i <= $maxi; $i++)
		{
			push(@a_definedseqdata, $a_m[$i]) if (defined($a_m[$i]));
			push(@a_definedseqdata, $a_p[$i]) if (defined($a_p[$i]));
		}
		#print "Compute mincoverage " . $o_param->Get('min_coverage_percent') . "% for " . $o_param->Get('outprefix') . ".position.gff3.";
		$min_coverage = &RnaSeq::GetMinCoverage(\@a_definedseqdata, $o_param->Get('min_coverage_percent'));
	}

    for (my $i = 1; $i <= $maxi; $i++)
    {
        if (defined($a_p[$i]) && $a_p[$i] >= $min_coverage)
        {
            &PrintGFFLine($seqname, $so, $f_out, 'depth', $tag, $i, $i, &Min($a_p[$i] / $scalef, $maxvalue),
                          '+', $cpt++);
        }
        if (defined($a_m[$i]) && $a_m[$i] >= $min_coverage)
        {
            &PrintGFFLine($seqname, $so, $f_out, 'depth', $tag, $i, $i, &Min($a_m[$i] / $scalef, $maxvalue),
                          '-', $cpt++);
        }
        if (defined($a_n[$i]) && $a_n[$i] >= $min_coverage)
        {
            &PrintGFFLine($seqname, $so, $f_out, 'depth', $tag, $i, $i, &Min($a_n[$i] / $scalef, $maxvalue) / 2,
                          '+', $cpt++);
            &PrintGFFLine($seqname, $so, $f_out, 'depth', $tag, $i, $i, &Min($a_n[$i] / $scalef, $maxvalue) / 2,
                          '-', $cpt++);
        }
    }
    $f_out->close;

    $o_param->SetUnlessDefined('min_exon_len', $MIN_EXON_LEN);
    $o_param->SetUnlessDefined('min_mean_coverage', ($min_coverage + $min_coverage / 10))
      ;    # => at least two overlapping hits if min_coverage == 1
    my $min_mean_coverage = $o_param->Get('min_mean_coverage');
    my $min_exon_len      = $o_param->Get('min_exon_len');

    $f_out = &GetStreamOut($o_param->Get('outprefix') . ".region.gff3");
    &DumpRegionsGFF3($seqname, $so, $f_out, 'region', $tag, $min_coverage, $min_mean_coverage, $min_exon_len, $scalef,
                     $maxvalue, '+', \@a_p, \@a_non_ig)
      if (scalar(@a_p) > 0);
    &DumpRegionsGFF3($seqname, $so, $f_out, 'region', $tag, $min_coverage, $min_mean_coverage, $min_exon_len, $scalef,
                     $maxvalue, '-', \@a_m, \@a_non_ig)
      if (scalar(@a_m) > 0);
    &DumpRegionsGFF3($seqname, $so, $f_out, 'region', $tag, $min_coverage, $min_mean_coverage, $min_exon_len, $scalef,
                     $maxvalue, '.', \@a_n, \@a_non_ig)
      if (scalar(@a_n) > 0);
    $f_out->close;

    if ($o_param->IsDefined('bedpe'))
    {
        $f_out = &GetStreamOut($o_param->Get('outprefix') . ".non_intergenic.gff3");
        &DumpRegionsGFF3($seqname, 'intergenic_region', $f_out, 'region', $tag, $min_coverage, $min_mean_coverage,
                         $min_exon_len, $scalef, $maxvalue, '.', \@a_non_ig, \())
          if (scalar(@a_non_ig) > 0);
        $f_out->close;
    }
}

sub DumpRegionsWIG
{
    my ($seqname, $f_out, $name, $desc, $min_coverage, $min_exon_len, $covmax, $strand, $ra_seq_cov_s1, $ra_seq_cov_s2)
      = @_;

    my @a_strand     = split('', $strand);
    my $cpt          = 1;
    my $region_start = 0;
    my $prev_cov     = -1;

    # for ucsc_wigToBigWig print $f_out "track type=wiggle_0 name=\"$name\" description=\"$desc\"\n";
    my $maxi = $#$ra_seq_cov_s1;
    if (defined($ra_seq_cov_s2))
    {
        $maxi = &Max($maxi, $#$ra_seq_cov_s2);
    }
    my $nb_strand = $#a_strand + 1;
    for (my $i = 1; $i <= $maxi; $i++)
    {
        if (defined($a_strand[0]))
        {
            if (defined($$ra_seq_cov_s1[$i]) && $$ra_seq_cov_s1[$i] >= $min_coverage)
            {
                $$ra_seq_cov_s1[$i] = $$ra_seq_cov_s1[$i] > $covmax ? $covmax : $$ra_seq_cov_s1[$i]
                  if (&IsNotNull($covmax) == &TRUE);
                my $score =
                  ($a_strand[0] eq '+' || $a_strand[0] eq '.') ? $$ra_seq_cov_s1[$i] : (-1 * $$ra_seq_cov_s1[$i]);
                printf $f_out "$seqname\t%d\t%d\t$score\n", $i - 1, $i;
            }
        }
        if (defined($a_strand[1]))
        {
            if (defined($$ra_seq_cov_s2[$i]) && $$ra_seq_cov_s2[$i] >= $min_coverage)
            {
                $$ra_seq_cov_s2[$i] = $$ra_seq_cov_s2[$i] > $covmax ? $covmax : $$ra_seq_cov_s2[$i]
                  if (&IsNotNull($covmax) == &TRUE);
                my $score =
                  ($a_strand[1] eq '+' || $a_strand[1] eq '.') ? $$ra_seq_cov_s2[$i] : (-1 * $$ra_seq_cov_s2[$i]);
                printf $f_out "$seqname\t%d\t%d\t$score\n", $i - 1, $i;
            }
        }
    }
}

sub DumpRegionsGFF3
{
    my (
        $seqname,      $so,                $f_out,        $type,   $tag,
        $min_coverage, $min_mean_coverage, $min_exon_len, $scalef, $scoremax,
        $strand,       $ra_seq_cov,        $ra_seq_non_ig
        ) = @_;

    my $cpt          = 1;
    my $region_start = 0;
    my $mean_cov     = 0;
    for (my $i = 1; $i <= ($#$ra_seq_cov + 1); $i++) # to be sure to test the last undef and then print out the last hit
    {
        if ($region_start)
        {
            if (!defined($$ra_seq_cov[$i]) || $$ra_seq_cov[$i] < $min_coverage)
            {
                my $len = ($i - 1) - $region_start + 1;    # $i-1 is the position of the last match
                                                           #my $mean_cov_raw = $mean_cov;
                $mean_cov = $mean_cov / $len;
                if ($len >= $min_exon_len && $mean_cov >= $min_mean_coverage)
                {
                    my $score = &Min($mean_cov / $scalef, $scoremax);
                    if ($so eq 'intergenic_region')
                    {
                        $score *= -1;
                    }
                    if ($strand eq '.')
                    {
                        &PrintGFFLine($seqname, $so,        $f_out, $type, $tag, $region_start,
                                      $i - 1,   $score / 2, '+',    $cpt++);
                        &PrintGFFLine($seqname, $so,        $f_out, $type, $tag, $region_start,
                                      $i - 1,   $score / 2, '-',    $cpt++);
                    }
                    else
                    {
                        &PrintGFFLine($seqname,      $so,    $f_out, $type,   $tag,
                                      $region_start, $i - 1, $score, $strand, $cpt++);
                    }

#printf $f_glint "$seqname\tglint\t$taggff\t$region_start\t$i\t%.1f\t.\t.\tID=glint$tagid.$cpt;length=$len;\n",$mean_cov;
#printf $f_glint "$seqname\tglint\t$taggff\t$region_start\t$i\t%.1f\t.\t.\tID=glint$tagid.$cpt.1;Parent=glint$tagid.$cpt;Target=target$tagid.$cpt 1 $len +;length=$len;\n",$mean_cov;
                }
                $cpt++;
                $region_start = 0;
                next;
            }
        }
        elsif (defined($$ra_seq_cov[$i]) && $$ra_seq_cov[$i] ne '' && $$ra_seq_cov[$i] >= $min_coverage)
        {
            $region_start = $i;
            $mean_cov     = 0;
        }
        if (defined($$ra_seq_cov[$i]) && $$ra_seq_cov[$i] ne '')
        {
            $mean_cov += $$ra_seq_cov[$i];
            $$ra_seq_non_ig[$i] = 0 if (defined($$ra_seq_non_ig[$i]));
        }
    }
}

sub PrintGFFLine
{
    my ($seqname, $so, $f_out, $type, $tag, $start, $end, $score, $strand, $cpt) = @_;

    my ($code) = ($strand eq '.') ? "u" : (($strand eq '+') ? 'p' : 'm');
    printf $f_out "$seqname\t$type$tag\t$so\t$start\t$end\t%.4f\t$strand\t.\tID=$seqname:$type$tag:$code$start.$cpt\n",
      $score;
}

sub Rpkm
{
    my ($o_param) = @_;

    $o_param->AssertFileExists('gff3', 'bedfile');
    $o_param->AssertInteger('num_of_mapped');
    my $gff3       = $o_param->Get('gff3');
    my %h_item_len = ();
    &GetInfo($gff3, \%h_item_len);
    my $patt_types = $o_param->Get('types');
    $patt_types .= " exon" if ($patt_types =~ /mRNA/);
    $patt_types =~ s/^\s+//;
    $patt_types =~ s/\s+$//;
    $patt_types =~ s/[ ,]+/\|/g;
    my $option    = $o_param->Get('option');                    # pairToBed/intersectBed options
    my $bedfile   = $o_param->Get('bedfile');
    my $nb_mapped = $o_param->Get('num_of_mapped') / 1000000;
    my $cde       = "";
    my $colid     = 0;
    my $compress  = &GetCompressPath($bedfile);

    if ($o_param->IsDefined('bedpe'))
    {
        my $exe_path = &GetExePath("pairToBed");
        if ($compress eq '')
        {
            $cde = "$exe_path $option -a $bedfile -b $gff3 |";    #| sort -u |";
        }
        else
        {
            $cde = "$compress -cd $bedfile | $exe_path $option -a stdin -b $gff3 |";    #| sort -u |";
        }

        # je ne mets pas le sort -u. La meme paire va donc compter pour 2 si elle mappe le meme exon
        # ce qui sera coherent avec le cas ou les deux extremites mappent avec des exons differents.
        # par contre si je mets le sort -u cela compterait 1 dans le premier cas et 2 dans le second
        # 20110712
        # rajout d'un hash pour qu'une paire ne compte que pour 1 pour chaque item
        # pour analyser les variants d'epissage il faut travailler au niveau exon
        # car desormais les PE sur plusieurs exons soumis a variation ne permettent plus de mettre
        # autant de "contraste" entre les niveaux d'expression des variants
        $colid = 6;
    }
    else
    {
        my $exe_path = &GetExePath("intersectBed");
        if ($compress eq '')
        {
            $cde = "$exe_path -wo $option -a $bedfile -b $gff3 |";    #| sort -u |";
        }
        else
        {
            $cde = "$compress -cd $bedfile | $exe_path -wo $option -a stdin -b $gff3 |";    #| sort -u |";
        }
        $colid = 3;
    }
    my $f_in          = new IO::File("$cde") or Carp::croak("Can't open >$cde<\n");
    my %h_tag_count   = ();
    my %h_tag_start   = ();
    my %h_tag_end     = ();
    my %h_tag_strand  = ();
    my %h_tag_seqname = ();
    my %h_tag_lab     = ();

#bedpe
#SMc     84432   84473   SMc     84436   84473   HWUSI-EAS691_0027:1:100:10030:20725#ATCACG      78      +       +       reg=84433-84473 len=41  SMc     iANT region   81374   87895   .       +       .       ID=region:SMc02938_MISC.92;Name=rrn,SMc02938;Note=
#bed

    my %h_pair_id = ();    # a pair must count one per biological item
    my @a_        = ();
    while (my $l = <$f_in>)
    {
        next if ($l =~ /^#/);
        chomp($l);
        $l =~ s/\s+$//;
        @a_ = split("\t", $l);
        if ($l =~ /($patt_types).+ID=([^\;\s]+)/)
        {
            my $type = $1;
            my $id   = $2;
            if (   $type eq 'exon'
                && $l =~ /Parent=([^\;\s]+)/)    # case: mRNA  count at the exon level for mRNA allows variant analysis
            {
                my $parent = $1;
                my $pair   = "$parent-$a_[$colid]";
                if (!defined($h_pair_id{$pair}))
                {
                    $h_tag_count{$parent}++;
                }
                $h_pair_id{$pair} = 1;
            }
            if ($type ne 'mRNA')
            {
                my $pair = "$id-$a_[$colid]";
                if (!defined($h_pair_id{$pair}))
                {
                    $h_tag_count{$id}++;
                }
                $h_pair_id{$pair} = 1;
            }
        }
    }
    $f_in->close;

    my $f_gff3 = $o_param->GetStreamIn('gff3');
    my $f_out  = $o_param->GetStreamOut('outprefix');
    my $cpt    = 0;
    while (my $line = <$f_gff3>)
    {
        chomp($line);
        $line =~ s/\s+$//;
        if ($line =~ /($patt_types).+ID=([^\;\s]+)/)
        {
            my $type = $1;
            my $id   = $2;
            $h_tag_count{$id} = 0 if (!defined($h_tag_count{$id}));
            my $rpkm = ($h_tag_count{$id}) / ($h_item_len{$id} / 1000);
            $rpkm /= ($nb_mapped);    # deja / 1000000 plus haut
            $line .= sprintf(";rpkm=%.2f;count=$h_tag_count{$id}", $rpkm);
            $line .= ";length=$h_item_len{$id}" if ($line !~ /length=/);
        }
        print $f_out "$line\n";
    }
    $f_out->close;
    $f_gff3->close;
}

sub GetInfo
{
    my ($gff3, $rh_item_len) = @_;

    my $f_gff3 = &GetStreamIn($gff3);
    my %h_line = ();
    my $nbexon = 0;
    my @a_     = ();
    while (my $line = <$f_gff3>)
    {
        chomp($line);
        @a_ = split(' ', $line);
        &ParseGFF3Line(\$line, 1, \%h_line);
        next if (!defined($h_line{'type'}));
        my $type = $h_line{'type'};
        next
          if ($type eq 'mRNA')
          ;    # la longueur doit etre mise a jour en concatenant les exons (pour pouvoir analyser les variants)
        if ($type eq 'exon')
        {
            $$rh_item_len{$h_line{'Parent'}} += $h_line{'len'};
            $nbexon++;
        }
        $$rh_item_len{$h_line{'ID'}} = $h_line{'len'};
    }
    $f_gff3->close;
}

sub Add
{
    my ($ra_, $min, $max, $count) = @_;

    for (my $i = $min; $i <= $max; $i++)
    {
        $$ra_[$i] += $count;
    }
}

sub Usage
{
    print STDERR<<END
$0 
	--help
	--verbose
	--bedfile filename
	--outprefix filename
	--mode    rpkm|depth|annotation|wig
	--bedpe                                     flag to specify paired-end input file

All modes
	--read_counts filename                      allows the rescue of expression levels when reads/pairs have been non-redundified before mapping (e.g sRNA)
                                                or the selection of a subset of the mapped reads

Rpkm mode
	--gff3    filename
	--option "pairToBed/intersectBed options" e [$DEFAULTBEDSEARCH] (e.g -s to enforce strandedness when finding overlaps.)
	--num_of_mapped integer                     Number of ends/pair-ends used for the normalisation (e.g # that can be mapped unambiously)
	--types  "quoted string of GFF3 types"      [$DEFAULT_TYPES]

Depth/annotation mode
	--scalingfactor   integer [$SCALINGFACTOR]
	--scoremax        integer [$SCOREMAX]	    (e.g set to $SCOREMAX, all scores > $SCOREMAX (after rescaling) will be set to 1)
	--tag
	--min_mean_coverage
	--min_coverage    integer	[1]
    --min_coverage_percent integer Minimum percentage of reads to keep [If defined, disable min_coverage parameter]
	--min_exon_len	  integer	[$MIN_EXON_LEN]
	--interval                                  this flag manages bedpe file as a single interval

Wig mode
	--tag
	--min_coverage    integer	[1]
	--interval                               this flag manages bedpe file as a single interval
	--covmax          integer  	    (e.g set to covmax, all scores > covmax (after rescaling) will be set to 1)
	--data_series_name string  []	        name of the data serie (e.g: "ArrayExpt1" )
	--data_series_desc string  []         description of the data serie (e.g: "20 degrees, 2 hr")

Warnings:

	when --bedpe is activated, the hit is considered as oriented when both ends are on the same strand
	when the strand is unknown, two hits on both strand are reported on the same region (half of the score for each strand)
	
END
}
