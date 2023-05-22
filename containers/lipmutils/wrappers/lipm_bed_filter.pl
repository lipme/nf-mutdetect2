#!/usr/bin/perl

#
# $Id: lipm_bed_filter.pl 1236 2013-08-23 15:11:17Z gouzy $
#

# Copyright INRA

# Erika.Sallet@toulouse.inra.fr
# Jerome.Gouzy@toulouse.inra.fr
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

=pod 

=head1 NAME


=head1 SYNOPSIS


=head1 DESCRIPTION


=cut

BEGIN
{
    $ENV{LANG} = "C";
    $VERSION = do {my @r = (q$Rev: 1236 $ =~ /\d+/g); sprintf "%d." . "%02d" x $#r, @r}
}

use strict;
use warnings;
use ParamParser;
use General;
use GeneralBioinfo;
my $SEP     = "£";
my $PE_TAGS = "/1|/2";

MAIN:
{
    my $o_param = New ParamParser(
                                  'GETOPTLONG',         \&Usage,         'in=s',            'out=s',
                                  'bedpe',              'score=i',       'insert_length=i', 'start=i',
                                  'end=i',              'specific',      'best',            'chrom=s',
                                  'pre_sorted_by_name', 'db_filename=s', 'db_column=s',     'db_insert',
                                  'keep',               'remove',        'help', 'nolog',   'oriented_end=i'
                                  );
    $o_param->AssertDefined('in', 'out');
    my $bedfmt = 1;
    my @a_sort_fields = ("4", "5nr");
    if ($o_param->IsDefined('bedpe'))
    {
        $bedfmt = 2;
        @a_sort_fields = ("7", "8nr");
    }
    my $insertid = ($o_param->IsDefined('db_insert')) ? &TRUE : &FALSE;
    my $specific = ($o_param->IsDefined('specific'))  ? &TRUE : &FALSE;
    my $best     = ($o_param->IsDefined('best'))      ? &TRUE : &FALSE;
    my $remove   = ($o_param->IsDefined('remove'))    ? &TRUE : &FALSE;
    my $keep     = ($o_param->IsDefined('keep'))      ? &TRUE : &FALSE;
    $o_param->SetUnlessDefined('score',         -1 * &MAXINT);
    $o_param->SetUnlessDefined('insert_length', 0);
    $o_param->SetUnlessDefined('start',         -1);
    $o_param->SetUnlessDefined('end',           &MAXINT);
    my $score         = $o_param->Get('score');
    my $insert_length = $o_param->Get('insert_length');
    my $debut         = $o_param->Get('start');
    my $fin           = $o_param->Get('end');
    my $chrom         = $o_param->GetIfDefined('chrom');

    my %h_select = ();
    if ($o_param->IsDefined('db_filename'))
    {
        $o_param->SetUnlessDefined('db_column', 1);
        my $db_col = $o_param->Get('db_column');
        $keep = &TRUE if (!$remove);
        &TabToHash($o_param->Get('db_filename'), \%h_select, $db_col, $db_col);
    }

    my %h_best   = ();
    my $f_out    = $o_param->GetStreamOut('out');
    my $tempfile = "lipm_bed_filter.tmp.$$." . time;
    my $f_in     = undef;
    if ($o_param->IsDefined('pre_sorted_by_name'))
    {
        $f_in = $o_param->GetStreamIn('in');
    }
    else
    {
        my $f_inraw = $o_param->GetStreamIn('in');
        &UnixSort($f_inraw, $tempfile, \@a_sort_fields);
        $f_inraw->close;
        $f_in = &GetStreamIn($tempfile);
    }
    my %h_res     = ();
    my %h_resprev = ();
    my $prev_tid  = '';
    my $tid_hits  = '';
	my $nb_reads_or_pairs=0;
	my $nb_specific_reads_or_pairs=0;
	my $nb_total_loci=0;
	my $oriented_end=$o_param->GetIfDefined('oriented_end');
	if ( defined($oriented_end) )
	{
		$o_param->AssertAllowedValue('oriented_end',1,2);
	}
    while (my $l = <$f_in>)
    {
        chomp($l);
        &ParseBedLine(\$l, $bedfmt, \%h_res);
        if ($h_res{name} ne $prev_tid && $prev_tid ne '' )
        {
			$nb_reads_or_pairs++;
            &DumpBedHits($f_out, $bedfmt, \$tid_hits, $specific,\$nb_specific_reads_or_pairs,\$nb_total_loci);
            $tid_hits = '';
        }
        $prev_tid = $h_res{name};
        if ($insertid)
        {
            my $insert = $h_res{name};
            $insert =~ s/($PE_TAGS)$//;
            next if ($keep   && !defined($h_select{$insert}));
            next if ($remove && defined($h_select{$insert}));
        }
        else
        {
            next if ($keep   && !defined($h_select{$h_res{name}}));
            next if ($remove && defined($h_select{$h_res{name}}));
        }
        next if ($h_res{score} < $score);
        if ($bedfmt == 2)
        {
            next if (defined($chrom) && ($h_res{chrom1} ne $chrom || $h_res{chrom2} ne $chrom));
            my ($start) = &Min($h_res{start1}, $h_res{start2});
            my ($end) = &Max($h_res{end1}, $h_res{end2});
            next if (!($start >= $debut && $end <= $fin));
            if ($insert_length)
            {
                my $delta = abs($end - $start) + 1;
                next if ($delta < $insert_length);
            }
			if (defined($oriented_end))
			{
				($h_res{strand1},$h_res{strand2}) = ( $oriented_end == 1 ) ? ($h_res{strand1},$h_res{strand1}) : ($h_res{strand2},$h_res{strand2});
				$l = "$h_res{chrom1}\t".($h_res{start1}-1)."\t$h_res{end1}\t$h_res{chrom2}\t".($h_res{start2}-1)."\t$h_res{end2}\t$h_res{name}\t$h_res{score}\t$h_res{strand1}\t$h_res{strand2}";
				$l .= "\t$h_res{misc}" if (defined($h_res{misc}) && $h_res{misc} ne '');
			}
        }
        else
        {
            next if (defined($chrom) && $h_res{chrom} ne $chrom);
            next if (!($h_res{start} >= $debut && $h_res{end} <= $fin));

			if (defined($oriented_end) && $oriented_end == 2)
			{
				$h_res{strand} = ( $h_res{strand} eq '+' ) ? '-' : '+';
				$l = "$h_res{chrom}\t".($h_res{start}-1)."\t$h_res{end}\t$h_res{name}\t$h_res{score}\t$h_res{strand}";
				$l .= "\t$h_res{misc}" if (defined($h_res{misc}) && $h_res{misc} ne '');
			}
        }
        if (!$specific && !$best)
        {
            $tid_hits .= $SEP . $l;
            next;
        }
        if ($tid_hits ne '')
        {
            &ParseBedLine(\$tid_hits, $bedfmt, \%h_resprev);
            if ($h_res{score} > $h_resprev{score})
            {
                $tid_hits = $l;
            }
            elsif ($h_res{score} == $h_resprev{score})
            {
                $tid_hits .= $SEP . $l;
            }
        }
        else
        {
            $tid_hits = $l;
        }
    }
	$nb_reads_or_pairs++;
	&DumpBedHits($f_out, $bedfmt, \$tid_hits, $specific,\$nb_specific_reads_or_pairs,\$nb_total_loci);
    $f_in->close;
    $f_out->close;
    if (! $o_param->IsDefined('nolog'))
	{
    	my $f_log    = &GetStreamOut($o_param->Get('out').".log");
		print $f_log "nb_total_loci\t$nb_total_loci\n";
		print $f_log "nb_total_reads_or_pairs\t$nb_reads_or_pairs\n";
		print $f_log "nb_specific_reads_or_pairs\t$nb_specific_reads_or_pairs\n";
		$f_log->close;
	}
    unlink($tempfile);
    exit(0);
}

sub DumpBedHits
{
my ($f_out, $bedfmt, $r_tid_hits, $specific,$r_nb_specific_reads_or_pairs,$r_nb_total_loci) = @_;

    $$r_tid_hits =~ s/^$SEP//;
    my @a_hits = split("$SEP", $$r_tid_hits);
	$$r_nb_specific_reads_or_pairs++ if ( scalar(@a_hits) == 1);
	$$r_nb_total_loci+=scalar(@a_hits);
    if ($specific)
    {
        print $f_out "$a_hits[0]\n" if (scalar(@a_hits) == 1);
    }
    else
    {
        for (my $i = 0; $i < scalar(@a_hits); $i++)
        {
            print $f_out "$a_hits[$i]\n";
        }
    }
}

sub Usage
{
    print STDERR<<END
$0
	--help
	--in   filename|stdin		bed or bedpe format
	--out  filename|stdout		bed or bedpe format

	--bedpe 					flag to specify bedpe format

	--best						flag to keep only best hit(s)
	--specific 					flag to keep only unambiguous best hit(s)

	--nolog                     do not create the log file

Hit features
	--insert_length	integer 	extract insert longer than this value (bedpe only)
	--score         integer     score cutoff

Extracting a region
	--start integer
	--end   integer
	--chrom string

Selection from a tabulated file
	--db_filename filename 
	--db_column   integer       column corresponding to sequence IDs
	--db_insert                 use this flag to select read IDs belonging to the insert IDs in the db_filename
	--keep						flag to keep the selection
	--remove            		flag to remove the selection

Strand modification (may be required with single end reads depending on protocols)
	--oriented_end integer      [1|2] specifies which end gives the orientation (modifies strand field)

Warning: input data are sorted in the current directory, a lot of space can be required depending on the size of the input
	--pre_sorted_by_name        this flag avoids the sorting step when data are already sorted (use it with care)

END
}

