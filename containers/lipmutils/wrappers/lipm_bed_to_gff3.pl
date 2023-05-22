#!/usr/bin/perl

#
# $Id: lipm_bed_to_gff3.pl 1447 2015-06-26 14:51:53Z sallet $
#

# Copyright INRA

# Erika.Sallet@toulouse.inra.fr
# Jerome.Gouzy@toulouse.inra.fr
# Sebastien.Carrere@toulouse.inra.fr
# Emmanuel.Courcelle@toulouse.inra.fr

# This software is a computer program whose purpose is to provide a
# web-based system for detecting, clustering and annotating non protein
# coding RNA genes.

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

lipm_bed_to_gff3.pl

=head1 SYNOPSIS

lipm_bed_to_gff3.pl --in file.bed --out file.gff3 [--bedpe]

=head1 DESCRIPTION

Convert bed or bedpe file to gff3 

=cut

BEGIN
{
    $ENV{LANG} = "C";
    $VERSION = do {my @r = (q$Rev: 1447 $ =~ /\d+/g); sprintf "%d." . "%02d" x $#r, @r}
}

use strict;
use warnings;
use Carp;
use ParamParser;
use General;
use GeneralBioinfo;

use constant SOURCE => ".";
use constant CODE   => "transcript_region";

MAIN:
{

    # Get parameters
    my $o_param =
      New ParamParser('GETOPTLONG', \&Usage, "in=s", "out=s", "source=s", "code=s", "help", 'bedpe', 'oriented_end=s','one_level');

    $o_param->AssertDefined('in', 'out');

    my $fh_in  = $o_param->GetStreamIn('in');
    my $fh_out = $o_param->GetStreamOut('out');

    $o_param->SetUnlessDefined('source', &SOURCE);
    $o_param->SetUnlessDefined('code',   &CODE);
    my $source = $o_param->Get('source');    # Second column of gff3 output
    my $code   = $o_param->Get('code');      # Third column of gff3 output

    my $oriented_end = undef;
    if ($o_param->GetIfDefined('oriented_end'))
    {
        $o_param->AssertAllowedValue('oriented_end', 1, 2, '/1', '/2');
        $oriented_end = $o_param->Get('oriented_end');
        $oriented_end =~ s/\///;
    }

    my $cpt = 0;                             # Use to set ID to gff3 features
    my $id_prefix;
    my $name;

    my $bedfmt = ($o_param->IsDefined('bedpe')) ? 2 : 1;
    my $onelevel = ($o_param->IsDefined('one_level')) ? &TRUE : &FALSE;

    print $fh_out "##gff-version 3\n";

    my %h_hit = ();
    my @a_    = ();
    while (my $line = <$fh_in>)
    {
        next if ($line =~ /^(\#|\;)/);
        chomp $line;
        &ParseBedLine(\$line, $bedfmt, \%h_hit);

        #use Data::Dumper;
        #print Dumper(%h_hit);
        $cpt++;
        if ($bedfmt == 1)
        {
            if (defined($h_hit{name}) && $h_hit{name} ne "")
            {
                $id_prefix = $h_hit{chrom} . "." . $h_hit{name} . "." . $cpt;
                $name      = $h_hit{name};
            }
            else
            {
                $id_prefix = "$h_hit{chrom}.$$.$cpt";
                $name      = "$h_hit{chrom}-$h_hit{start}-$h_hit{end}";
            }
            $h_hit{score}  = 0   if (!defined($h_hit{score})  || $h_hit{score} eq "");
            $h_hit{strand} = '.' if (!defined($h_hit{strand}) || $h_hit{strand} eq "");
            print $fh_out
              "$h_hit{chrom}\t$source\t$code\t$h_hit{start}\t$h_hit{end}\t$h_hit{score}\t$h_hit{strand}\t.\tID=$code:$id_prefix;Name=$name\n";
        }
        elsif ($bedfmt == 2)
        {
            next if ($h_hit{chrom1} ne $h_hit{chrom2});    # unable to manage multiple sequence hits in GFF3
            my $id = "$h_hit{name}.$cpt";
            @a_ = ($h_hit{start1}, $h_hit{start2}, $h_hit{end1}, $h_hit{end2});
            my ($start, $end)  = &MinMaxArray(\@a_);
            my ($min1,  $max1) = &MinMax($h_hit{start1}, $h_hit{end1});
            my ($min2,  $max2) = &MinMax($h_hit{start2}, $h_hit{end2});

            # in GFF3, both ends must be on the same strand
            my $strand = '.';
			my $misc = '';
            if (defined($oriented_end))
            {
                $strand = ($oriented_end == 1) ? $h_hit{strand1} : $h_hit{strand2};
            }
			if ( defined($h_hit{misc}) )
			{	
				$misc = ";".$h_hit{misc};
				$misc =~ s/\s+/\;/g;
			}
            print $fh_out "$h_hit{chrom1}\t$source\tmatch\t$start\t$end\t.\t$strand\t.\tID=$id;Name=$h_hit{name}$misc\n";
			next if ( $onelevel );
            print $fh_out
              "$h_hit{chrom1}\t$source\tmatch_part\t$min1\t$max1\t.\t$strand\t.\tID=$id.1;Name=$h_hit{name}.1;Parent=$id\n";
            print $fh_out
              "$h_hit{chrom1}\t$source\tmatch_part\t$min2\t$max2\t.\t$strand\t.\tID=$id.2;Name=$h_hit{name}.2;Parent=$id\n";
        }
    }
    $fh_in->close();
    $fh_out->close();
}

=item Usage

	$o_param->Usage();
	Print the usage of the program

=cut

sub Usage
{
    my $code   = &CODE;
    my $source = &SOURCE;

    print STDERR<<END
Usage:
$0	

Convert bed/bedpe file into gff3

   --help                  this message

[Mandatory]
   --in      filename       bed file
   --out     filename       gff3 formatted output file

[Optional]
   --source  string      algorithm [$source]
   --code   string       region type [$code]

specific to bedpe format
   --bedpe               flag to activate bedpe mode, default is bed
   --oriented_end string specify which end gives the orientation (1, 2, /1, /2)
   --one_level           do not report the details of the hit (reports match butnot match_part)
END
}

