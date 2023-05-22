#!/usr/bin/perl

#
# $Id: lipm_fastaq2fastaq.pl 752 2011-11-22 12:23:05Z gouzy $
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

lipm_fastaq2fastaq.pl

=head1 SYNOPSIS

Convert a fastq file in fasta format

lipm_fastaq2fastaq.pl --in file.fastq.xz --out file.fa.xz

+ extract a substring from position 5 to 65

	lipm_fastaq2fastaq.pl --in file.fastq.xz --out file.fa.xz --start 5 --end 65

+ extract a substring of length 61 starting at position 5

	lipm_fastaq2fastaq.pl --in file.fastq.xz --out file.fa.xz --start 5 --length 61

+ remove all reads with 'N'  (considering the substring when start/end/length are specified)

	lipm_fastaq2fastaq.pl --in file.fastq.xz --out file.fa.xz --qualfilter remove

+ trim read sequences after the first N  (considering the substring when start/end/length are specified)

	lipm_fastaq2fastaq.pl --in file.fastq.xz --out file.fa.xz --qualfilter trim

+ mask (with a N) read sequences after the first N (considering the substring when start/end/length are specified)

	lipm_fastaq2fastaq.pl --in file.fastq.xz --out file.fa.xz     --qualfilter mask

=head1 DESCRIPTION

Convert/edit fasta/fastq read files

=cut

BEGIN
{
	use File::Basename;
	my $dirprg=dirname($0);
	unshift @INC, "$dirprg/../corelib";
	$VERSION = do { my @r = (q$Rev: 752 $ =~ /\d+/g); sprintf "%d."."%02d" x $#r, @r }
}

use strict;
use warnings;
use Carp;
use ParamParser;
use General;
use GeneralBioinfo;

MAIN:
{
    my $o_param = New ParamParser('GETOPTLONG', \&Usage, "in=s", "out=s", "help", "start=i", "end=i", "length=i", "outfastaq=s", "minlen=i","qualfilter=s");

	$o_param->AssertDefined('in', 'out'); 
	$o_param->AssertAllowedValue('qualfilter','trim','remove','mask','insert') if ( $o_param->IsDefined('qualfilter'));
	
	my $fh_in      = $o_param->GetStreamIn('in');
	my $fh_out     = $o_param->GetStreamOut('out');

	my $length   = \undef; # \ is required for running substr
	if ( $o_param->IsDefined('start') && $o_param->IsDefined('end') )
	{
		$length = $o_param->Get('end')-$o_param->Get('start')+1;
	}
	elsif ( $o_param->IsDefined('end') )
	{
		$o_param->SetUnlessDefined('length',$o_param->Get('end'));	
		$length = $o_param->Get('end');
	}
	elsif ($o_param->IsDefined('length'))
	{
		$length = $o_param->Get('length');
	}
	
	$o_param->SetUnlessDefined('start',1);	
	$o_param->SetUnlessDefined('end',&MAXINT);	
	$o_param->SetUnlessDefined('minlen',0);	

	my $start = $o_param->Get('start');
	my $end   = $o_param->Get('end');

	Carp::croak("Error: start value ($start) must be lower than end value ($end)\n") if ( $start > $end );

	my $fh_outqual = undef;
	if ( $o_param->IsDefined('outfastaq') ) 
	{
		$fh_outqual = $o_param->GetStreamOut('outfastaq');
	}

	my $qualfilter   = $o_param->Get('qualfilter');
	my $minlen       = $o_param->Get('minlen');

	my ($id,$seq,$qual)  = ('','','');
	while(my $l=<$fh_in>)
	{
		chomp($l);
		if ( $l =~ /^\@(\S+).*/ )
		{
			$id = $1;
			$seq  = <$fh_in>;
			my $tmp  = <$fh_in>;
			$qual = <$fh_in>;
			chomp($qual);	
			$qual = substr($qual,$start-1,$length) if (defined($fh_outqual));
		}
		elsif ( $l =~ /^>(\S+)/ )
		{
			$id = $1;
			$seq  = <$fh_in>;
		}
		chomp($seq);
		$seq  = substr($seq,$start-1,$length);
		if ( $qualfilter )
		{
			if ( $qualfilter eq 'remove' )
			{
				next if ($seq =~ /N/i );
				next if ($qual=~ /B/  );
			}
			elsif ( $qualfilter eq 'mask' )
			{
				$seq =~ s/(N.*)/'N' x length($1)/e;
			}
			elsif ( $qualfilter eq 'insert' )
			{
				my $insert = 'N';
				while( $seq =~ /([^N]+)/g)
				{
					$insert = $1 if ( length($1) > length($insert) );
				}
				$seq = $insert;
			}
			elsif ( $qualfilter eq 'trim' )
			{
				$seq =~ s/N.*//i;
				$seq = 'N' if ($seq eq '');
				$qual = substr($qual,0,length($seq)) if ( defined($fh_outqual) );
			}
		}
		my $len = length($seq);
		next if ( $len < $minlen );

		print $fh_out ">$id len=$len\n$seq\n";


		next if ( ! defined($fh_outqual) );
		print $fh_outqual ">$id len=$len\n";
		my @a_ = split('',$qual);
		for(my $i=0;$i<=$#a_;$i++)
		{
			my $sq = $a_[$i];
			my $Q = 10 * log(1 + 10 ** (ord($sq) - 64) / 10.0) / log(10);
			$Q = 255 if ( $Q > 255 );
			printf $fh_outqual "%d",int($Q);
			print $fh_outqual " " if ( $i != $#a_ );
		}
		print $fh_outqual "\n";
	}

	$fh_in->close();
	$fh_out->close();
	$fh_outqual->close() if ( defined($fh_outqual) );
}

sub Usage
{
    print STDERR<<END

Convert/edit fasta/fastq read files

Usage:
$0	
    --help
Mandatory
    --in         filename|stdin	name of a fasta or fastq stream
    --out        filename|stdout

	compressed filenames are allowed

Optional
    --start      integer
    --end        integer
    --length     integer
    --outfastaq  filename|stdout  report quality values in a fastaq formatted file
    --qualfilter  string       accepted values are 
                                 remove : remove reads containing 'N'
                                 mask   : replace by a 'N' all nt found after a 'N'
                                 trim   : remove all nt found after a 'N'
                                 insert : identify the longest insert sequence (i.e longest unmasked sequence)
    --minlen  integer          if defined, reports only reads longer than minlen          

WARNING 1 : one-line sequence only
WARNING 2 : parameter "minlen" and "qualfilter remove" change the sequence number in the file 	
	=> the fasta files cannot be used (directly) for assembling paired-end
	=> some analytic software expect the same length for all reads

END
}

