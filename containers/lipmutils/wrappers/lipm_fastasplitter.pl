#!/usr/bin/perl

#
# $Id$
#

# Copyright INRA

# Sebastien.Carrere@toulouse.inra.fr
# Erika.Sallet@toulouse.inra.fr
# Ludovic.Cottret@toulouse.inra.fr
# Ludovic.Legrand@toulouse.inra.fr
# Jerome.Gouzy@toulouse.inra.fr

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

lipm_fastasplitter.pl

=head1 SYNOPSIS

lipm_fastasplitter.pl --in fastafile --outprefix string [--num_per_slice integer] [--maxlength_per_slice integer]


=head1 DESCRIPTION


=cut


BEGIN
{

		$ENV{LANG} = "C";
        $VERSION = do { my @r = (q$Rev: 897 $ =~ /\d+/g); sprintf "%d."."%02d" x $#r, @r }
}

use strict;
use warnings;
use Carp;
use ParamParser;
use General;
use GeneralBioinfo;

my $NUM_PER_SLICE = 1000000;

MAIN:
{
	my $o_param = New ParamParser('GETOPTLONG',\&Usage,'help','in=s','outprefix=s','num_per_slice=i','maxlength_per_slice=i');
	$o_param->AssertDefined('in');
	$o_param->SetUnlessDefined('outprefix',$o_param->Get('in').".slice");
	$o_param->SetUnlessDefined('num_per_slice',$NUM_PER_SLICE);
	my $maxlen=0;
	if ( $o_param->IsDefined('maxlength_per_slice') )
	{
		$o_param->AssertInteger('maxlength_per_slice');
		$maxlen=$o_param->Get('maxlength_per_slice');
	}
	my $f_in = $o_param->GetStreamIn('in');
	my $cpt=0;
	my $slice=1;
	my $outp = $o_param->Get('outprefix');
	my $num  = $o_param->Get('num_per_slice');
	my $out = &GetOut($outp,$slice++);
	my $f_out = &GetStreamOut($out);
	my $sum=0;
	while(my $l=<$f_in>)
	{
		if ( $l =~ /^>/ )
		{
			if ( ( ++$cpt > $num ) || ( $maxlen && $sum > $maxlen ) )
			{
				$out = &GetOut($outp,$slice++);
				$f_out->close;
				$f_out = &GetStreamOut($out);
				$cpt=1;
				$sum=0;
			}
		}
		else
		{
			$sum+=(length($l)-1);
		}
		print $f_out $l;
	}
	$f_in->close;
	$f_out->close;
}

sub GetOut
{
my($outp,$slice)=@_;

	return sprintf("$outp%05d",$slice);
}

sub Usage
{
print STDERR<<END
$0
	--help
Mandatory
	--in            fastafilename         compressed format allowed
Optional
	--outprefix      string               output file prefix 
	--num_per_slice  integer              [$NUM_PER_SLICE] number of sequence per slice
	--maxlength_per_slice  integer        [0] maximun number of nucleotides
END
}
