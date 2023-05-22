#!/usr/bin/perl

#
# $Id: lipm_bed_split_by_sequence.pl 1074 2013-02-25 17:24:42Z gouzy $
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

lipm_bed_split_by_sequence.pl

=head1 SYNOPSIS

lipm_bed_split_by_sequence.pl --in filename

=head1 DESCRIPTION

splits a bed/bedpe. It creates one file for each sequence names (sorted by start position)

=cut


BEGIN
{
		$ENV{LANG} = "C";
        $VERSION = do { my @r = (q$Rev: 1074 $ =~ /\d+/g); sprintf "%d."."%02d" x $#r, @r }
}


use strict;
use warnings;
use ParamParser;
use General;
use File::Basename;

MAIN:
{
	my $o_param = New ParamParser('GETOPTLONG',\&Usage,'in=s','mode=s');
	$o_param->AssertFileExists('in');
	my $in = $o_param->Get('in');
	# 1) tri par sequence pour ne pas ouvrir/fermer en permance les fichiers
	&UnixSort($in,"$in.sorted.gz","1");
	my $seqid= '';
	my $f_in = &GetStreamIn("$in.sorted.gz");
	my $bin = basename($in);
	my $f_out= undef;
	my %h_files = ();
	my $pid = $$;
	while(my $l=<$f_in>)
	{
		chomp($l);
		next if ( $l =~ /^#/ );
		my @a_=split(' ',$l);
		if ( ! defined($f_out) || $a_[0] ne $seqid )
		{
			$f_out->close if ( defined($f_out) );
			$seqid = $a_[0];
			$f_out = &GetStreamOut("$seqid.raw$pid.$bin");
			$h_files{"$seqid.raw$pid.$bin"}=1;
		}	
		print $f_out "$l\n";
	}
	$f_in->close;
	$f_out->close if (defined($f_out) );

	foreach my $raw_file (sort keys(%h_files))
	{
		my $file = $raw_file;
		$file =~ s/\.raw$pid\././;
		# tri par position croissante
		&UnixSort($raw_file,$file,"2n");
		unlink $raw_file;
	}
	unlink "$in.sorted.gz";
}


sub Usage
{
print STDERR<<END
$0 
	--help
	--in    filename      bed or bedpe formatted file
END
}
