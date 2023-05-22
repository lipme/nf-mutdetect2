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
use ParamParser;
use General;
use GeneralBioinfo;
use Digest::MD5 qw/md5_hex/;
use Data::Dumper;

MAIN:
{
	my $o_param = New ParamParser('GETOPTLONG', \&Usage, 'in=s','out=s','verbose','help');
	$o_param->AssertFileExists('in');
	
	my $in = $o_param->Get('in');
	my $clean = "$in.tmp";
	system("cut -d ' ' -f1 $in | sed 's/\|/-/g' > $in.tmp");
	my %h_seq = ();
	&FastaToHash("$in.tmp", \%h_seq);
	my %h_nr = ();
	foreach my $seqid (keys %h_seq)
	{

		my $seq = uc($h_seq{$seqid}->{sequence});
		chomp $seq;
		if ($seq eq '')
		{
			print "[WARNING]\t$seqid has zero length sequence - SKIPPED\n" if ($o_param->IsDefined('verbose'));
			next;
		}
		
		my $md5 = md5_hex($seq);
		print STDERR "[DEBUG]$seqid\t$md5\t$seq\n";
		
		$h_nr{$md5} = {seq => $seq, ids => []} unless defined $h_nr{$md5};
		push(@{$h_nr{$md5}->{ids}}, $seqid);
	}
	
	my $fh_out = $o_param->GetStreamOut('out');
	foreach my $id (keys %h_nr)
	{
		print $fh_out ">$id ids=" , join(';',@{$h_nr{$id}->{ids}}) , "\n", $h_nr{$id}->{seq}, "\n";
	}
	$fh_out->close();
	
	unlink ("$in.tmp");
	exit 0;
	
}

sub Usage
{
	print STDOUT <<END;
	
	$0 --in <fasta> --out <fasta>
	
END
}
