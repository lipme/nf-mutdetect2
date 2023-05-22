#!/usr/bin/perl

#
# $Id: lipm_dbprot_remove_repbase.pl 1444 2015-06-11 11:28:52Z sallet $
#

# Copyright INRA

# Jerome.Gouzy@toulouse.inra.fr
# Erika.Sallet@toulouse.inra.fr
# Sebastien.Carrere@toulouse.inra.fr
# Emmanuel.Courcelle@toulouse.inra.fr
# Ludovic.Legrand@toulouse.inra.fr
# Ludovic.Cottret@toulouse.inra.fr
# 

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


use strict;
use warnings;
use ParamParser;
use General;
use GeneralBioinfo;

my $NCPUS          = 8;
my $REPBASE        = "/db/license/repbase/repbase_aaSeq_cleaned_TE.fa";
my $BLASTPARAM     = " -evalue 1e-6 -max_target_seqs 5 -outfmt 6  -seg yes";
my $USEARCHPARAM   = " -evalue 1e-6 -qmask seg  ";
my $BLASTPRG      = "blastp";
my $FASTAFILTERPRG = "lipm_fastafilter.pl";

MAIN:
{
	my $o_param = New ParamParser('GETOPTLONG',\&Usage,'in=s','out=s','repbase=s','blast_prg=s','blastparam=s', 'no-clean', 'cpus=i');
	$o_param->AssertFileExists('in');
	my $in  = $o_param->Get('in');
	$o_param->SetUnlessDefined('out',"$in-norep");
	$o_param->SetUnlessDefined('repbase', $REPBASE);
	$o_param->AssertFileExists('repbase');

	$o_param->SetUnlessDefined('blastparam', $BLASTPARAM);
	$o_param->SetUnlessDefined('cpus', $NCPUS);
	$o_param->SetUnlessDefined('blast_prg', $BLASTPRG);
	$o_param->AssertAllowedValue('blast_prg', 'blastp', 'usearch-blastx','usearch-blastp');
	
	$o_param->SetUnlessDefined('fastafilter_prg', $FASTAFILTERPRG);
	my $blastpg       = $o_param->Get('blast_prg');
	my $fastafilterpg = $o_param->Get('fastafilter_prg');
	$NCPUS = $o_param->Get('cpus');
	#&GetExePath($blastpg);
	&GetExePath($fastafilterpg);
	my $out        = $o_param->Get('out');
	my $repbase    = $o_param->Get('repbase');
	my $blastparam = $o_param->Get('blastparam');
	my $outtmp     = "$in.vs.repeatdb";

	# Blast+ command line
	my $cde = "";
	my $dbcol = 1;
	if ($blastpg eq 'usearch-blastx')
	{
		$dbcol = 2;
		$cde = "usearch -ublast $repbase -db $in.udb -blast6out $outtmp $USEARCHPARAM";
	}
	elsif ($blastpg eq 'usearch-blastp')
	{
		$cde = "usearch -ublast $in -db $repbase.udb -blast6out $outtmp $USEARCHPARAM";
	}
	else
	{
		$cde = "$blastpg -query $in -db $repbase $blastparam -num_threads $NCPUS -out $outtmp";
	}
	print "[CMD]\t$cde\n";
	system($cde);
		
	$cde = "$fastafilterpg --in $in --out $out.raw.gz --db_filename $outtmp --db_col $dbcol --remove";
	print "[CMD]\t$cde\n";
	system($cde);

	# suppression des | dans les ID
	my $f_in  = &GetStreamIn("$out.raw.gz");
	my $f_out = &GetStreamOut($out);
	my $ok=&FALSE;
	while(my $l=<$f_in>)
	{
		if ($l=~/^>/ )
		{
			$l =~ s,^>\S+\|,>,;
			$l =~ s,\|\S*,,;
		}
		print $f_out $l;
	}
	$f_in->close;
	$f_out->close;
	if ( ! $o_param->IsDefined('no-clean') )
	{
		unlink("$out.raw.gz");
		unlink("$outtmp");
	}
	exit 0;
}

sub Usage
{
print STDERR<<END
Remove proteins putatively related to transposable elements
$0
	--help
	--in           fastafile
[Optional]
	--out          fastafile [in.norep]
	--repbase      filename  [$REPBASE] blast indexed repeat database
	--blastparam   string    [$BLASTPARAM]
	--blast_prg	   string	 [$BLASTPRG]
	--cpus [$NCPUS]
END
}

