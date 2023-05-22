#!/usr/bin/perl

#
# $Id: lipm_get_promoter_regions.pl 1517 2016-09-27 11:18:37Z gouzy $
#

# Copyright INRA

# Jerome.Gouzy@toulouse.inra.fr
# Erika.Sallet@toulouse.inra.fr
# Sebastien.Carrere@toulouse.inra.fr
# Emmanuel.Courcelle@toulouse.inra.fr
# Ludovic.Legrand@toulouse.inra.fr
# Ludovic.Cottret@toulouse.inra.fr

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
    use File::Basename;
    my $dirprg = dirname($0);
    unshift @INC, "$dirprg/../corelib";
    $VERSION = do {my @r = (q$Rev: 1517 $ =~ /\d+/g); sprintf "%d." . "%02d" x $#r, @r}
}

use strict;
use warnings;
use ParamParser;
use General;
use GeneralBioinfo;
use File::Basename;
use Data::Dumper;
use URI::Escape;


my $DEFAULT_UPSTREAM_NT   = 1500;
my $DEFAULT_DOWNSTREAM_NT = 0;
my $DEFAULT_MINLEN = 50;

MAIN:
{
    my $o_param = New ParamParser('GETOPTLONG',\&Usage,'genome=s','output=s','help','gff3=s','upstream_nt=i','downstream_nt=i','type=s','kingdom=s','debug','minlen=i','keepID');
	$o_param->AssertFileExists('genome','gff3');
	$o_param->AssertDefined('output');
	$o_param->SetUnlessDefined('upstream_nt',$DEFAULT_UPSTREAM_NT);
	$o_param->SetUnlessDefined('downstream_nt',$DEFAULT_DOWNSTREAM_NT);
	$o_param->SetUnlessDefined('type','gene');
	$o_param->SetUnlessDefined('kingdom','eukaryota');
	$o_param->SetUnlessDefined('minlen',$DEFAULT_MINLEN);
	$o_param->AssertAllowedValue('kingdom','eukaryota');
	$o_param->AssertAllowedValue('type','gene','mRNA','pre_miRNA','ncRNA');

	my $type   = $o_param->Get('type');
	my $upnt   = $o_param->Get('upstream_nt');
	my $dwnt   = $o_param->Get('downstream_nt');
	my $minlen = $o_param->Get('minlen');

	my $debug = $o_param->IsDefined('debug');

	my %h_seq=();
	&FastaToHash($o_param->Get('genome'),\%h_seq);

	my $f_out = $o_param->GetStreamOut('output');
	my %h_res=();
	# je recupere le product et je l'attche au Parent, comme ca je n'ai toujours
	my %h_product=();
	my $f_gff3 = $o_param->GetStreamIn('gff3');
	while(my $line=<$f_gff3>)
	{
		chomp($line);
		my %h_res = ();
		next if ($line =~ /^#/ );
		&ParseGFF3Line(\$line,2,\%h_res);
		if ( defined($h_res{product}) )
		{
			if  ( defined($h_res{Parent})  )
			{
				$h_res{Parent}=~s/^\w+://;
				$h_product{$h_res{Parent}}=$h_res{product};
			}
			$h_product{$h_res{Name}}=$h_res{product};
		}
	}
	$f_gff3->close;

	$f_gff3 = $o_param->GetStreamIn('gff3');
	while(my $line=<$f_gff3>)
	{
		chomp($line);
		my %h_res = ();
		next if ($line =~ /^#/ );
		&ParseGFF3Line(\$line,2,\%h_res);
		next if ( $h_res{'type'} !~ /$type/i );

		my $chrlen = $h_seq{$h_res{seqid}}{len};
		my $start_reg=0;
		my $len=0;
		my $upnt_actual = $upnt;
		if ( $h_res{strand} eq '-' )
		{
			$start_reg = $h_res{end}+1-$dwnt;
		}
		else
		{
			$start_reg = $h_res{start}-$upnt; 
		}
		$start_reg--; # on base en base 0
		if ( $start_reg < 0 )
		{
 			$upnt_actual = $upnt-abs($start_reg);
			$start_reg = 0;
		}
		$len = $dwnt+$upnt_actual;
		if ( $start_reg+$len-1 > $chrlen )
		{
			$len = $chrlen-$start_reg;
		}
		my $raw = substr($h_seq{$h_res{seqid}}{sequence},$start_reg,$len);
		my $upstream = (defined($raw)) ? $raw : '';
		if ( $h_res{strand} eq '-' )
		{
			$upstream = &RevCompSeq(\$raw);
		}
		my $id = $h_res{ID};
		if( defined($h_res{locus_tag}) )
		{
			if ( $o_param->IsDefined('keepID') && $id !~/$h_res{locus_tag}/)
			{
				$id =~ s/^\w+://;
				$id .= ",$h_res{locus_tag}";
			}
			else
			{
				$id = $h_res{locus_tag};
			}
		}
		elsif( defined($h_res{Name}) )
		{
			if ( $o_param->IsDefined('keepID') && $id !~/$h_res{Name}/)
			{
				$id =~ s/^\w+://;
				$id .= ",$h_res{Name}";
			}
			else
			{
				$id = $h_res{Name};
			}
		}
		my $product='';
		if ( defined($h_product{$h_res{Name}}) )
		{
			$product = " ".uri_unescape($h_product{$h_res{Name}});
		}
		if (  length($upstream)  > $minlen )
		{
			if ( $debug )
			{
				printf $f_out ">${id}_%d_%d_$h_res{strand}  $type:$h_res{start}-$h_res{end} chrlen=$chrlen$product\n".&FormatSeq(\$upstream)."\n",$start_reg+1,$start_reg+1+$len-1;
			}
			else
			{
				printf $f_out ">${id}_%d_%d_$h_res{strand}$product\n".&FormatSeq(\$upstream)."\n",$start_reg+1,$start_reg+1+$len-1;
			}
		}
	}
	$f_gff3->close;
	$f_out->close;
}



sub Usage
{
    print STDERR<<END
$0
	--help

Mandatory
	--genome         fastafile
	--gff3           filename
	--output         fastafile

Optional
	--upstream_nt    integer    [$DEFAULT_UPSTREAM_NT]
	--downstream_nt  integer    [$DEFAULT_DOWNSTREAM_NT]
	--type           string     [gene], allowed values are gene, mRNA, pre_miRNA, ncRNA
	--minlen         integer    [$DEFAULT_MINLEN] minimum length of the region
	--kingdom        string     [eukaryota] 
	--keepID                    keep ID in the fasta header
	--debug                     dump more comprehensive fasta headers
	
Warning: seule la version eukaryota a ete implemente pour le moment (pour les prokaryotes il faudrait gerer les operons)
END
}

