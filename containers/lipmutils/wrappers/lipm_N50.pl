#!/usr/bin/perl

#
# $Id: lipm_N50.pl 1446 2015-06-24 08:20:29Z gouzy $
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

lipm_N50.pl

=head1 SYNOPSIS

lipm_N50.pl --in fastafile [--noN] [--histo filename --histostep integer] [--histomax integer]


=head1 DESCRIPTION


=cut


BEGIN
{
		$ENV{LANG} = "C";
        $VERSION = do { my @r = (q$Rev: 1446 $ =~ /\d+/g); sprintf "%d."."%02d" x $#r, @r }
}

use strict;
use warnings;
use Carp;
use ParamParser;
use General;
use GeneralBioinfo;

MAIN:
{
	my $o_param = New ParamParser('GETOPTLONG',\&Usage,'help','in=s','noN','Nstats','histo=s','histomax=i','histostep=i','out=s','minlen=i','format=s');
	if ( ! $o_param->IsDefined('in') && defined($ARGV[0]) )
	{
		$o_param->Set('in',$ARGV[0]);
	}
	$o_param->AssertDefined('in');
	$o_param->SetUnlessDefined('out','stdout');
	$o_param->SetUnlessDefined('minlen',0);
	$o_param->SetUnlessDefined('format','fasta');
	$o_param->AssertAllowedValue('format','fasta','fastq');
	my $f_in = $o_param->GetStreamIn('in');
	my $f_out = $o_param->GetStreamOut('out');
	my $noN     = ( $o_param->IsDefined('noN') ) ? 1 : 0;
	my $f_histo = undef;
	my $f_histo_noN = undef;
	if ( $o_param->IsDefined('histo') )
	{
		$o_param->AssertDefined('histostep');
		$f_histo=$o_param->GetStreamOut('histo');
		$f_histo_noN=&GetStreamOut($o_param->Get('histo').'noN');
	}
	
	my %h_len = ();
	my %h_len_noN = ();
	my $id = ();
	my $ra_nnn_stretches  = [];
	my $one_line_seq = '';
	my @a_nnn_stretches = ();

	my $chr = '>';
	my $s = 2;
	if ( $o_param->Get('format') eq 'fastq' )
	{
		$chr = '@';
		$s= 4;
	}

	my $seq_ok = 0;
	my $cpt_line = -1;
	while(my $l=<$f_in>)
	{
		chomp($l);
		$cpt_line++;
		if ( $l =~ /^$chr(\S+)/ && ($s != 4 || ($cpt_line % $s) == 0) ) # c'est ok si ce n'est pas fastq ou si c'est mod 4
		{
			if ( $o_param->IsDefined('Nstats') )
			{
				push (@a_nnn_stretches,@{&GetNNNStretches($one_line_seq)}) if (defined $id && $id ne '' && $one_line_seq ne '');
			}
			$id = $1;
			$one_line_seq = '';
			$seq_ok=1;
		}
		elsif ( $seq_ok )
		{
			$h_len{$id} += length($l);
			$one_line_seq .= $l;
			$seq_ok=0 if ( $s == 4 ); # pour les fastq on ne gere que les simples lignes
			next if ( ! $noN );
			$l =~ s/[NX\s]//ig;
			$h_len_noN{$id} += length($l);
		}
	}
	$f_in->close;

	if ( $o_param->IsDefined('Nstats') )
	{
		push (@a_nnn_stretches,@{&GetNNNStretches($one_line_seq)}) if (defined $id && $id ne '' && $one_line_seq ne '');
	}

	
	
	
	
	my $minlen = $o_param->Get('minlen');
	if ( $minlen  > 0 )
	{
		foreach my $id (keys(%h_len))
		{
			if ( $h_len{$id} < $minlen )
			{
				delete($h_len{$id});
				delete($h_len_noN{$id});
			}
		}
	}
	my @a_len     = sort { $b <=> $a } values (%h_len);
	
	
	&Print('',\@a_len,$f_out);
	
	my @a_len_noN = sort { $b <=> $a } values (%h_len_noN);
	if ( $#a_len_noN != -1)
	{
		print $f_out "=========================\n";
		&Print("-noN",\@a_len_noN,$f_out);
	}
	&PrintNNNStats(\@a_nnn_stretches,$f_out) if ( $o_param->IsDefined('Nstats') );
	
	exit if ( ! $o_param->IsDefined('histo') );
	my $step    = $o_param->Get('histostep');
	my $max     = ( $o_param->IsDefined('histomax') ) ? $o_param->Get('histomax') : $a_len[0];
	my @a_class = ();
	my $len = 0;
	while ( $len < $max )
	{
		push(@a_class,$len);
		$len+=$step;
	}
	push(@a_class,$len);
	&Histo(\@a_len,\@a_class,$f_histo);
	if ( $#a_len_noN != -1)
	{
		&Histo(\@a_len_noN,\@a_class,$f_histo_noN);
	}
	
	
}

sub PrintNNNStats
{
	my ($ra_nnn_stretches,$f_out) = @_;
	
	my @a_nnn_stretches = sort {$b cmp $a} @{$ra_nnn_stretches};
	my $num_nnn = scalar @a_nnn_stretches;
	my $total_nnn = 0;
	$total_nnn = length (join ('',@a_nnn_stretches)) if ($num_nnn > 0);
	my $mean_nnn = 0;
	$mean_nnn = sprintf ("%.2f", $total_nnn / $num_nnn ) if ($num_nnn > 0);
	my $median_nnn = 0;
	$median_nnn = length ( $a_nnn_stretches[($num_nnn/2)]) if ($num_nnn > 0);
	
	
	
	my $max_nnn = 0;
	$max_nnn = length(shift (@a_nnn_stretches)) if ($num_nnn > 0);
	my $min_nnn = 0;
	$min_nnn = length(pop (@a_nnn_stretches)) if ($num_nnn > 0);
	
	print $f_out "=========================\n";
	print $f_out "NUM-NNN\t$num_nnn\n";
	print $f_out "MAX-NNN\t$max_nnn\n";
	print $f_out "MIN-NNN\t$min_nnn\n";
	print $f_out "MEAN-NNN\t$mean_nnn\n";
	print $f_out "MEDIAN-NNN\t$median_nnn\n";
	print $f_out "TOTAL-NNN\t$total_nnn\n";
	
	
	
	return;

}


sub GetNNNStretches
{
	my $one_line_seq = shift;
	my @a_nnn_stretches = ();
	while ($one_line_seq =~ /(n+)/gi)
	{
		my $stretch = $1;
		push (@a_nnn_stretches, $stretch);
	}
	return \@a_nnn_stretches;
}

sub Print
{
my($tag,$ra_len,$f_out)=@_;
	
	my $nb = scalar(@$ra_len);
	return if ( ! $nb );
	print $f_out "NUM$tag\t$nb\n";
	print $f_out "MIN$tag\t".$$ra_len[$nb-1]."\n";
	print $f_out "MAX$tag\t".$$ra_len[0]."\n";
	my $sum=0;
	for(my $i=0;$i<$nb;$i++)
	{
		$sum+=$$ra_len[$i];
	}
	my $limit = $sum/2;
	my $n50 = 0;
	my $nb_n50 = 0;
	my $misum = 0;
	foreach my $len (@$ra_len)
	{
		last if( $misum > $limit );
		$misum += $len;
		$n50 = $len;
		$nb_n50++;
	}
	
	$limit = 0.9 * $sum;
	my $n90 = 0;
	my $nb_n90 = 0;
	$misum = 0;
	foreach my $len (@$ra_len)
	{
		last if( $misum > $limit );
		$misum += $len;
		$n90 = $len;
		$nb_n90++;
	}
	
	print $f_out "N50 BP$tag\t$n50\n";
	print $f_out "N50 NUM$tag\t$nb_n50\n";
	print $f_out "N90 BP$tag\t$n90\n";
	print $f_out "N90 NUM$tag\t$nb_n90\n";
	
	printf $f_out "MEAN$tag\t%d\n",$sum/$nb;
	print $f_out "MEDIAN$tag\t".$$ra_len[$nb/2]."\n";
	print $f_out "BP$tag\t$sum\n";
}



sub Usage
{
print STDERR<<END
$0
	--help
Mandatory
	--in        fastafilename         stream or compressed format allowed
Optional
	--out       filename              [stdout]
	--histo     filename              file name of the distribution of lengths
	--histomax integer                max of the distribution
	--histostep integer               size of the intervals
	--noN                             flag, do the same analyses without considering Ns
	--minlen    integer               do not consider sequences smaller than this value [0]
	--Nstats                          flag, compute analyses on NNN stretches
END
}
