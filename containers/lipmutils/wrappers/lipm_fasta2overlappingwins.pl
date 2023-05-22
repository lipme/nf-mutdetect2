#!/usr/bin/perl

# $Id: lipm_fasta2overlappingwins.pl 1468 2015-10-02 12:56:14Z  $


BEGIN
{
    $ENV{LANG} = "C";
}




use strict;
use warnings;
use General;
use GeneralBioinfo;
use ParamParser;

my $DEF_MAXRDLEN = 1999;

MAIN:
{
	my $o_param = New ParamParser('GETOPTLONG',\&Usage,'in=s','out=s','maxlen=i','step=i');
	$o_param->AssertFileExists('in');
	$o_param->AssertDefined('out');
	$o_param->SetUnlessDefined('maxlen',$DEF_MAXRDLEN);
	my $maxrdlen = $o_param->Get('maxlen');
	$o_param->SetUnlessDefined('step',$maxrdlen/2);
	my $step = $o_param->Get('step');

	my %h_seq = ();
	&FastaToHash($o_param->Get('in'),\%h_seq);
	my $f_out = $o_param->GetStreamOut('out');
	my $seq = undef;
	foreach my $seqid (sort keys(%h_seq))
	{
		my $seq = $h_seq{$seqid}{sequence};
		&Dump($f_out,$seqid,\$seq,$maxrdlen,$step);
	}
	$f_out->close;
}

sub Dump
{
my($f_out,$id,$r_seq,$maxrdlen,$step)=@_;
	
	my $seqlen = length($$r_seq);
	if ( $seqlen <= $maxrdlen )
	{
		my $begin = sprintf('%09d',1);
		print $f_out ">$id:1-$seqlen\n$$r_seq\n";
	}
	else
	{
		my $cpt=1;
		for(my $i=0;($i+$maxrdlen-1)<$seqlen;$i+=$step)
		{
			my $sseq = substr($$r_seq,$i,$maxrdlen);
			my $begin =  sprintf('%09d',$i +1) ;
			my $end = $begin + $maxrdlen - 1 ;
			print $f_out ">$id:$begin-$end\n$sseq\n";
			$cpt++;
		}
		# je veux que la fin de seq ait aussi un ovl max
		#SC:20220207 Il me manquait la derniere base d'ou le $maxrdlen + step
		my $sseq = substr($$r_seq,$seqlen-$maxrdlen-1,$maxrdlen + $step);
		my $begin = sprintf('%09d',$seqlen-$maxrdlen) ;
		my $end = $seqlen ;
		print $f_out ">$id:$begin-$end\n$sseq\n";
	}
}

sub Usage
{
print STDERR<<END
$0
	--help
	--in fastafile   input
	--out fastafile  output
[Optional]
	--maxlen  integer  maximum length of a sequence
	--step    integer  overlap length
END
}
