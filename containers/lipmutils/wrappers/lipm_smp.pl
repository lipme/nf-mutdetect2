#!/usr/bin/perl

# $Id$

BEGIN
{
	$ENV{LANG}='C';
	$VERSION = do { my @r = (q$Rev: 1121 $ =~ /\d+/g); sprintf "%d"."%02d" x $#r, @r }
}

=head1 NAME

lipm_smp.pl

=head1 SYNOPSIS

lipm_smp.pl --in filename --ncpus integer 

=head1 DESCRIPTION

Run command lines in parallel (one line is an atomic command)

=cut


use strict;
use warnings;
use General;
use GeneralBioinfo;
use ParamParser;


my $NB_CPUS=8;


MAIN:
{
	my $o_param = New ParamParser('GETOPTLONG',\&Usage,'in=s','ncpus=i');

	$o_param->AssertFileExists('in');
	$o_param->SetUnlessDefined('ncpus',$NB_CPUS);

	my @a_tasks=();
	my $f_task=$o_param->GetStreamIn('in');
	while(my $l=<$f_task>)
	{
		chomp($l);
		push(@a_tasks,$l);
	}
	$f_task->close;

	&MultithreadingBlock(\&System,$o_param->Get('ncpus'),\@a_tasks);
}

sub System
{
my ($imin,$imax,$ra_tasks,$ra_inputs) = @_;

	for(my $i=$imin;$i<=$imax;$i++)
	{
		system("$$ra_tasks[$i]");
	}
}

sub Usage
{
print STDERR<<END
$0
	--help
	--in	filename
	--ncpus integer  [$NB_CPUS]
END
}

