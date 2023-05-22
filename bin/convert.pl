#!/usr/bin/perl

use strict;
use ParamParser;
use General;
use JSON;
MAIN:
{
	my $o_param = New ParamParser('GETOPTLONG', \&Usage, 'vcf=s', 'level=s');
	
	$o_param->AssertFileExists('vcf');
	
	my $vcf = $o_param->Get('vcf');
	$o_param->SetUnlessDefined('level', 'MODIFIER');
	$o_param->AssertAllowedValue('level', 'HIGH', 'MODERATE','LOW','MODIFIER');
	
	my $level = $o_param->Get('level');
	
	my @a_levels = ();

	if ($level eq 'HIGH') 	{		@a_levels = ('HIGH') ;};
	if ($level eq 'MODERATE') 	{		@a_levels = ('HIGH', 'MODERATE') ;};
	if ($level eq 'LOW') 	{		@a_levels = ('HIGH', 'MODERATE', 'LOW') ;};
	if ($level eq 'MODIFIER') 	{		@a_levels = ('HIGH', 'MODERATE', 'LOW', 'MODIFIER'); };
	
	
	
	
	my @a_lines = `bcftools query -f '%ANN\[\t%SAMPLE=%GT]\n' $vcf`;
	chomp @a_lines;
	
	my %h_gene_mutants = ();
	my %h_mutant_genes = ();
	foreach my $l (@a_lines)
	{
		$l =~ s/\t\S+=[0\.]\/[0\.]//g;
		$l =~ s/=[01]\/\d+//g;
		#print "$l\n";
		my ($ann,@a_mutants) = split ("\t", $l);
		my @a_ann = split (',', $ann);
		foreach my $ann_item (@a_ann)
		{
			my @a_ann_fields = split ('\|', $ann_item);
			
			my $eff_level = $a_ann_fields[2];
			
			if (grep (/$eff_level/, @a_levels))
			{
				my $gn = $a_ann_fields[3];
				$h_gene_mutants{$gn} = {} unless defined $h_gene_mutants{$gn};
				
				foreach my $mutant (@a_mutants)
				{
					$h_mutant_genes{$mutant} = {} unless defined $h_mutant_genes{$mutant};
					$h_mutant_genes{$mutant} ->{$gn} = 1;
					$h_gene_mutants{$gn}->{$mutant} = 1;
				}
			}
		}
	}

	my $fh_out = &GetStreamOut("$level.gene-mutants.tsv");
	foreach my $g (sort keys %h_gene_mutants)
	{
		print $fh_out "$g\t" , join (',', sort keys %{$h_gene_mutants{$g}}), "\n";
	}
	 $fh_out->close();
	 
	 $fh_out = &GetStreamOut("$level.mutant-genes.tsv");
	foreach my $m (sort keys %h_mutant_genes)
	{
		print $fh_out "$m\t" , join (',', sort keys %{$h_mutant_genes{$m}}), "\n";
	}
	 $fh_out->close();
	 
	#print to_json(\%h_gene_mutants);
	exit 0;
}

sub Usage {
	
}
