#!/usr/bin/perl
#
# $Id: lipm_transfer_gff3_attributes.pl 1426 2015-03-20 15:29:34Z sallet $
#

# Copyright INRA

# Jerome.Gouzy@toulouse.inra.fr
# Erika.Sallet@toulouse.inra.fr
# Sebastien.Carrere@toulouse.inra.fr
# Ludovic.Cottret@toulouse.inra.fr
# Ludovic.Legrand@toulouse.inra.fr


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

lipm_transfer_gff3_attributes.pl 

=head1 SYNOPSIS

>lipm_transfer_gff3_attributes.pl  --in_gff3 in.gff3 --out_gff3 in.gff3.with_annot --ref_gff3 prev.gff3 

>lipm_transfer_gff3_attributes.pl  --in_gff3 in.gff3 --out_gff3 in.gff3.with_annot --ref_gff3 prev.gff3  --prefix_new BACTOCHE_

=head1 DESCRIPTION

Search the overlapped regions (same strand, fraction overlap reciprocal (-r -s) ) between input and reference GFF3 files. 
Transfer GFF3 attributes from the reference to the input and save the results in output file.
Remove operons features
Change the type if tRNA or rRNA ovelapping.

WARNING: 
* Tested only with bacteria annotations
* Tested only with EuGene formatted input files

=cut

BEGIN
{
	$ENV{LANG} = "C";
}


use FindBin;

use strict;
use warnings;
use ParamParser;
use General;
use GeneralBioinfo;
use Cwd;

my $DEBUG   = 0;
my $clean   = 1;
my $OVERLAP = 0.8;
my $DEFAULT_ATTRS      = "Name;Note;protein_id;locus_tag;product";
my $DEFAULT_RM_ATTRS   = "nb_exon;length";
my $DEFAULT_TYPES      = "CDS;tRNA;rRNA;ncRNA";
my $DEFAULT_DIGIT_NB   = 6;

MAIN:
{
    my $o_param = New ParamParser('GETOPTLONG',\&Usage, 'in=s', 'ref=s', 'out=s', 'min_overlap=f',
                                  'attrs=s', 'rm_attrs=s', 'prefix_new=s', 'digit_new=i', 'intersectbed_type=s', 'debug', 'rm_operon');

	$DEBUG = 1 if ($o_param->IsDefined('debug'));

	$o_param->AssertFileExists('in', 'ref');
	$o_param->AssertDefined('out');

	GetExePath('intersectBed');

	$o_param->SetUnlessDefined('min_overlap',       $OVERLAP);
	$o_param->SetUnlessDefined('attrs',             $DEFAULT_ATTRS);
	$o_param->SetUnlessDefined('rm_attrs',          $DEFAULT_RM_ATTRS);
	$o_param->SetUnlessDefined('prefix_new',        "");
	$o_param->SetUnlessDefined('digit_new',    $DEFAULT_DIGIT_NB);
	$o_param->SetUnlessDefined('intersectbed_type', $DEFAULT_TYPES);

	my $intersectbed_type = $o_param->Get('intersectbed_type');
	my $in      = $o_param->Get('in');
	my $ref     = $o_param->Get('ref');
	my $out     = $o_param->Get('out');
	my $overlap = $o_param->Get('min_overlap');
	
	my $workdir = getcwd;
	my $prefix   = $workdir."/seq".$$;

	my $in_reg   = $prefix.".in.gff3";
	my $ref_reg  = $prefix.".ref.gff3";
	my $ref_gene = $prefix.".ref.gene.gff3";

	# Extract genes from the reference and get some attributes
	my %h_ref_gene;
	&ExtractLines($ref, $ref_gene, "gene");
	&GetInfo($ref_gene, \%h_ref_gene);	

	# Extract CDS and ncRNA regions from input gff3
	&ExtractLines($in, $in_reg,   $intersectbed_type);
	# Extract CDS and ncRNA regions from reference GFF3
	&ExtractLines($ref, $ref_reg, $intersectbed_type);
	
	
	# Compute the intersection between cds
	my %h_info;
	&Intersect($in_reg, $ref_reg, \%h_ref_gene, \%h_info, $overlap, $clean);

	&Transfer($in, $out, \%h_info, $o_param);

	if ($clean)
	{
		unlink($in_reg, $ref_reg, $ref_gene);
	}
}

=head2 procedure GetInfo

  Title        : GetInfo
  Usage        : &GetInfo($gff3file, \%h_attr);
  Prerequisite : -
  Function     : Parse a GFF3 file and save data in a hashtable ($rh_gene). 
                 Keys=GFF3 ID, Value = ref to a hash 
  Returns      : -
  Args         : $gff3      File                    GFF3 file
				 $rh_gene   Ref to a hashtable      
  Globals      : none

=cut
sub GetInfo
{
	my ($gff3, $rh_gene) = @_;	


	my $f_ = GetStreamIn($gff3);
	while (my $l = <$f_>)
	{
		next if $l =~ /^#/;
		
		chomp($l);
		my %h_gff3;
		&GeneralBioinfo::ParseGFF3Line(\$l, 1, \%h_gff3);
		$rh_gene->{$h_gff3{ID}} = \%h_gff3;
	}
	$f_->close;

	return;
}



=head2 procedure ExtractLines

  Title        : ExtractLines
  Usage        : &ExtractLines($ingff3, $outgff3, "CDS;ncRNA")
  Prerequisite : -
  Function     : Read a GFF3 and extract only lines in specified types
  Returns      : -
  Args         : $in   GFF3 file
				 $out  GFF3 file
                 $type string      List of valid GFF3 types (separated by ;)      
  Globals      : none

=cut
sub ExtractLines
{
	my ($in, $out, $types) = @_;

	my $type_string .= ";".$types.";";

	my $f_in  = &GetStreamIn($in);
	my $f_out = &GetStreamOut($out);
	
	while (my $line = <$f_in>)
	{
		chomp($line);
		next if ($line =~ /^#/);
		my @a_f  = split(/\t/, $line);
		my $type = $a_f[2];
		next if ($type_string !~ /;$type;/);
		
		print $f_out "$line\n";
	}
	$f_in->close;
	$f_out->close;

	return;
}


=head2 procedure Intersect

  Title        : Intersect
  Usage        : &Intersect($gff3_1, $gff3_2, \%h_transfer, \%h_complement_info, "0.8", $clean) 
  Prerequisite : -
  Function     : 
  Returns      : -
  Args         : $in   GFF3 file
				 $out  GFF3 file
                 $type string      List of valid GFF3 types (separated by ;)      
  Globals      : none

=cut
sub Intersect
{
	my ($a, $b, $rh_ref_gene, $rh_info, $overlap, $clean) = @_;
	
	my $intersecout = $a.".intersec";
	my $cmd = "intersectBed -a $a -b $b -f $overlap -wo -s -r > $intersecout ";
	#print $cmd;
	system($cmd);

	my %h_ref_view;
	my %h_in_view;
	
	my $f_in  = &General::GetStreamIn($intersecout);
	while (my $line = <$f_in>)
	{
		chomp($line);
		
		my @a_ = split(/\t/, $line);
		my %h_1;		
		my %h_2;
		&GeneralBioinfo::ParseGFF3Line(\join("\t", @a_[0..8]), 1, \%h_1);
		&GeneralBioinfo::ParseGFF3Line(\join("\t", @a_[9..17]), 1, \%h_2);
		
		my $a_id = &GetId($h_1{ID}, $h_1{type});
		
		# ID=cds17;Name=NP_766658.1;Parent=gene17;Note=RSbeta%7E transposase%3B RSbeta%7E transposase;
        # Dbxref=Genbank:NP_766658.1,GeneID:1049764;gbkey=CDS;product=transposase;protein_id=NP_766658.1

		# For CDS, assert the stop position is the same
		if ( $h_1{type} eq "CDS" && $h_2{start} !=  $h_1{start} && $h_2{end} !=  $h_1{end})
		{
			print STDERR "NO TRANSFER: $line\n" if ($DEBUG);
			next;
		}

		$h_2{protein_id} = $h_2{Name};
		# complete with information from the gene parent
		if ( defined ($h_2{Parent}) && defined($rh_ref_gene->{$h_2{Parent}}))
		{
			$h_2{locus_tag}  = $rh_ref_gene->{$h_2{Parent}}->{locus_tag} if (defined ($rh_ref_gene->{$h_2{Parent}}->{locus_tag}));
			$h_2{Name}       = $rh_ref_gene->{$h_2{Parent}}->{Name}      if (defined ($rh_ref_gene->{$h_2{Parent}}->{Name}));
		}


		# check if already intersection
		if ( exists($h_ref_view{$h_2{ID}} ))
		{
			print STDERR "WARNING: the reference entry ". $h_2{ID}. " overlaps many input entries (Information are transfered on ".$h_ref_view{$h_2{ID}} ." only).\n";
			next;
		}
		if ( exists $h_in_view{$a_id} )
		{
			print STDERR "WARNING: the input entry $a_id overlaps many reference entries. (Only information of " . $h_in_view{$a_id} . " were transfered).\n>$line<\n" ;
			next;
		}
		
		# save that input id and reference id already view
		$h_in_view{$a_id}     = $h_2{ID};
		$h_ref_view{$h_2{ID}} = $a_id;
		
		# Save all the attributes
		$rh_info->{$a_id} = \%h_2;
	}
	$f_in->close;
	
	print STDERR scalar(keys(%$rh_info)) ." reference genes were transfered to input entrie\n";

	unlink($intersecout) if ($clean == 1);

	return;
}

=head2 function GetId

  Title        : GetId
  Usage        : $cleanid = &GetId($egn_id, "CDS");
  Prerequisite : -
  Function     : Compute the ID (remove eugene prefix and accroding tu type, the final prefix ".\d") 
  Returns      : The cleaning ID
  Args         : $in  string   Ex: "three_prime_UTR:NC_20001.1"
				 $out string   Ex: "NC_20001"
  Globals      : none

=cut
sub GetId
{
	my ($string, $type) = @_;
	
	my $name = $string;
	$name =~ s/^(gene|CDS|mRNA|exon|ncRNA|five_prime_UTR|three_prime_UTR)://;
	$name =~ s/\.\d$// if ($type eq "CDS" || $type eq "exon" || $type =~ /_prime_UTR/);
	#print "$string ==> $name\n";#ID : "$name

	return $name;
}


=head2 procedure Transfer

  Title        : Transfer
  Usage        : &Transfer($in_gff3, $res_gff3, \%h_info, $o_param)
  Prerequisite : -
  Function     : Transform $in gff3 file:
                 Add or substitute attributes previously collected
                 Remove attributes (operon_parent,nb_exon,length, Ontology_term)
                 Remove operon entries if 'rm_operon' option
                 if tRNA (resp rRNA), change type ncRNA to tRNA (resp. rRNA)
  Returns      : -
  Args         : $in  File   GFF3 file
				 $out File   GFF3 file
                 $rh_info     Ref to a hash
                 $o_param
  Globals      : none

=cut
sub Transfer
{
	my ($in, $out, $rh_info, $o_param) = @_;
	
	my $attrs      = $o_param->Get('attrs');
	my $rm_attrs   = $o_param->Get('rm_attrs');
	my $prefix_new = $o_param->Get('prefix_new');
	my $nb         = $o_param->Get('digit_new');

	my @a_attrs    = split(/;/, $attrs);
	my @a_rm_attrs = split(/;/, $rm_attrs);
	my $gene_cpt   = 0;

	my $f_in  = &GetStreamIn($in);
	my $f_out = &GetStreamOut($out);

	while (my $line = <$f_in>)
	{
		chomp($line);
		
		if ($line =~ /^#/)
		{
			print $f_out $line."\n";
			next;
		}
		
		my %h_1;
		&GeneralBioinfo::ParseGFF3Line(\$line, 1, \%h_1);
		
		next if ($h_1{type} eq "operon" && $o_param->IsDefined('rm_operon'));
		&RmAttributes(\%h_1, 'operon_parent') if ( $o_param->IsDefined('rm_operon'));

		foreach my $rm_att (@a_rm_attrs)
		{
			&RmAttributes(\%h_1, $rm_att);
		}

		my $id = &GetId($h_1{ID}, $h_1{type});

		if (exists($rh_info->{$id}))
		{
			foreach my $att (@a_attrs)
			{
				&AddAttributeFromHash(\%h_1, $att, $rh_info->{$id}); 
			}
			
			$h_1{type} = "tRNA" if ($h_1{type} eq "ncRNA" && $rh_info->{$id}->{type} eq "tRNA");
			$h_1{type} = "rRNA" if ($h_1{type} eq "ncRNA" && $rh_info->{$id}->{type} eq "rRNA");
		}
		else
		{
			$gene_cpt++ if ($h_1{type} eq "gene");
			if ($prefix_new ne "")
			{
				my $new_name = sprintf("$prefix_new%0${nb}d", $gene_cpt);
				&AddAttribute(\%h_1, 'Name',      $new_name);
				&AddAttribute(\%h_1, 'locus_tag', $new_name);
			}
		}
		my @a_att = ($h_1{seqid}, $h_1{source}, $h_1{type}, $h_1{start}, $h_1{end}, $h_1{score}, $h_1{strand}, $h_1{phase}, $h_1{attributes});
		my $newline = join("\t", @a_att);
		
		print $f_out "$newline\n";
	}

	$f_in->close;
	$f_out->close;
	

	print STDERR "$gene_cpt input genes do not overlap a reference gene.\n";
	return;
}




sub AddAttributeFromHash
{
	my ($rh_1, $name, $rh_info) = @_;
	
	if (exists $rh_info->{$name})
	{
		&AddAttribute($rh_1, $name, $rh_info->{$name});
	}
	return;
}


sub AddAttribute
{
	my ($rh_1, $name, $val) = @_;

	return if (!defined($val) || $val eq "");

	if ($$rh_1{attributes} =~ /;?$name=/)
	{
		$$rh_1{attributes} =~ s/$name=[^;]+/$name=$val/;
	}
	else
	{
		$$rh_1{attributes} .= ";$name=$val";
	}

	return;
}



sub RmAttributes
{
	my ($rh_1, $name) = @_;
	
	$$rh_1{attributes} =~ s/;?$name=[^;]+//;
	return;
}

sub Usage
{
    print STDERR<<END
$0
Search overlapped regions (same strand, fraction overlap reciprocal) between input and reference GFF3 files. 
Transfer GFF3 attributes from the reference to the input and save the results in output file.
WARNING: 
* Tested only with bacteria annotations
* Tested only with EuGene formatted input files

[Mandatory]
    --in                GFF3 file      Input file
    --ref               GFF3 file      Reference file
    --out               GFF3 file      Output file

[Optional]
    --attrs             quoted string  List of GFF3 attributes to transfer
                                       [$DEFAULT_ATTRS]
    --rm_attrs          quoted string  List of GFF3 attributes to remove
                                       [$DEFAULT_RM_ATTRS]    
    --min_overlap       float          [$OVERLAP] Minimum overlap required between features. 
                                       (Equivalent to -f parameter of intersectBed)
    --intersectbed_type quoted string  List of GFF3 types (third column). Search overlaps between features from these types. 
                                       [$DEFAULT_TYPES]
    --prefix_new        string         If defined, rename the new genes with this prefix, followed by a gene number.
    --digit_new         int            [$DEFAULT_DIGIT_NB] Number of digits for new gene number.
                                       (Only make sense with 'prefix_new')
    --rm_operon                        To remove operon features
    --debug
    --help
END
}
