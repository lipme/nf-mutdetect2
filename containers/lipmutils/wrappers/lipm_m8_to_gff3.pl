#!/usr/bin/perl

#
# $Id$
#

# Copyright INRA

# Erika.Sallet@toulouse.inra.fr
# Jerome.Gouzy@toulouse.inra.fr
# Sebastien.Carrere@toulouse.inra.fr
# Emmanuel.Courcelle@toulouse.inra.fr
# Ludovic.Cottret@toulouse.inra.fr
# Ludovic.Legrand@toulouse.inra.fr

# This software is a computer program whose purpose is to provide a
# web-based system for detecting, clustering and annotating non protein
# coding RNA genes.

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

=head1

=head1 NAME

    lipm_m8_to_gff3.pl

=head1 SYNOPSIS

   ./lipm_m8_to_gff3.pl --infile file.m8

=head1 DESCRIPTION

	Convert m8 file in gff3 format. Link match parts of a protein with match/match_part.

=head1 MANDATORY

 --infile        filename       m8 file

=head1 OPTIONAL

  --help this message

  --outfile         filename 	[infile.gff3] Result GFF3 file
  --source          string      [blastx] Name of the databank from which the subject sequence comes
  --type            string      [match|EST_match|protein_match|cDNA_match]. Type of the match (corresponding to a SO term)
   --seq_column     integer     The column of the sequence ID     (warning: the first column is 1, default is 1)
  --max_intron_len  integer     [default 5000]

=head1

=cut


BEGIN
{
	$ENV{LANG} = "C";
}

use strict;
use warnings;
use ParamParser;
use General;
use GeneralBioinfo;


my $MAX_INTRON_LEN_DEF=5000;
my $BUFFER_SIZE='1G';

MAIN:
{
    my $o_param = New ParamParser;
    $o_param->SetUsage(my $usage = sub {&Usage();});
    $o_param->Update( "GETOPTLONG", "I", ( "infile=s", "outfile=s", "help:s", "source=s", "type=s", "max_intron_len=s", "seq_column=i"));
    $o_param->SetBehaviour('assert_empty_file_allowed') ;

    $o_param->AssertFileExists('infile');
    my $m8file = $o_param->Get('infile');
    $o_param->SetUnlessDefined('outfile', $m8file . ".gff3");
    my $gff3file = $o_param->Get('outfile');
    my ($name_seq, $num_blast) = $m8file =~ /.*\/(.+)\.blast(\d)/;
    ($name_seq, $num_blast) = $m8file =~ /(.+)\.blast(\d)/ if (!defined($name_seq));
    if (defined  $num_blast)
    {
	$o_param->SetUnlessDefined('source', "blastx" . $num_blast);
    }
    else
    {
	$o_param->SetUnlessDefined('source', "blastx");
    }
    $o_param->SetUnlessDefined('seq_column',      1);
    $o_param->SetUnlessDefined('max_intron_len', $MAX_INTRON_LEN_DEF);
    $o_param->AssertInteger('max_intron_len');
    $o_param->SetUnlessDefined('type', "match");
    $o_param->AssertAllowedValue('type',       ('match', 'EST_match', 'protein_match', 'cDNA_match'));
    $o_param->AssertAllowedValue('seq_column', (1,2));
    
    my $source         = $o_param->Get('source');
    my $type           = $o_param->Get('type');
    my $max_intron_len = $o_param->Get('max_intron_len');
    my $seq_column     = $o_param->Get('seq_column');
    my $sub_type       = "match_part";
    my $ok = 0;
    #  -- Sort the input file --
    # 1 Supprime les lignes redondantes (meme query, meme target, meme position sur la query, meme position sur la target)
    # 2 Tri sur le nom de la proteine, puis sur la position de debut sur le genome, 
    # puis sur la position de debut de la proteine : correspond a ce qui est fait dans le pipeline EGNP
    my $m8filetmp =  $gff3file. time.$$;
    #&General::UnixSort($m8file, $m8filetmp, ["2", "7n", "9n"]); 
    my $cmd = "sort --buffer-size=$BUFFER_SIZE -u -k1,2 -k7,10 $m8file| sort --buffer-size=$BUFFER_SIZE -n -k9,9| sort --buffer-size=$BUFFER_SIZE -s -n -k7,7| sort --buffer-size=$BUFFER_SIZE -s -k2,2 >  $m8filetmp";
    system($cmd);

    #my %h_contigs_list;
    #my $rh_contigs = \%h_contigs_list;    # reference sur le hachage
    #my @contig_hit;
    
    my $hr_subject;

    #my $nb_subject_num       = -1;
    my $subject_previous_end = 0;
    my $strand_previous      = 0;
    my $subject_previous     = "";
    my $tstart_previous = 0;
    my $tend_previous   = 0;

    my $num_dafgroup = 0;
    my $seqname;
    
    my $f_gff3 = &General::GetStreamOut($gff3file);
    my $f_m8   = &General::GetStreamIn($m8filetmp);
    while (my $line = <$f_m8>)
    {
        # -- exemple :
	#GMI11495-Rm2011G.a	Q930Y0_CH603_RHIME	100.00	544	0	0	67828	66197	1	544	0.0	1182
	#GMI11495-Rm2011G.a	Q930Y0_CH603_RHIME	78.57	308	66	0	400000	399077	115	422	0.0	 529
        #GMI11495-Rm2011G.a	Q930Y0_CH603_RHIME	53.85	117	54	0	399103	398753	415	531	0.0	 144
	#GMI11495-Rm2011G.a	A6UH06_CH604_SINMW	90.07	544	54	0	67828	66197	1	544	0.0	1078
	#GMI11495-Rm2011G.a	A6UH06_CH604_SINMW	77.27	308	70	0	400000	399077	115	422	5e-146	 526
	
        chomp $line;
	next if ($line =~ /^#/);
	
	my %h_res;
	&GeneralBioinfo::ParseM8Line(\$line, $seq_column, \%h_res);
	$h_res{target} =~ s/=|;/_/g;
	
	# new match
	if ( ($h_res{target} ne $subject_previous)                           || 
	     (abs($h_res{qstart}-$subject_previous_end) >=  $max_intron_len) || 
	     ($h_res{strand} ne $strand_previous )                           ||
	     ($tstart_previous != 0 && $h_res{strand} eq "+" && (  ($h_res{tstart}   <  $tstart_previous)     || 
								   ($h_res{tend}     <= $tend_previous  )     || 
								   (($tend_previous - $h_res{tstart}) >= 20) ) 
	     ) ||
	     ($tstart_previous != 0 && $h_res{strand} eq "-" && (  ($tstart_previous <  $h_res{tstart})       ||
								   ($tend_previous   <= $h_res{tend}  )       ||
								   (($h_res{tend} - $tstart_previous) >= 20) )
	     )
	   )
	{
	    
	    # save previous 
	    if ($ok == 1)
	    {		
		$num_dafgroup++;
		&PrintMatch($hr_subject, $num_dafgroup, $f_gff3, $source, $type, $sub_type);
	    }
	    
	    # Cleaning before next match
	    
            #$nb_subject_num++; # utile???
	    $hr_subject->{contig_name} = $h_res{query};
	    $hr_subject->{subject_name} = $h_res{target};
	    $hr_subject->{subject_lgth} = "";
	    $hr_subject->{subject_sens} = $h_res{strand};
	    $hr_subject->{subject_start} = 10000000000;
	    $hr_subject->{subject_stop} = 0;
	    $hr_subject->{match} = [];

	    
            $subject_previous = $h_res{target};
            $strand_previous =  $h_res{strand}; 
	    $tstart_previous = 0;
	    $tend_previous   = 0;
        }
	

	$ok = 1;
        # on remplit le tableau C car nouveau match
        push @{$hr_subject->{match}},
		{
            contig_start  => $h_res{qstart},
            contig_stop   => $h_res{qend},
            evalue        => $h_res{"e-value"}, 
            sens          => $h_res{strand},
            subject_start => $h_res{tstart},
            subject_stop  => $h_res{tend},
            phase         => "0",
            score_hit     => $h_res{score}
		};
	
	$subject_previous_end =  $h_res{qend};
	$tstart_previous      = $h_res{tstart};	
	$tend_previous        = $h_res{tend};
	
        # determination min(contig_start) et max(contig_stop) pour un subject
        if ($h_res{qstart} < $hr_subject->{'subject_start'})
        {
            $hr_subject->{'subject_start'} = $h_res{qstart};
        }
        if ($h_res{qend} > $hr_subject->{'subject_stop'})
        {
            $hr_subject->{'subject_stop'} = $h_res{qend};
        }
    }

    if ($ok)
    {
	$num_dafgroup++;
	&PrintMatch($hr_subject, $num_dafgroup, $f_gff3, $source, $type, $sub_type);
    }
    
    $f_m8->close;
    unlink($m8filetmp); # remove temporary file
    $f_gff3->close();
    # print "File created : $gff3file \n";
}




sub PrintMatch
{
    my ($hr_subject, $num_dafgroup, $f_gff3, $source, $type, $sub_type) = @_; 

    #$num_dafgroup++;                    #a faire avant!!!
    my $num_dnadnaalign = 0;
    my $seqname = $hr_subject->{'contig_name'};
    # on ecrit une ligne correspondant au "Parent" du match (subject)
    print $f_gff3 $seqname . "\t" 
	. $source . "\t"
	. $type . "\t"
	. $hr_subject->{'subject_start'} . "\t"
	. $hr_subject->{'subject_stop'} . "\t.\t"
	. $hr_subject->{'subject_sens'}
    . "\t.\tID="
	. $source . ":"
	. $seqname . "."
	. $num_dafgroup
	. ";Dbxref="
	. $source . ":"
	. $hr_subject->{'subject_name'} 
    . ";Name=".$hr_subject->{'subject_name'}. "\n";
    
    foreach my $hr_match (@{$hr_subject->{'match'}})    # pour chaque match
    {
	$num_dnadnaalign++;
	my $lgth_match = $hr_match->{'subject_stop'} - $hr_match->{'subject_start'} + 1;
	
	print $f_gff3 $seqname . "\t" 
	    . $source . "\t"
	    . $sub_type . "\t"
	    . $hr_match->{'contig_start'} . "\t"
	    . $hr_match->{'contig_stop'} . "\t"
	    . $hr_match->{'evalue'} . "\t"
	    . $hr_match->{'sens'}
	. "\t.\tID="
	    . $source . ":"
	    . $seqname . "."
	    . $num_dafgroup . "."
	    . $num_dnadnaalign
	    . ";Dbxref="
	    . $source . ":"
	    . $hr_subject->{'subject_name'}
	. ";Target="
	    . $hr_subject->{'subject_name'} . " "
	    . $hr_match->{'subject_start'} . " "
	    . $hr_match->{'subject_stop'}
	. " +;Gap=M"
	    . $lgth_match
	    . ";frame_hit="
	    . $hr_match->{'phase'}
	. ";score_hit="
	    . $hr_match->{'score_hit'}
	. ";Parent="
	    . $source . ":"
	    . $seqname . "."
	    . $num_dafgroup . "\n";
    }
    
    return;
}


# sort m8 file according to target name, genomic position then target position, 
# and save result in m8sortfile
sub Sort
{
	my ($m8file, $m8sortfile) = @_;
	
	my @a_2sort;
	my $f_m8 = &General::GetStreamIn($m8file);
    while (my $line = <$f_m8>)
    {
		chomp($line);
		push(@a_2sort, $line);
	}
	$f_m8->close;
	my @a_sort_line = sort m8sort @a_2sort;

	my $f_m8sort = &General::GetStreamOut($m8sortfile);
	foreach my $line (@a_sort_line)
	{
		print $f_m8sort "$line\n";
	}
	$f_m8sort->close;	
	
	return;
}


# tri sur le nom de la proteine, puis sur la position de debut sur le genome, 
# puis sur la position de debut de la proteine : correspond a ce qui est fait dans le pipeline EGNP
sub m8sort
{
	my %A; 
	my %B; 
	&GeneralBioinfo::ParseM8Line(\$a, 1, \%A);
	&GeneralBioinfo::ParseM8Line(\$b, 1, \%B);
	
    #GMI11495-Rm2011G.a	Q930Y0_CH603_RHIME	100.00	544	0	0	67828	66197	1	544	0.0	1182
	#GMI11495-Rm2011G.a	Q930Y0_CH603_RHIME	78.57	308	66	0	400000	399077	115	422	0.0	 529
	if ($A{target} eq $B{target}) # If same protein
	{
		if ($A{qstart} == $B{qstart}) # same genomic start position
		{
			return $A{tstart} <=> $B{tstart}; # sort by protein start position
		}
		else
		{
			return $A{qstart} <=> $B{qstart}; # sort by genomic start position
		}
	}
	else
	{
		return $A{target} cmp $B{target}; # sort by protein name 
	}

	return;	
}



=item Usage

	$o_param->Usage();
	Print the usage of the program

=cut

sub Usage
{
    print STDERR<<END
$0     
Convert into GFF3 format a m8 file. 
Works only with a unique genomic sequence

=head1 MANDATORY

 --infile        filename       m8 file

=head1 OPTIONAL

  --help this message

  --outfile         filename 	[infile.gff3] Result GFF3 file
  --source          string      [blastx] Name of the databank from which the subject sequence comes
  --type            string      [match|EST_match|protein_match|cDNA_match]. Type of the match (corresponding to a SO term)
  --seq_column      integer     [1] The column of the sequence ID  (warning: the first column is 1)
  --max_intron_len  integer     [default 5000]
END
;

}

