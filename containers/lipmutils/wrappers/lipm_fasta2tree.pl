#!/usr/bin/perl

#
# $Id: lipm_fasta2tree.pl 604 2011-01-07 15:18:27Z gouzy $
#

# Copyright INRA

# Erika.Sallet@toulouse.inra.fr
# Jerome.Gouzy@toulouse.inra.fr
# Sebastien.Carrere@toulouse.inra.fr
# Emmanuel.Courcelle@toulouse.inra.fr

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

=head1 NAME

lipm_fasta2tree.pl

=head1 SYNOPSIS

./lipm_fasta2tree.pl --sequence=/path/to/multifasta.file --tree=/path/to/output/directory --type=dna

=head1 DESCRIPTION

Splits a multifasta file into unique sequence subdirectories and files.

=cut

BEGIN
{
	$ENV{LANG} = "C";
	$VERSION = do { my @r = (q$Rev: 604 $ =~ /\d+/g); sprintf "%d."."%02d" x $#r, @r }
}

use strict;
use warnings;
use ParamParser;
use General;
use GeneralBioinfo;
use Carp;
use Cwd;


#lipm_fasta2tree.pl --max_seq_per_subdir <integer> --sequence <multifasta> --line_length --type <dna|rna|proteic> --tree <directory>
#ajouter len= si pas present
#tester si sequence du bon type GeneralBioinfo.pm ajouter type=  CleanSeq
#sous-repertoires 
#gerer les duplicated ID (modifier le header fasta .$count)
#generer un fichier listing des fichiers crees entries.fof.
#un fichier entries.stat (nb ID dupliques, N50, Min Max, nb de seq, nt, aa)

# Ecrire sous forme d'un module
#	GetDuplicated
#	GetMin etc

my %H_DEFAULTS = (max_subdirectories => 5000, line_length => -1, extension => '');
MAIN:
{
		&Main();
}

sub Main
{
	my $o_param = New ParamParser('GETOPTLONG',\&Usage,'sequence=s','line_length=i','type=s','clean', 'remove_empty_sequences', 'outdir=s','max_subdirectories=i','extension=s','sort','sort_by_len','sort_by_id','help','verbose');
	#$o_param->SetUnlessDefined('line_length', &GeneralBioinfo::LINE_LENGTH);
	$o_param->SetUnlessDefined('extension', $H_DEFAULTS{extension});
	$o_param->SetUnlessDefined('max_subdirectories', $H_DEFAULTS{max_subdirectories});
	$o_param->SetUnlessDefined('line_length', $H_DEFAULTS{line_length});
	
	$o_param->AssertFileExists('sequence');
	$o_param->AssertDefined('outdir');
	
	
	my $clean            = &FALSE;
	my $remove_empty_seq = &FALSE;
	
	if ($o_param->IsDefined('clean'))
	{
		$clean = &TRUE;
		$o_param->AssertAllowedValue('type','dna','rna','proteic') ;
	}
	# Flag to remove empty sequences and sequences only composed of N
	if ($o_param->IsDefined('remove_empty_sequences'))
	{
		$remove_empty_seq = &TRUE;
	}
	
	my ($ra_fastaordered_ids,$rh_sequences,$rh_duplicated,$rh_stats) = &LoadMultiFasta($o_param->GetStreamIn('sequence'),$clean, $o_param->Get('type'), $remove_empty_seq);
	my $sort = "";
	$sort = "id"  if ( $o_param->IsDefined('sort') || $o_param->IsDefined('sort_by_id'));
	$sort = "len" if ( $o_param->IsDefined('sort_by_len') ); 
	my $line_length = $o_param->Get('line_length');
	
	my $ra_sequence_files = &BuildTree($ra_fastaordered_ids,$rh_sequences,$o_param->Get('outdir'),$o_param->Get('max_subdirectories'),$line_length,$sort,$o_param->Get('extension'));
	&PrintInfoFiles($ra_sequence_files,$rh_duplicated,$rh_stats,$o_param->Get('outdir'));
	
}



=head2 procedure PrintInfoFiles

  Title        : PrintInfoFiles
  Usage        : &PrintInfoFiles($ra_sequence_files,$rh_duplicated,$rh_stats,$outdir);
  Prerequisite : NONE
  Function     : Print list of files, stats and duplicated entry ids in separate files 
  Returns      : NONE
  Args         : $ra_sequence_files : ref to the list of fasta files produced
				 $rh_duplicated : hashref of duplicated entries
				 $rh_stats : hash ref of basic statistics on sequences
				 $outdir : output directory
  Globals      : none

=cut

sub PrintInfoFiles
{
	my ($ra_sequence_files,$rh_duplicated,$rh_stats,$outdir) = @_;
	Carp::croak("ERROR - List of written files should be provided as an array reference") unless (&General::SafeIsArray($ra_sequence_files));
	Carp::croak("ERROR - Duplicated ids should be provided through a hash ref") unless (&General::SafeIsHash($rh_duplicated));
	Carp::croak("ERROR - Statistics should be provided through a hash ref") unless (&General::SafeIsHash($rh_stats));
	
	my $fh_out_entries_list = &General::GetStreamOut("$outdir/entries.fof");
	foreach my $file (@{$ra_sequence_files})
	{
			print $fh_out_entries_list $file ."\n";
	}
	$fh_out_entries_list->close;
	
	my $fh_out_duplicated_list = &General::GetStreamOut("$outdir/entries.duplicated");
	foreach my $duplicated_id (sort (keys %{$rh_duplicated}))
	{
			print $fh_out_duplicated_list $duplicated_id . "\t" . $rh_duplicated->{$duplicated_id} ."\n";
	}
	$fh_out_duplicated_list->close;
	
	my $fh_out_stats = &General::GetStreamOut("$outdir/entries.stats");
	foreach my $stat_label (sort (keys %{$rh_stats}))
	{
			print $fh_out_stats $stat_label . "=" . $rh_stats->{$stat_label} ."\n" unless($stat_label =~ /^histogram/);
	}
	$fh_out_stats->close;
	
		
	my $class = 0;
	my $step = 50;
	my @a_class = ();
	while ($class < $rh_stats->{max})
	{
		push (@a_class,$class);
		$class += $step;
	}
	my $fh_out_histo = &General::GetStreamOut("$outdir/entries.histo");
	&GeneralBioinfo::Histo($rh_stats->{histogram_data},\@a_class,$fh_out_histo);
	$fh_out_histo->close;
	
}

=head2 function BuildTree

  Title        : BuildTree
  Usage        : my $ra_files = &BuildTree($ra_fastaordered_ids,$rh_sequences,$outdir,$max_subdirectories,$line_length,$sort,$extension);
  Prerequisite : NONE
  Function     : Print sequences in subdirectories in the same order as the original fasta file or alpha-numericaly sorted
  Returns      : $ra_files a ref to the list of created file pathes
  Args         : $ra_fastaordered_ids : a ref to a list of sequence ids in the same order of the multifasta one
				 $rh_sequences : sequences hash ref
				 $outdir : output directory
				 $max_subdirectories: max number of subdirectories to create per directory
				 $line_length : sequence line length if reformat wished (-1 means one-line sequence)
				 $sort : "" (no sort), "len" (sequence length sort) or "id" (alphanumeric sort on sequence ids)
				 $extension : sequence file extension
  Globals      : none

=cut


sub BuildTree
{
		my ($ra_fastaordered_ids,$rh_sequences,$outdir,$max_subdirectories,$line_length,$sort,$extension) = @_;
		my @a_sequence_files = ();

		Carp::croak("ERROR - Multifasta ids should be provided as an array reference") unless (&General::SafeIsArray($ra_fastaordered_ids));
		Carp::croak("ERROR - Multifasta ids should be provided as an array reference") unless (&General::SafeIsHash($rh_sequences));
		
		my $raw_outdir = $outdir;
		$outdir = Cwd::abs_path($outdir) ;
		#Carp::croak("ERROR - >$outdir< is an existing file") if (-e $outdir);
		
		#$extension =~ s/^\.//;
		
		mkdir $outdir if (! -d $outdir);
		
		my @a_sequence_ids = @{$ra_fastaordered_ids};

		if ($sort eq "id")
		{
			@a_sequence_ids = sort {$a cmp $b} (keys %{$rh_sequences});
		}
		elsif ($sort eq "len")
		{
			@a_sequence_ids = sort {$$rh_sequences{$b}->{length} <=> $$rh_sequences{$a}->{length} } (keys %{$rh_sequences});
		}

		my $subdirectory = 1;
		$subdirectory = sprintf ("%04d", $subdirectory);
		
		mkdir "$outdir/$subdirectory";
		my $subdirectory_count = 0;
		foreach my $sequence_id (@a_sequence_ids)
		{
			if ($subdirectory_count == $max_subdirectories)
			{
					$subdirectory_count = 0;
					$subdirectory++;
					$subdirectory = sprintf ("%04d", $subdirectory);
					mkdir "$outdir/$subdirectory";
			}
			
			
			my $sequence_directory = "$outdir/$subdirectory/$sequence_id";
			my $sequence_raw_directory = "$raw_outdir/$subdirectory/$sequence_id";
			
			Carp::croak("ERROR - >$sequence_directory< is an existing directory") if (-d $sequence_directory);
			mkdir $sequence_directory;
			my $fh_out = &General::GetStreamOut("$sequence_directory/$sequence_id"."$extension");
			my $format_seq = $rh_sequences->{$sequence_id}->{sequence};
			Carp::croak("ERROR - sequence not found for >$sequence_id<") unless (&General::IsNotNull($format_seq));
			$format_seq = &FormatSeq(\$format_seq,$line_length) if ($line_length > 0);
			chomp $format_seq;
			print $fh_out ">" . $rh_sequences->{$sequence_id}->{header} ."\n" .$format_seq ."\n";
			$fh_out->close;
			$subdirectory_count++;
			push (@a_sequence_files,"$sequence_raw_directory/$sequence_id"."$extension");
		}
		return \@a_sequence_files;
}



=head2 function LoadMultiFasta

  Title        : LoadMultiFasta
  Usage        : my ($ra_fastaordered_ids,$rh_sequences,$rh_duplicated,$rh_stats) = &LoadMultiFasta($ra_fastaordered_ids,$rh_sequences,$rh_duplicated,$rh_stats) = &LoadMultiFasta($fh_sequence,$clean,$type);
  Prerequisite : input is a file handle
  Function     : Parse a multifasta file stream and build a perl structure (hashref)
					key : sequence_id (unique, so if diuplicated entries exist in multifasta file, .$num is added to id)
					value : hashref { header => '>header fasta', sequence => 'raw sequence', length => 'sequence length'}
  Returns      : $ra_fastaordered_ids a ref to a list of sequence ids in the same order of the multifasta one
				 $rh_sequences  sequences in hashref
                 $h_duplicated	hashref of list of duplicated entries with number of duplication
                 $rh_stats		a hashref containings stats (nb of duplicated entries, max len, min, mean, number of sequences)
  Args         : $fh_in			a multifasta filehandle
				 $clean			boolean: true is to clean sequence according to the type
				 $type			sequence type (for cleaning: dna/rna/proteic)
				 $remove_empty_seq boolean: true is to remove sequences only composed of 'N|X' or empty
  Globals      : none

=cut

sub LoadMultiFasta
{
	my ($fh_sequence,$clean,$type, $remove_empty_seq) = @_;

	my %h_sequences        = ();
	my %h_stats        = ();
	my %h_duplicated = ();
	my ($min,$max,$sum,$number) = (undef,undef,0,0);
	#pour garder le meme ordre que dans le fichier fasta
	my @a_fastaordered_ids = ();
	my @a_lengths = ();
    my $ori_separator     = $/;
    $/ = "\n>";

	Carp::croak("ERROR - input file not opened") unless defined $fh_sequence;
    if (defined $fh_sequence)
    {
        while (my $entry = <$fh_sequence>)
        {
            chomp $entry;
            if ($entry ne '')
            {
                my @a_lines = split("\n", $entry);
                my $header = $a_lines[0];
                my ($sequence_id) = ($header =~ /^>?(\S+)/);
                $sequence_id =~ s/^lcl\|//;
                
                # If header FASTA with GenBank identifier (gi) get the accession number (field 4)
               	# Ex : >gi|49175990|ref|NC_000913.2|
                if ($sequence_id =~ /^\|?gi\|/)
                {
				    $sequence_id =~ s/^\|//; # remove the first pipe if exist
				    my @a_fields = split(/\|/, $sequence_id);
				    if (defined($a_fields[3] ) && $a_fields[3] ne "" )
					{
						$sequence_id = $a_fields[3] ;
						# Remove the version number  NC_000913.2
						if ($sequence_id =~ /NC_\d+\.\d$/)
						{
							$sequence_id =~ s/\.\d+$//; # NC_000913
						}
					}
                }
               	# Substitute '|' characters by '_'
                $sequence_id =~ s/\|/_/g;

                Carp::croak("ERROR - sequence_id could not be determined for >$a_lines[0]<") unless defined $sequence_id;
                shift @a_lines;
                my $sequence = join('', @a_lines);
                my $length = length($sequence);

				if ($length <= 0)
				{
					print STDERR "WARNING - $sequence_id ignored because length=$length\n";
					next;
				}
                
                if ($remove_empty_seq == &TRUE && $sequence =~ /^([NnXx]+)$/)
                {
                    print STDERR "REMOVE $sequence_id because only composed of 'N'|'X'\n";
                    next;
                }
                
                if (exists $h_sequences{$sequence_id})
                {
					$h_duplicated{$sequence_id} = 1 unless exists $h_duplicated{$sequence_id};
					$h_duplicated{$sequence_id} ++;
					$sequence_id = $sequence_id . '.' . $h_duplicated{$sequence_id};
				}
				
				$header =~ s/^\S+/$sequence_id len=$length/;
				
				if ($clean == &TRUE)
				{
					print STDERR "CLEANING $type";
					&GeneralBioinfo::CleanSeq(\$sequence,$type) ;
				}
                
                
                
				$h_sequences{$sequence_id} = {header => $header, sequence => $sequence, length => $length};
				push (@a_fastaordered_ids,$sequence_id);
				$min = $length unless defined $min;
				$max = $length unless defined $max;
				
				$min = &General::Min($min,$length);
				$max = &General::Max($max,$length);
				$sum += $length;
				push (@a_lengths,$length);
				$number++;
            }
        }
        $fh_sequence->close;
    }
	
    $/ = $ori_separator;
    $h_stats{min} = $min;
    $h_stats{max} = $max;
    $h_stats{mean} = $sum / $number;
	$h_stats{sum} = $sum;
	$h_stats{number} = $number;
	$h_stats{duplicated} = scalar keys (%h_duplicated);
	$h_stats{histogram_data} = \@a_lengths;
    return (\@a_fastaordered_ids,\%h_sequences,\%h_duplicated,\%h_stats);
}

sub Usage
{
	print STDERR <<END;

Description:
Splits a multifasta file into unique sequence subdirectories and files.
	
Usage:
$0 --sequence=<multifasta file> --outdir=<output directory> <<END

Mandatory
	
  --sequence    		filename    Path to a multifasta file 
  --outdir     			dirname     Path to the output directory

Optional

  --line_length            integer     [sequence in one line] Length of sequence line
  --max_subdirectories     integer     [$H_DEFAULTS{max_subdirectories}] maximum number of subdirectories to be created (depending on filesystem)
  --extension              string      [$H_DEFAULTS{extension}] sequence file extension (givig xz|gzip|bz2 will compress automatically your files)
  --sort_by_len            OFF         switched on, your files will be sorted according to sequence length
  --sort_by_id, --sort     OFF         switched on, your files will be sorted according to sequence ids; otherwise multifasta order will be followed
  --clean                  OFF         switched on, your sequence will be cleaned according to their biological type (--type argument)
  --remove_empty_sequences OFF         switched on, the sequences only composed of N|X, or empty are removed
  --type                   NULL        sequence biological type : dna,rna or proteic
  
END
	exit 1;
}
