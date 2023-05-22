#!/usr/bin/perl

#
# $Id: lipm_fastafilter.pl 1445 2015-06-22 14:14:26Z sallet $
#

# Copyright INRA

# Sebastien.Carrere@toulouse.inra.fr
# Emmanuel.Courcelle@toulouse.inra.fr
# Erika.Sallet@toulouse.inra.fr
# Jerome.Gouzy@toulouse.inra.fr

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

lipm_fastafilter.pl

=head1 SYNOPSIS


# remove sequences that are not related to db_filename
 
lipm_fastafilter.pl --in multifastafile --out multifastafile --keep --db_filename tabulated_filename

# remove sequences related to db_filename

lipm_fastafilter.pl --in filename --out filename --remove --db_filename tabulated_filename

# remove sequences matching the given keyword

lipm_fastafilter.pl --in filename --out filename --keyword string

# extract sequences based on their lengths

lipm_fastafilter.pl --in filename --out filename --min_length integer --max_length integer

# extract subsequences

lipm_fastafilter.pl --in filename --out filename --substr --db_filename m8_tabulated_file --db_column number --cols number,number

=head1 DESCRIPTION

filter out any file in fasta format 

=cut

BEGIN
{
    $ENV{LANG} = "C";
    $VERSION = do {my @r = (q$Rev: 1445 $ =~ /\d+/g); sprintf "%d." . "%02d" x $#r, @r}
}

use strict;
use warnings;
use Carp;

use ParamParser;
use IO::File;
use Data::Dumper;
use GeneralBioinfo;
use constant SEP     => "%_";
use constant SEP_POS => ",";

MAIN:
{
    my $o_param = New ParamParser(
                                  'GETOPTLONG',  "in=s",     "out=s",    "db_filename=s",
                                  "db_column=i", "min_length=i", "minlen=i","max_length=i", "maxlen=i","keyword=s", "bothends=s",
                                  "cols=s",      "keep",         "remove",       "substr", "mask"
                                  );

    $o_param->SetUsage(my $usage = sub {&Usage();});
	$o_param->SetUnlessDefined('out','stdout');
    $o_param->AssertDefined('in','out');
    my $keyword   = $o_param->Get('keyword');
	if ( $o_param->IsDefined('maxlen') )	# compat
	{
		$o_param->SetUnlessDefined('max_length',$o_param->Get('maxlen'));
	}
	if ( $o_param->IsDefined('minlen') )	# compat
	{
		$o_param->SetUnlessDefined('min_length',$o_param->Get('minlen'));
	}
    my $max_len   = $o_param->Get('max_length');
    my $min_len   = $o_param->Get('min_length');
    $keyword   = undef if ($keyword   eq '');
    $min_len   = undef if ($min_len   eq '');
    $max_len   = undef if ($max_len   eq '');
    my $start_col = undef;
    my $stop_col  = undef;
    my $tgt_col   = undef;

    my $sep = &SEP_POS;

    # To extract subsequences, column positions have to be defined
    if ($o_param->IsDefined('substr') || $o_param->IsDefined('mask'))
    {
        $o_param->AssertDefined('cols');
        my $cols = $o_param->Get('cols');
        if ($cols !~ /^\d+$sep\d+$sep\d+$/)
        {
            Carp::croak("Error: bad --cols format. Example: 8,9,2 => start is column 8, end is at 9 and the target ID is column 2\n");
        }
        my @a_col = split(/$sep/, $cols);
        $start_col = $a_col[0]-1;
        $stop_col  = $a_col[1]-1;
        $tgt_col  = $a_col[2]-1;
    }

    my %h_list_of_ids = ();
    my %h_cols        = ();    # Hash. ID: name, Value: positions of substrings

	my $bothends = undef;
	if ( $o_param->IsDefined('bothends') )
	{
		$bothends = $o_param->Get('bothends');
	}

    if ($o_param->IsDefined('db_filename'))
    {
        $o_param->AssertFileExists('db_filename');
        my $f_src = $o_param->GetStreamIn('db_filename');
        while (my $l = <$f_src>)
        {
            chomp($l);
            my @a_ = split(/\t/, $l);
            my $name = '';
            if ($l =~ /^>(\S+)/ || $l =~ /^(\S+)$/)    # fasta format ||  list of ids (one per line)
            {
                $name = $1;
            }
            else                                       # m8 or m9
            {
                next if ($l =~ /^#/);
                $o_param->AssertInteger('db_column');
                my $db_column = $o_param->Get('db_column');
				$db_column--;
                $name = $a_[$db_column];
               	#Clean du namespace blast lcl|
               	$name =~ s/^lcl\|//;
                if ($o_param->IsDefined('substr') || $o_param->IsDefined('mask'))
                {
					my $tgt = $a_[$tgt_col];
					#Clean du namespace blast lcl|
					$tgt =~ s/^lcl\|//;
                    # Get the positions of the subseq to extract and save it in the hash
                    my ($start, $stop,$strand) = &MinMaxStrand($a_[$start_col], $a_[$stop_col]);
                    $h_cols{$name} .= $start . &SEP_POS . $stop . &SEP_POS. $tgt . &SEP_POS. $strand. &SEP;
                }
            }
			if (defined($bothends))
			{
				my ($nm) = $name =~ /$bothends/o;
				$name = $nm;
			}
			#SC20150226
			#Clean du namespace blast lcl|
			$name =~ s/^lcl\|//;
            $h_list_of_ids{$name} = 1;
        }
        $f_src->close;
    }

    my $mode = undef;
    if ($o_param->Get('keep'))    # keep > remove
    {
        $mode = 'keep';
    }
    elsif ($o_param->Get('remove'))
    {
        $mode = 'remove';
    }
    elsif ($o_param->Get('substr'))
    {
        $mode = 'substr';
    }
    elsif ($o_param->Get('mask'))
    {
        $mode = 'mask';
    }

    my %h_seq_length = ();
    
    my %h_seq = ();
    &FastaToHash($o_param->GetStreamIn('in'),\%h_seq);
    
    my $f_outfile = $o_param->GetStreamOut('out');
    # For each sequence of the fasta file
    foreach my $name (keys %h_seq)
    {
        #SC20150226
        #Clean du namespace blast lcl|
        
        $name =~ s/^lcl\|//;
		my $nameinsert = $name;
		if (defined($bothends))
		{
			($nameinsert) = $nameinsert =~ /$bothends/o;
		}
		
        my $desc = $h_seq{$name}->{header};
        my $len  = $h_seq{$name}->{len};
        next if (defined($min_len) && ($len < $min_len));
        next if (defined($max_len) && ($len > $max_len));
        next if (defined($keyword) && ($name !~ /$keyword/ && $desc !~ /$keyword/));

		# 1) keep/remove filters
        if (defined($mode))
        {
            next if ($mode eq 'keep'   && ! defined($h_list_of_ids{$nameinsert}));
            next if ($mode eq 'remove' && defined($h_list_of_ids{$nameinsert}));
        }

        if (defined($mode) && $mode eq 'substr')
        {
            # Next if no substring positions defined
            next if (!exists($h_cols{$name}));

            $sep = &SEP;
			$h_cols{$name} =~ s/$sep$//;

            # get the positions of the substrings
            my $allsubstr = $h_cols{$name};    # "10,20_30,45"
            my @a_substr = split(/$sep/, $allsubstr);    # [10,20][30,45]
                                                         # For each substring position
            foreach my $one_substr (@a_substr)
            {
                $sep = &SEP_POS;
                my @a_pos = split(/$sep/, $one_substr);    #[10,20]
                my $fullseq = $h_seq{$name}->{sequence};
                my $seq = &ExtractSeq (\$fullseq,$a_pos[0], $a_pos[1]);
				$seq = &RevCompSeq(\$seq) if ( $a_pos[3] eq '-' );
                printf $f_outfile ">$a_pos[2]:$a_pos[0]:$a_pos[1]:$a_pos[3] len=%d\n$seq\n",$a_pos[1]-$a_pos[0]+1;
            }
        }
		elsif (defined($mode) && $mode eq 'mask')
		{
            # Next if no substring positions defined
			my $currseq = $h_seq{$name}->{sequence};
            if (exists($h_cols{$name}))
			{
				# get the positions of the substrings
				my $allsubstr = $h_cols{$name};    # "10,20_30,45"
				$sep = &SEP;
				my @a_substr = split(/$sep/, $allsubstr);    # [10,20][30,45]
															 # For each substring position
				my @a_seq = split('',$currseq);
				foreach my $one_substr (@a_substr)
				{
					$sep = &SEP_POS;
					my @a_pos = split(/$sep/, $one_substr);    #[10,20]
					for(my $i=$a_pos[0]-1;$i<=$a_pos[1]-1;$i++)
					{
						$a_seq[$i]='N';
					}
				}
				$currseq = join('',@a_seq);
			}
            print $f_outfile ">$name $desc\n$currseq\n";
		}
        else
        {
            print $f_outfile ">$name $desc\n" . $h_seq{$name}->{sequence}, "\n";
        }
    }
    $f_outfile->close;
}

sub Usage
{
    print STDERR<<END

Filter out a fasta file using various criteria/ways (sequence length, m8 file, keyword, etc.)

$0
    --help
    --in              filename|stdin	fasta formatted file
    --out             filename|stdout

[Optional]
From an intrinsic criterium

    --min_length integer
    --max_length integer
    --keyword string

From an external source
    --db_filename     filename	a file (fasta,list,tabulated) containing the list of selected ids
    --db_column       integer	the column of the ID that is in db_filename	(warning: the first column is 1)

    --keep                      flag to keep only hits related to sequences in db_filename
    or 
    --remove                    flag to remove hits related to sequences in db_filename
    --bothends        quoted string    pattern defining the insert (ex: '(.+)\\/' for solexa paired-ends (tested with --keep, not yet with --remove)

    --substr                    flag to extract the subsequences 
	or
	--mask                      flag to mask (N) the subsequences
    --cols            string    number of the columns which contain start,stop positions and the ID of the subsequence to extract. Syntax: col1,col2,col3 (warning: the first column is 1)

Compressed filename are allowed for in, out and db_filename parameters

END
}
