#!/usr/bin/perl

#$Id

# Copyright INRA

# Jerome.Gouzy@toulouse.inra.fr
# Erika.Sallet@toulouse.inra.fr
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



BEGIN
{
 	$ENV{LANG} = "C";
 	our $VERSION = do {my @r = (q$Rev: 747 $ =~ /\d+/g); $r[0]};
}


use strict;
use warnings;
use ParamParser;
use Carp;
use General;
use GeneralBioinfo;
use RnaSeq;


my $SCALINGFACTOR = 1;
my $SCOREMAX      = &MAXINT;
my $MIN_EXON_LEN  = 18;
my $MIN_COVERAGE =  1;
MAIN:
{
    my $o_param = New ParamParser('GETOPTLONG', \&Usage, "wigfile=s","outprefix=s", 
                                  "strand=s", 
								  "scalingfactor=i","scoremax=i", "tag=s",
								  "min_mean_coverage=i","min_coverage=s", 
								  "min_coverage_percent=i",
								  "min_exon_len=i", "verbose", "split");
	
	$o_param->SetUnlessDefined('tag',"SRCTAG");

	$o_param->AssertFullPath('outprefix');
	my $outprefix = $o_param->Get('outprefix');
	my ($outpath, $outname) = $outprefix =~ /(.+)\/(.+)/;
	$o_param->Set('outpath', $outpath);
	$o_param->Set('outname', $outname);

	$o_param->AssertDefined('strand');
	$o_param->AssertAllowedValue('strand', "forward", "reverse");
	$o_param->AssertFileExists('wigfile');
	$o_param->SetUnlessDefined('scalingfactor', $SCALINGFACTOR);
	$o_param->SetUnlessDefined('scoremax',      &MAXINT);
	$o_param->SetUnlessDefined('min_coverage', $MIN_COVERAGE);	
	$o_param->SetUnlessDefined('min_exon_len', $MIN_EXON_LEN);

	
	&Depth($o_param);
	
}



sub LoadCov
{
	my($o_param, $rh_, $ra_seq) = @_;

	my $seqname       = "";
	my $well_formated = 0;
	my $track_type    = "";
	my $span;	
	my $step;
	my $current_pos;

	my $f_wig = $o_param->GetStreamIn('wigfile');
	while(my $l = <$f_wig>)
	{
		chomp($l);
		$l =~ s/^\s+//;
		$l =~ s/\s+$//;

		if ($l =~ /^track/)
		{
			# New track
			$well_formated = 1;
			next;
		}
		if ($l =~ /^variableStep/)
		{
			# Assert line data
			if ($l !~ /chrom=/)
			{
				Carp::croak("Line variableStep has to contain \"chrom=\" value. Check your wig file\n");
			}
			#print "variable";
			$well_formated = ($well_formated == 1) ? 2 : 0;
			$track_type    = "variable";
			($seqname)     = $l =~ /chrom=(\S+)/;
			$span          = ($l =~ /span=(\d+)/) ? $1 : 1;
			push(@$ra_seq, $seqname);
			next;
		}
		
		if ($l =~ /^fixedStep/)
		{
			if ($l !~ /chrom=/ || $l !~ /start=\d+/)
			{
				Carp::croak("Line fixedStep has to contain \"chrom=\" and \"start=\" values. Check your wig file\n");
			}
			
			$well_formated = ($well_formated == 1) ? 2 : 0;
			$track_type    = "fixed";
			($seqname)     = $l =~ /chrom=(\S+)/;
			($current_pos) = $l =~ /start=(\d+)/;
			$span          = ($l =~ /span=(\d+)/) ? $1 : 1;
			$step          = ($l =~ /step=(\d+)/) ? $1 : 1;
			push(@$ra_seq, $seqname);
		}
		
		Carp::croak("File have to start by track line then fixed|variableStep line. Check your wig file.\n") if ($well_formated != 2);

		if ($l =~ /^(\d+)\s+(-?[\d\.]+)$/ && $track_type eq "variable")
		{
			my $pos = $1;
			my $val = $2;
			&UpdateCov($pos, $pos+$span-1, $val, $rh_, $seqname);
		}
		elsif ($l =~ /^(-?[\d\.]+)$/ && $track_type eq "fixed")
		{
			my $val = $1;

			&UpdateCov($current_pos, $current_pos+$span-1, $val, $rh_, $seqname) if ($val != 0);
			$current_pos += $step;	
		}
		else
		{
			Carp::croak("Not seem to be a wig file. Check the format\n") if ($well_formated != 2);
		}
	}
	$f_wig->close;

	return;
}


sub Depth
{
	my($o_param) = @_;

	my %h_cov; # key=seqname, value = ref to an array of cov
	my @a_seq; # List of sequename

	&LoadCov($o_param, \%h_cov, \@a_seq);

	my $so                = "transcript_region";	
	my $cpt               = 1;
	my $scalef            = $o_param->Get('scalingfactor');
	my $maxvalue          = $o_param->Get('scoremax');	
	my $tag               = $o_param->Get('tag');
	my $min_mean_coverage = $o_param->Get('min_mean_coverage') if ($o_param->IsDefined('min_mean_coverage'));     
	my $min_exon_len      = $o_param->Get('min_exon_len');
	my $strand            = ($o_param->Get('strand') eq "forward") ? "+" : "-";
	my $min_coverage      = $o_param->Get('min_coverage');;
	my $f_position_out;
	my $f_region_out;
	
	if (!$o_param->IsDefined('split'))
	{
		$f_position_out = &GetStreamOut($o_param->Get('outprefix').".position.gff3");
		$f_region_out   = &GetStreamOut($o_param->Get('outprefix').".region.gff3");
	}

	foreach my $seqname (@a_seq)
	{
		my $ra_seqdata = $h_cov{$seqname}; # cov for this sequence
		my $maxi       = scalar(@$ra_seqdata) -1; 
		
		if ($o_param->IsDefined('min_coverage_percent'))
		{
			my @a_definedseqdata;
			for(my $i=1; $i<=$maxi; $i++)
			{
				push(@a_definedseqdata, $$ra_seqdata[$i]) if (defined($$ra_seqdata[$i]));
			}
			$min_coverage = &RnaSeq::GetMinCoverage(\@a_definedseqdata, $o_param->Get('min_coverage_percent'));
		}		
        # => at least two overlapping hits if min_coverage == 1
		$min_mean_coverage = ($min_coverage+$min_coverage/10) if (!$o_param->IsDefined('min_mean_coverage'));
		

		if ($o_param->IsDefined('split'))
		{
			$f_position_out = &GetStreamOut($o_param->Get('outpath')."/$seqname.".$o_param->Get('outname').".position.gff3");
			$f_region_out   = &GetStreamOut($o_param->Get('outpath')."/$seqname.".$o_param->Get('outname').".region.gff3");
		}
		
		for(my $i=1; $i<=$maxi; $i++)
		{
			if ( defined($$ra_seqdata[$i]) && $$ra_seqdata[$i] >= $min_coverage )
			{
				&PrintGFFLine($seqname, $so, $f_position_out, 'depth', $tag, $i, $i, 
							  &Min( $$ra_seqdata[$i]/$scalef,$maxvalue),$strand, $cpt++);
			}
		}
		
		if ( scalar(@$ra_seqdata) > 0 )
		{
			&DumpRegionsGFF3($seqname, $so, $f_region_out, 'region', $tag, $min_coverage, $min_mean_coverage,
							 $min_exon_len, $scalef, $maxvalue, $strand, $ra_seqdata);
		}

		if ($o_param->IsDefined('split'))
		{
			$f_position_out->close;
			$f_region_out->close;
		}	
	}

	if (!$o_param->IsDefined('split'))
	{
		$f_position_out->close;
		$f_region_out->close;
	}

	return;
}


sub DumpRegionsGFF3
{
	my ($seqname, $so, $f_out, $type, $tag, $min_coverage, $min_mean_coverage,
	$min_exon_len, $scalef, $scoremax, $strand, $ra_seq_cov) = @_;

	my $cpt          = 1;
	my $region_start = 0;
	my $mean_cov     = 0;

	for(my $i=1;$i<=($#$ra_seq_cov+1);$i++) # to be sure to test the last undef and then print out the last hit
	{
		if ( $region_start )
		{
			# end of a region! We write the region and reinit the variables
			if ( !defined($$ra_seq_cov[$i]) || $$ra_seq_cov[$i] < $min_coverage)
			{
				my $len = ($i-1) - $region_start + 1; 
				#my $mean_cov_raw = $mean_cov;
				$mean_cov = $mean_cov/$len;
				if ( $len >= $min_exon_len && $mean_cov >= $min_mean_coverage )
				{	
					my $score = &Min($mean_cov/$scalef,$scoremax);
					if ( $so eq 'intergenic_region' )
					{
						$score *= -1;
					}
					if ( $strand eq '.' )	
					{
						&PrintGFFLine($seqname,$so,$f_out,$type,$tag,$region_start,$i-1,$score/2,'+',$cpt++);
						&PrintGFFLine($seqname,$so,$f_out,$type,$tag,$region_start,$i-1,$score/2,'-',$cpt++);
					}
					else
					{
						&PrintGFFLine($seqname,$so,$f_out,$type,$tag,$region_start,$i-1,$score,$strand,$cpt++);
					}
				}
				
				$cpt++;
				$region_start = 0;
				next;
			}
		}
		elsif ( defined($$ra_seq_cov[$i]) && $$ra_seq_cov[$i] ne '' && $$ra_seq_cov[$i] >= $min_coverage)
		{
			# Start of a region
			$region_start = $i;
			$mean_cov     = 0;
		}

		if ( defined($$ra_seq_cov[$i]) && $$ra_seq_cov[$i] ne '' )
		{
			$mean_cov += $$ra_seq_cov[$i];
		}
	}

	return;
}



sub PrintGFFLine
{
	my($seqname, $so, $f_out, $type, $tag, $start, $end, $score, $strand, $cpt) = @_;

	my ( $code  ) = ( $strand eq '.' ) ? "u" : (( $strand eq '+' ) ? 'p' : 'm' );
	printf $f_out "$seqname\t$type$tag\t$so\t$start\t$end\t%.4f\t$strand\t.\tID=$seqname:$type$tag:$code$start.$cpt\n",$score;

	return;
}



sub UpdateCov
{
	my($min, $max, $count, $rh_, $seqname) = @_;
	
	my $ra_seq;
	if (defined $$rh_{$seqname})
	{
		#print "yess : $seqname\n";
		$ra_seq = $$rh_{$seqname};
	}
	else
	{
		my @a_data = ();
		$$rh_{$seqname} = \@a_data;
		$ra_seq = $$rh_{$seqname};
	}
	
	for (my $i=$min; $i<=$max; $i++)
	{
		$$ra_seq[$i] += $count;
	}

	return;
}



sub Usage
{
print STDERR<<END
$0 
    --help
    --verbose
    --wigfile   filename
    --outprefix filename
    --strand    forward|reverse

[Optional]
    --scalingfactor        integer [$SCALINGFACTOR]
    --scoremax             integer [$SCOREMAX]	    (e.g set to $SCOREMAX, all scores > $SCOREMAX (after rescaling) will be set to 1)
    --tag                  string  [SRCTAG] Value of the column 3 in GFF3 file results
    --min_coverage         integer [1] Minimum coverage at each position 
    --min_coverage_percent integer Minimum percentage of reads to keep [If defined, disable min_coverage parameter]  
    --min_mean_coverage    integer Minimum mean coverage of regions
    --min_region_len       integer	[$MIN_EXON_LEN] Minimum region length
    --split                     Split by sequence in different files.

	
END
}
