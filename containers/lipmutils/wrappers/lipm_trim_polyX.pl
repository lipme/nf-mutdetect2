#!/usr/bin/perl

use strict;
use warnings;
use General;
use GeneralBioinfo;
use ParamParser;

# svn propset svn:keywords "Id Rev" *.pl

# isoseq data (ccs+lima => oriented, anchored on 3' ==> polyA tail at the 3' end)
# extended to the 4 nt in order to facilitate gmap mapping

my $DEF_MINLEN = 150;

MAIN:
{
	my $o_param = New ParamParser('GETOPTLONG',\&Usage,'in=s','out=s','bothends','split-on-strech=i','minlen=i');

	$o_param->AssertFileExists('in');
	$o_param->AssertDefined('out');
	my $bothends = $o_param->IsDefined('bothends');
	$o_param->SetUnlessDefined('split-on-strech',0);
	my $splitnt = $o_param->Get('split-on-strech');
	$o_param->SetUnlessDefined('minlen',$DEF_MINLEN);
	my $minlen = $o_param->Get('minlen');

	my %h_seq = ();

	&FastaToHash($o_param->Get('in'),\%h_seq);

	my $f_out = $o_param->GetStreamOut('out');
	foreach my $id (sort keys(%h_seq))
	{
		my $seq = $h_seq{$id}{sequence};
		foreach my $nt ('A','T','G','C')
		{
			my $prevlen = 1;
			my $len     = 0;
			while( $prevlen != $len )
			{
				$prevlen = length($seq);
				$seq =~ s/$nt{10,}$//i;  			# removing polyX tail anchored at the sequence end
				$seq =~ s/$nt{10,}.{0,10}$//i;		# now add some tolerance (max 10nt)
				if ( defined($bothends) )
				{
					$seq =~ s/^$nt{10,}//i;  			# removing polyX head anchored at the sequence start
					$seq =~ s/^.{0,10}$nt{10,}//i;		# now add some tolerance (max 10nt)
				}
				$len = length($seq);
			}
		}
		my @a_slices = ();
		if ( $splitnt )
		{
			$seq=~s/A{$splitnt,}/N/ig;
			$seq=~s/T{$splitnt,}/N/ig;
			$seq=~s/C{$splitnt,}/N/ig;
			$seq=~s/G{$splitnt,}/N/ig;
			my @a_slices = split(/N+/,$seq);
			if ( $#a_slices > 0 )
			{
				my $cpt=1;
				my $n = length($#a_slices); # compute digit number required to display all slices
				foreach my $slice (@a_slices)
				{
					my $len = length($slice);
					next if ( $len < $minlen );
					printf $f_out ">$id.%0${n}d len=$len\n",$cpt;
					print  $f_out "$slice\n";
					$cpt++;
				}
				next;
			}
		}
		my $len = length($seq);
		next if ( $len < $minlen );
		print $f_out ">$id len=$len\n";
		print $f_out "$seq\n";
	}
	$f_out->close;
	exit 0;
}

sub Usage
{
print STDERR<<END
$0 - remove 3' N streches
	--help
[Mandatory]
	--in                           fastafile
	--out                          fastafile
[Optional]
	--bothends                     remove 5' N streches too
	--split-on-strech integer      split sequences containing internal <integer>-nucleotide strech
	--minlen                       [$DEF_MINLEN] sequence length cut-off
END
}

