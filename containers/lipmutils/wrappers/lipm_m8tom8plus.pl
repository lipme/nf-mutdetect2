#!/usr/bin/perl

BEGIN
{
	$VERSION = do { my @r = (q$Rev: 1293 $ =~ /\d+/g); sprintf "%d."."%02d" x $#r, @r };
	$ENV{LANG} = "C";
}

use strict;

use Data::Dumper;
use ParamParser;
use General;
use File::Basename;
use GeneralBioinfo;
use Cwd 'abs_path';


MAIN:
{
	my $o_param = New ParamParser('GETOPTLONG', \&Usage, 'm8=s', 'm8plus=s','qfasta=s', 'sfasta=s','help');

	
	
	
	$o_param->AssertFileExists('m8');
	$o_param->AssertDefined('m8plus','qfasta','sfasta');
	
	my $ra_qfasta = $o_param->GetListOfFiles('qfasta');
	my $ra_sfasta = $o_param->GetListOfFiles('sfasta');
	
	my $rh_fasta = {};
	foreach my $sfasta (@{$ra_sfasta})
	{
		&LocalFastaToHash($sfasta, $rh_fasta);
	}
	foreach my $qfasta (@{$ra_qfasta})
	{
		&LocalFastaToHash($qfasta, $rh_fasta);
	}
	
	
	
	&Enrich($o_param->GetStreamIn('m8'),$o_param->GetStreamOut('m8plus'),$rh_fasta);
	exit;
}

sub LocalFastaToHash
{
		my ($file_in, $rh_sequences) = @_;
		
		my $fh_in = &GetStreamIn($file_in);
		my $id  = '';
		my $cpt = 0;
		while (my $line = <$fh_in>)
		{
			chomp($line);
        if ($line =~ /^>(\S+)\s*(.*)/)
        {
            $id = $1;
            my $header = (defined($2)) ? $2 : "";
            $id =~ s/lcl\|//;
            $cpt++;

            if (defined($rh_sequences->{$id}->{len}))
            {
                Carp::carp("WARNING - sequence $id already loaded\n");
            }
            $rh_sequences->{$id}->{header} = $header;
            $rh_sequences->{$id}->{db} = basename($file_in);
            while ($header =~ /(\S+)=(\S*)/g)
            {
                $rh_sequences->{$id}->{$1} = $2;
            }
            if ($header =~ /def=(.+)/)
            {
				
                $rh_sequences->{$id}->{def} = $1;
            }
            $rh_sequences->{$id}->{len} = 0;
        }
        elsif (defined $rh_sequences->{$id})
        {
            #$rh_sequences->{$id}->{sequence} .= $line;
            $rh_sequences->{$id}->{len} += length($line);
        }
		}
		$fh_in->close;
		
		return;
}

sub Enrich
{
    my ($fh_m8,$fh_m8plus,$rh_fasta) = @_;
    
    while (my $line = <$fh_m8>)
    {
		chomp $line;
		if ($line !~ /^#/)
		{
			my $rh_m8_line = {};
			&ParseM8Line(\$line, 1, $rh_m8_line);
			
			my $qhsplen = $rh_m8_line->{qend} - $rh_m8_line->{qstart} + 1;
			if ($rh_m8_line->{qend} <  $rh_m8_line->{qstart})
			{
					$qhsplen = $rh_m8_line->{qstart} - $rh_m8_line->{qend} + 1;
			}
			
			my $query = $rh_m8_line->{query};
			$query =~ s/lcl\|//;
			
			die ('ERROR - ' . $query . ' NOT FOUND IN FASTA FILES' ) if (! defined $rh_fasta->{$query});
			
			my $thsplen = $rh_m8_line->{tend} - $rh_m8_line->{tstart} + 1;
			if ($rh_m8_line->{tend} <  $rh_m8_line->{tstart})
			{
					$thsplen = $rh_m8_line->{tstart} - $rh_m8_line->{tend} + 1;
			}
			my $target = $rh_m8_line->{target};
			$target =~ s/lcl\|//;
			
			die ('ERROR - ' . $target . ' NOT FOUND IN FASTA FILES' ) if (! defined $rh_fasta->{$target});


			my $query_len = $rh_fasta->{$query}->{len};
			my $query_def = $rh_fasta->{$query}->{def};
			my $query_db = $rh_fasta->{$query}->{db};
			
			my $target_len = $rh_fasta->{$target}->{len};
			my $target_def = $rh_fasta->{$target}->{def};
			my $target_db = $rh_fasta->{$target}->{db};
			
			
			my $pcq = ($qhsplen / $query_len ) * 100;
			my $pcs = ($thsplen / $target_len ) * 100;
			
			
			print $fh_m8plus $line . "\t" . sprintf("%.2f",$pcq) . "\t" . sprintf("%.2f",$pcs) . "\t" . $query_def ."\t" . $target_def ."\t" . $query_db ."\t" . $target_db ."\n";
		}
		elsif ($line =~ /Fields:/)
		{
			print $fh_m8plus $line . ', q.percent, s.percent, q.definition, s.definition, q. databank, s. databank' . "\n";		
		}
		else
		{
			print $fh_m8plus $line	. "\n";		
		}
		
	}
	$fh_m8->close;
	$fh_m8plus->close;
	
    return;
}

sub Usage
{
	
	print STDOUT <<END;
Add Query/Subject alignment percentage, definition and databank origin to m8 formatted file
	
	$0 --m8=<m8 or m9 formatted file> --qfasta=<Query multi fasta file or list of files (coma separated)> --sfasta=<Subject multi fasta file or list of files (coma separated)> --m8plus=<Enriched m8 or m9 file>

END
    
}
