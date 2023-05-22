#!/usr/bin/perl

use strict;
use General;

my $mpileup = $ARGV[0];
my $pos = $ARGV[1];


my $fh_pos = &GetStreamIn($pos);

my %h_pos = ();
while (my $l = <$fh_pos>)
{
	chomp $l;
	$l=~ s/\t/:/;
	$h_pos{$l} = 1;
}
$fh_pos->close();

my $fh_mpileup = &GetStreamIn($mpileup);
while (my $l = <$fh_mpileup>)
{
	chomp $l;
	my ($id) = ($l =~ /^(\S+\t\d+)/);
	$id=~ s/\t/:/;
	if (defined $h_pos{$id})
	{
		print "$l\n";
		delete $h_pos{$id};
	}
}

$fh_mpileup->close();

exit 0;
