package RnaSeq;

use strict;
use warnings;

# Copyright INRA

# Erika.Sallet@toulouse.inra.fr
# Jerome.Gouzy@toulouse.inra.fr
# Thomas.Schiex@toulouse.inra.fr

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

=head1 NAME

=head1 DESCRIPTION

=cut
BEGIN
{
    use Exporter ();
    use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

    @ISA = qw(Exporter);
    @EXPORT = qw( GetMinCoverage);
}

=head2 function GetMinCoverage

	Title      : GetMinCoverage
	Usage      : &GetMinCoverage(\@a_wig,90)	
	Prerequiste: 	
	Function   : Return the minimum coverage value to keep at least $percentage % 
                 of the total count of the data
	Returns    : The minimum coverage value
	Args       : $ra_data: ref to an array of values. each value is the expression depth 
                 at a given position
	             $percentage: percentage of data to keep
	globals    : none

=cut

sub GetMinCoverage
{
	my ($ra_data, $percentage) = @_;
	
	my $totalcount = 0; # count the total of the expression
	my $posnb      = 0;

	foreach my $val (@$ra_data)
	{ 
		$totalcount += $val;
		$posnb++; # Number of position
	}
	
	my $threshold_count  = int($totalcount*$percentage/100);
	my $threshold_pos_nb = int($posnb*$percentage/100);

	# descending sort of the expressions
	my @a_sorted_data = sort { $b <=> $a } @$ra_data;
	my $data_len      = scalar(@a_sorted_data);

	my $count = 0;
	my $i     = 0;
	# Count the expression while the threshold is not reached
	for ($i = 0; $count < $threshold_count && $i < $data_len; $i++)
	{
		$count += $a_sorted_data[$i];
	}
	# Here we cumulate at least $percentage % of the total expression:
	# Get the last expression value
	my $min_coverage =  $a_sorted_data[$i];

	my $realcount  = 0;
	my $othercount = 0;
	foreach my $v (@a_sorted_data)
	{
		$othercount += $v if ( $v >  $min_coverage );
		$realcount  += $v if ( $v >= $min_coverage );
	}
	
	#print "Min Coverage=$min_coverage for $percentage %. ($realcount/$totalcount. If min_coverage = ($min_coverage+1): $othercount)\n";
	#print "[Min value to get $percentage % of the position=" . $a_sorted_data[($threshold_pos_nb+1)] . "]\n";

	return $min_coverage;
}


1;    # package RnaSeq




