#!/usr/bin/perl

# ISAblatFilter.pm
# Dr. Xiaolin Wu & Daria W. Wells

# This module contains the subroutines that process the demultiplexed and pair-matched 
# ISA data and generates fasta files for alignment.

use 5.18.4;

package ISAblatFilter; 

use strict;
use warnings;
use Exporter qw(import);

our @EXPORT_OK = qw(blat_filter1 blat_filter2 blat_filter3 blat_filter4 blat_filter5 blat_filter6 blat_filter7);


#Blat output format
#psLayout version 3
#
#match	mis- 	rep. 	N's	Q gap	Q gap	T gap	T gap	strand	Q        	Q   	Q    	Q  	T        	T   	T    	T  	block	blockSizes 	qStarts	 tStarts
#     	match	match	   	count	bases	count	bases	      	name     	size	start	end	name     	size	start	end	count
#---------------------------------------------------------------------------------------------------------------------------------------------------------------
#80	0	0	0	0	0	0	0	+	FTF2IS401C9OTF_LTR	123	42	122	chr1	197195432	43516678	43516758	1	80,	42,	43516678,
#71	0	0	0	0	0	0	0	-	FTF2IS401A43VI_LTR	77	0	71	chr1	197195432	82334892	82334963	1	71,	6,	82334892,


# This subroutine filters out alignments from the BLAT output that start more than 3 bases
# into the ISA sequence and are shorter than 20 bases.
sub blat_filter1 {
	# Create the input and output file names.
	my ($read, $folder) = @_;
	my $file2open = "$folder/$read\_blat_results.txt";
	my $file2out = "$folder/$read\_blat_results_high_quality.txt";

	my @lines; # The passing alignments will be stored in @lines.

	open my $data_in, '<', $file2open or die "Unable to open $file2open in BLAT filter 1: $!\n";

	my $num_passed_reads = my $num_rejected_reads = 0;
	# Test each line.
	while (my $line=<$data_in>) {
		next unless $line =~ /^\d/;  # Skip lines that don't start with a digit.
		my @eles=split(/\t/, $line);
		#start <4
		#Qgap[5] <10
		#Tgap[7] <10
		#identity >95% [0]/([12]-[11])
		if($eles[11]<4 and $eles[5]<10 and $eles[7]<10 and $eles[0]>19 and $eles[0]/($eles[12]-$eles[11])>0.95) {
			push(@lines, $line);
			$num_passed_reads++;
		}else{$num_rejected_reads++;}
	}
	close $data_in;

	open my $data_out, ">$file2out" or die "Unable to print data in BLAT filter 1: $!\n";

	# Sort @lines to organize the matches for the next filter step.
	my @sorted_lines = sort {
		my @a_fields = split /\t/, $a;
		my @b_fields = split /\t/, $b;

		$a_fields[9] cmp $b_fields[9]  # sort by seq ID
		||
		$b_fields[0] <=> $a_fields[0]  # then sort by match length
	} @lines;
	
	foreach my $line (@sorted_lines) {
		print $data_out "$line";
	}
	close $data_out;
}


##########################################################################################
# This subroutine iterates through the sorted data from sub blat_filter1 to determine which
# integration read pairs can be mapped to a unique genomic location (single-mappers) and which 
# map roughly as well to more than one location (multi-mappers).

sub blat_filter2 {
	$DB::single=2;
	my ($folder) = @_;
	
	# Create the input file names
	my $read1_in_file = "$folder/read1_blat_results_high_quality.txt";
	my $read2_in_file = "$folder/read2_blat_results_high_quality.txt";
		
	# Create the output file names.
	my $read1_out_best_file = "$folder/read1_blat_results_best_alignments.txt";
	my $read1_out_multi_file = "$folder/read1_blat_results_repetitive_regions.txt";
	my $read2_out_best_file = "$folder/read2_blat_results_best_alignments.txt";
	my $read2_out_multi_file = "$folder/read2_blat_results_repetitive_regions.txt";
	
	# Create hashes to store the single- and multi-mapping reads. Parse the data.
	my ($best_read1, $repeat_read1) = parse_alignments($read1_in_file);
	my ($best_read2, $repeat_read2) = parse_alignments($read2_in_file);
	
	# In some instances, alignments to the same integration site will be considered a single mapper if the read
	# is long enough but a multi-mapper if it's short. Change these multis to singles if the integration
	# sites are identical and the length of the longer read is at least 2x + 20 bases longer than
	# the shorter reads.
	($best_read1, $repeat_read1) = compare_read_lengths($best_read1, $repeat_read1);
	($best_read2, $repeat_read2) = compare_read_lengths($best_read2, $repeat_read2);
	
	# In some instances, one read will map only once while the other maps multiple times. Change the multi-mapper
	# to a single-mapper if the other read appears in the alignment data only once.
	($best_read1, $repeat_read1, $best_read2, $repeat_read2) = crosscheck_for_singles($best_read1, $repeat_read1, $best_read2, $repeat_read2);
	
	
	# Create hashes to store the best_reads, but organized by integration site.
	# my ($read1_best_by_IS, $read2_best_by_IS) = best_to_by_IS($best_read1, $best_read2);
	
	# Print the results
	print_best($best_read1, $read1_out_best_file);
	print_best($best_read2, $read2_out_best_file);
	print_repeat($repeat_read1, $read1_out_multi_file);
	print_repeat($repeat_read2, $read2_out_multi_file);

	
	return;
}


##########################################################################################

# This subroutine reads the filtered BLAT data and converts it into the preferred 
# results format. The read1 and read2 paired-end read coordinates go into their respective 
# hashes with the sequence ID as the key and are printed at the end to two files. One for further
# filtering and another for Uclust. Read1 sequences without corresponding read2's are 
# included in the output for Uclust.

sub blat_filter3 {
	$DB::single=2;
	my $folder = shift;
	
	my $read1_file = "$folder/read1_blat_results_best_alignments.txt";
	my $read2_file = "$folder/read2_blat_results_best_alignments.txt";
	open my $in_L, '<', $read1_file or die "Unable to open $read1_file in blat_filter3: $!\n";
	open my $in_R, '<', $read2_file or die "Unable to open $read1_file in blat_filter3: $!\n";
	
	my %ID2corL; # Hash to store the left reads.
	my %ID2corR; # Hash to store the right reads.
	my %id2count; # Hash to count how many times a set of coordinates is seen.
	
	# @parameters consists of two array references. Each one stores the left/right input
	# file handle and a reference to the read hash. 
	my @parameters = ([$in_L, \%ID2corL],[$in_R, \%ID2corR]);
	
		foreach my $p (@parameters) { # Perform this analysis for each array reference in @parameters.
			my $fh = $p->[0]; # $p->[0] is the file handle to read from.
			while (my $line = <$fh>) { # Iterate through the input file
				next if $line !~ m/\w/;
				chomp($line);
				my @elements = split(/\t/, $line); # Make an array of each tab-separated line element.
				my $id = $elements[9]; # The sequence ID with the unique pair count.
				my ($idwocount, $MID, $count) = (split /\#/,$id); # Separate the ID from the count.
				$id = $idwocount.'#'.$MID;
				$id2count{$id} = $count; # Store the count.
			
				my $chrcor; # The sequence coordinates.
				if($elements[8] eq "+") { # The line elements for the coordinates depend on the strand.
					$chrcor="$elements[13]"."\t+\t"."$elements[15]"; # Example: chr17	+	3618824
				}else{
					$chrcor="$elements[13]"."\t-\t"."$elements[16]";
				}
			
				$p->[1]{$id} = $chrcor; # Store the coordinates. 
			}	
		}
	close $in_L;
	close $in_R;
		$DB::single=2;
	my @bestmatches; # Create a list of the read 1 best match sequence IDs and coordinates.
	# Iterate through the read 1 hash
	while (my ($id, $cor)=each(%ID2corL)) {
		my $bestmatch = "$id\t$cor"; # Create the first portion of the output line 
		push (@bestmatches, $bestmatch); # Add the portion to @bestmatches
	}
	
	# Sort @bestmatches.
	my @sorted_bestmatches = sort{
	# Split the string into two parts
	my @a_fields = split /\t/, $a;
	my @b_fields = split /\t/, $b;
		
	$a_fields[1] cmp $b_fields[1]  # First sort by chromosome.
	||
	$a_fields[3] <=> $b_fields[3]  # Then by ascending nucleotide coordinate.
	} @bestmatches;
	
	# Print the coordinates for every left-hand read along with the corresponding
	# right-hand read if there is one.
	# Example 9	M02560:63:000000000-ACWWU:1:1101:17574:2571#9	chr17	+	3618824	chr17	-	3618893
	my $matched_pair_file = "$folder/combined_unfiltered_match_pairs.txt";
	my $uclust_file = "$folder/all_read1_coordinates_for_Uclust.txt";
	open my $matched_out, '>', $matched_pair_file or die "Cannot open $matched_pair_file for printing in blat_filter3: $!\n";
	open my $uclust_out, '>', $uclust_file or die "Cannot open $uclust_file for printing in blat_filter3: $!\n";
	
	# Print the sorted data
	foreach my $RIS (@sorted_bestmatches) {
		
		my @elements=split /\t/, $RIS; # Split the line into its component parts
		my $id = $elements[0]; # Obtain the sequence ID
		# If this sequence ID has two coordinates, print it to both files. If it has only one,
		# only print it to the Uclust file
		if (exists $ID2corR{$id}) {
			print $matched_out "$id2count{$id}\t$RIS\t$ID2corR{$id}\n";
			print $uclust_out "$id2count{$id}\t$RIS\t$ID2corR{$id}\n";
		}else{
			print $uclust_out "$id2count{$id}\t$RIS\n";
		}
	}

	close $matched_out;
	close $uclust_out;
	
	# Open repeat match to print multi-map meta data to the log file
	my $repeat_read1 = "$folder/read1_blat_results_repetitive_regions.txt";
	my $repeat_read2 = "$folder/read2_blat_results_repetitive_regions.txt";
	my %seen_repeats;
	# Open each file handle separately. Count the number of unique sequence IDs seen
	foreach ($repeat_read1, $repeat_read2) {
		open my $handle, "<", $_ or die "Unable to open $_ in blat_filter3: $!\n";
		while (my $line = <$handle>) {
			(my $id) = (split /\t/, $line)[10];
			$seen_repeats{$id} = "" unless exists $seen_repeats{$id};
		}
		close $handle;
	}
	my $repeat_count = keys %seen_repeats;
	# Return the number of unique sequence IDs for use in sub blat_filter4
	return($repeat_count);
}

##########################################################################################

# This subroutine removes integration sites from the data if the paired-reads aren't on the
# same chromosome, are on the same strand, or are farther than 10,000 bases away from
# each other.

sub blat_filter4 {
	$DB::single=2;
	my ($folder, $repeat_count) = @_;
	my $input_file = "$folder/combined_unfiltered_match_pairs.txt";
	my $output_file = "$folder/filtered_out_unmatched_read_loci.txt";
	open my $data_in, '<', $input_file or die "Unable to open $input_file in blat_filter4: $!\n";
	open my $data_out, '>', $output_file or die "Unable to open $output_file for printing in blat_filter4: $!";
	my $total_unique_pairs = my $total_raw_pairs = my $all_pairs = 0;
	
	# Iterate through the data
	while (my $line=<$data_in>) {
		next unless $line =~ m/\w/; # Skip any blank lines, like the last one.
		chomp($line);
		$all_pairs++; # Count ever line analyzed
		my @eles = split /\t/, $line; # Split the line into its component parts
		if( # If the two coordinates meet the criteria, print and count them. Otherwise,
			# exclude them from the next report.
			($eles[2] eq $eles[5]) and 
			($eles[3] ne $eles[6]) and
			(abs($eles[4]-$eles[7]) <10000)
		) {
			print $data_out "$line\n";
			# Count the number of pairs that passed and the number of raw reads
			# associated with each passing pair
			$total_unique_pairs += 1;
			$total_raw_pairs += $eles[0];
		}
			
			
	}

	close $data_in,;
	close $data_out;
	
	open my $log_fh, ">>$folder/blat_log.txt" or die "Unable to print to log file in blat_filter5 $!\n";
	print $log_fh "Filter unmatched read1 and read2 coordinates\nPairs analyzed:\t$all_pairs\n\tUnique pairs mapped to 1 location:\t$total_unique_pairs\n";
	print $log_fh "\tRaw pairs mapped to 1 location:\t$total_raw_pairs\n\tRaw pairs with at least one read\n\t\tmapping to multiple locations:\t$repeat_count\n\n";
	close $log_fh;
	
	return;
}

##########################################################################################

# This subroutine determines which reads should be checked for potential merging.
# Integration sites within 10 bases are sent to sub merged_lines for further review.

sub blat_filter5 {
	$DB::single=2;
	my $folder = shift;

	my $input_file =  "$folder/filtered_out_unmatched_read_loci.txt";
	my $output_file = "$folder/merged_integration_sites.txt";
	open my $fh, '<', $input_file or die "Unable to open $input_file in blat_filter5: $!\n";

	# Read all the lines. The lines are sorted by integration site.
	my @data_lines = <$fh>;
	my $mapped_unique_pairs = @data_lines;
	
	my $line = shift @data_lines; # Define the first line to be compared
	chomp($line);
	 
	# Split the line into individual scalars. The integration site coordinates will be compared.
	my ($junction_chr, $junction_base) = (split "\t", $line)[2,4];
	
	my @lines_to_merge; # Seqs that may potentially be merged will be stored here.
	my @merged_data; # The final merged lines will be stored here.
	
	push @lines_to_merge, $line; # Always add the first instance of an integration site to @lines_to_merge.
	
	# Iterate through the data, sending lines that may be merged to subroutine merge_lines.
	foreach my $next_line (@data_lines) {
		chomp($next_line);
	
		# First split the next line into individual scalars
		my ($next_junction_chr, $next_junction_base) = (split "\t", $next_line)[2,4];
		
		# Sequences with integration sites on the same chromosome and within 10 bases of each
		# other may potentially be merged. Store them in @lines_to_merge for now.
		if ($junction_chr eq $next_junction_chr and abs($junction_base - $next_junction_base) <= 10) {
			push @lines_to_merge, $next_line;
			if ($next_line eq $data_lines[-1]) {
				push @merged_data, merge_lines(\@lines_to_merge);
			}	
	 	}else { # Move on to the next integration site if the integration sites aren't near each other.
			# If there is only one instance of an integration site, add it to the final data.
			# If not, send @lines_to_merge to &merged_lines.
			@lines_to_merge == 1 ? push @merged_data, @lines_to_merge : push @merged_data, merge_lines(\@lines_to_merge);
			
			if ($next_line eq $data_lines[-1]) {
				push @merged_data, $next_line;
			}
				
			# Reset the variables for the next comparisons
			($junction_chr, $junction_base) = (split "\t", $next_line)[2,4];
			$line = $next_line;
			@lines_to_merge = ($line);
		}	
	}
	
	# Separate the MID from the sequence ID and append the fragment length to the end 
	# of each line of data. The fragment length is the breakpoint coordinate minus the 
	# integration site coordinate. Print the results to a tab-delimited text file.
	# Tally the final number of read pairs.
	my $final_pair_tally = my $merged_unique_pairs = my $HIV_internal_pairs = my $host_pairs = 0;
	
	open my $out, '>', $output_file or die "Unable to open $output_file for printing in blat_filter5: $!\n";
	# Iterate through the post-merging data
	foreach my $line (@merged_data) {
		$merged_unique_pairs++; # Count the number of entries after merging
	 	my @elements = split "\t", $line; # Split each line into its component parts
	 	$final_pair_tally += $elements[0]; # Count the number of raw pairs after merging. None should be lost
	 	$elements[2] =~ m/HIV/ ? $HIV_internal_pairs += $elements[0] : $host_pairs += $elements[0];
	 	# Remove the MID from the sequence ID for the output. It's now one of the output line elements
	 	my ($sequence_id, $mid) = ($elements[1] =~ m/(^.*?)\#([ATGCN]{10})/);
	 	# Append the fuzz designation back to the sequence ID. It will be handled later
	 	$sequence_id .= "#$1" if $elements[1] =~ m/(fuzz\?{0,1})/;
	 	my $line_out = "$elements[0]\t$sequence_id\t$mid\t".join("\t", @elements[2..7])."\t".abs($elements[7] - $elements[4])."\n";
		print $out $line_out;
	}
	close $out;
	
	 
	open my $log_fh, ">>$folder/blat_log.txt" or die "Unable to print to log file in blat_filter6 merge sites. $!\n";
	print $log_fh "Merge mapped unique pairs\nUnique pairs before merging:\t$mapped_unique_pairs\n\tUnique pairs after merging:\t$merged_unique_pairs\n\tRaw pairs mapped:\t$final_pair_tally\n";
	print $log_fh "\t\tInternal HIV pairs:\t$HIV_internal_pairs\n\t\tHost integration site pairs:\t$host_pairs\n\n";
	close $log_fh;
}
	
##########################################################################################
# This subroutine merges integration sites determined to have come from the same read.	
# Sites are merged based if they meet all of the following criteria:
#	1.	The breakpoints are within 10 bases of each other.
#	2.	The MIDs differ by no more than two bases. (Hamming distance)
#	3.	The site with the larger count is greater than double the smaller count plus 1. (count score)
# Sites are then labelled as fuzz or potential fuzz if:
#	1.	If the breakpoints are within 5 base pairs of one another AND the smaller read count is 5% or 
#		less of the larger count, label it "fuzz"
#	2.	If the breakpoints are within 5 base pairs of one another AND the smaller read count is between
#		5% and 10% of the larger read count, label it "fuzz?"
#
	
sub merge_lines {
	my $lines_to_merge = shift;

	my %all_lines; # Unsorted lines wil go here
	my @sorted_lines; # Will hold the lines sorted by descending read count
	my @merged_lines; # Will hold the the merged lines from each pass
	my %already_merged; # Will indicate which lines have been merged on a previous pass.
	my @returned_lines; # The final array that will be returned
	my $changed_list = undef; # 
	
	# Each line will be stored as an array in %all_lines with $current_seq_id as the key
	foreach my $line (@{$lines_to_merge}) {
		my @split_line = split "\t", $line;
		$all_lines{$split_line[1]} = [@split_line];
	}

	# Sort the lines by descending read count.
	@sorted_lines = sort{
		$all_lines{$b}->[0] <=> $all_lines{$a}->[0]
	}keys %all_lines;


	# Merge the integration sites if all the following criteria are 1:
			#	1.	The breakpoints are within 10 bases of each other.
			#	2.	The MIDs differ by no more than two bases. (Hamming distance)
			#	3.	The site with the larger count is greater than double the smaller count plus 1. (count score)
			
	my $last_seq_ID = $sorted_lines[-1]; # Used to indicate the final array element.
		
	# Begin searching for matches.
	while (@sorted_lines != 0) {	
		
		# Compare the first (most reads) array element to the rest of the array.
		my $current_seq_id = shift @sorted_lines; # Remove the first element
		
		# Skip to the next loop if $current_seq_id has already been merged.
		next if exists $already_merged{$current_seq_id};
		
		# If $current_seq_id is the last element of the array AND it has already been
		# merged with a previous line, end the loop. If not, add it to the results.
		if ($current_seq_id eq $last_seq_ID) {
			exists $already_merged{$current_seq_id} ? last : push @merged_lines, $all_lines{$current_seq_id};
		}
	
		# Compare $current_seq_id to the rest of the sorted lines.
	 	foreach my $next_seq_id (@sorted_lines) {

			# First, check if the next array element has already been merged.
		 	if (exists $already_merged{$next_seq_id}) {
		 	
		 		# If the next element is the last element of the array and $current_seq_id
		 		# has already been merged in a previous loop, add current line to the results
		 		# array. If this is not the last element of the array, move on to the next loop.
		 		$next_seq_id eq $last_seq_ID ? (push @merged_lines, $all_lines{$current_seq_id} and next) : next;
		 	}
		 	
		 	# Begin the comparison.
	 		(my $mid) = ($all_lines{$current_seq_id}[1] =~ m/#([ATGCN]{10})$/); # Obtain the MID from the current sequence ID
	 		(my $next_mid) = ($all_lines{$next_seq_id}[1] =~ m/#([ATGCN]{10})$/); # Obtain the MID from the next sequence ID

	 		my $break_diff = abs($all_lines{$current_seq_id}[7] - $all_lines{$next_seq_id}[7]); # Calculate the distance between the two breakpoints.
	 		my $hamming = ($mid ^ $next_mid) =~ tr/\0//c; # Calculate the Hamming distance between the two MIDs.
	 		my $score = 2 * $all_lines{$next_seq_id}[0] + 1; # Calculate the count score of the comparison.
		
			# Merge and store the data_lines if all three criteria are met
			if (
				$break_diff == 0
				or
				$mid eq $next_mid
				or
				($break_diff == 1 and $all_lines{$next_seq_id}[0]/$all_lines{$current_seq_id}[0] < 0.15)
				or
				($break_diff <= 4 and $all_lines{$next_seq_id}[0]/$all_lines{$current_seq_id}[0] < 0.02)
				or
				($break_diff <= 3 and $all_lines{$current_seq_id}[0] <=5 and $all_lines{$next_seq_id}[0] <=5)
				or 
				($break_diff <= 10 and $hamming <= 2 and $all_lines{$current_seq_id}[0] > $score)
			) {
				
				$all_lines{$current_seq_id}[0] += $all_lines{$next_seq_id}[0]; # Add the smaller count to the larger count.
				$already_merged{$next_seq_id} = ""; # Indicate that the smaller count line has been merged.
				
				# This was the last element, add the current line to the results array and
				# indicate that it has been merged.
				if ($next_seq_id eq $last_seq_ID) {
					push @merged_lines, $all_lines{$current_seq_id};
					$already_merged{$current_seq_id} = "";
				}
			}elsif($next_seq_id eq $last_seq_ID) {
				push @merged_lines, $all_lines{$current_seq_id} unless exists $already_merged{$current_seq_id};
				$already_merged{$current_seq_id} = "";
			}
		}
	}
	
	# After all lines are merged that can be, check to see if any of the remaining reads are
	# considered fuzz.
	# Sites are then labelled as fuzz or potential fuzz if:
	#	1.	If the breakpoints are within 5 base pairs of one another AND the smaller read count is 5% or 
	#		less of the larger count, label it "fuzz"
	#	2.	If the breakpoints are within 5 base pairs of one another AND the smaller read count is between
	#		5% and 10% of the larger read count, label it "fuzz?"
	while (my $line = shift @merged_lines) {
		push @returned_lines, join("\t", @{$line});
		my @next_loop;
			while (my $next_line = shift @merged_lines) {
				if (abs($line->[7] - $next_line->[7]) < 5) {
				my $diff = $next_line->[0] / $line->[0];
					if ($next_line->[0] / $line->[0] <= 0.05) {
						$next_line->[1] .= "#fuzz";
						push @returned_lines, join("\t", @{$next_line});
					}elsif ($next_line->[0] / $line->[0] <= 0.1) {
						$next_line->[1] .= "#fuzz?";
						push @returned_lines, join("\t", @{$next_line});
					}else {
						push @returned_lines, join("\t", @{$next_line});
					}
				} else {
					push @next_loop, $next_line;
				}
			}
		@merged_lines = @next_loop;
	}	
	
 	
 	return @returned_lines;
}


####################################################################################################
# Iterate through each blat results file. Determine which reads can be mapped to a single location and
# which may map to repetitive regions.
sub parse_alignments {
	
	my ($input_file) = @_;
	open my $data_in, '<', $input_file or die "Unable to open $input_file in parse_alignments: $!\n";
	
	# Initialize the variables that will be used to track the scores from
	# the first line of the input file.
	chomp(my $alignment = <$data_in>);

	my ($sequence_ID) = (split /\t/, $alignment)[9];
	
	# @seq_ID_entries will temporarily store the lines associated with each sequence ID
	# while checking for matche to repetitive elements. It's reset when the script moves on
	# to another sequence ID.
	my @seq_ID_entries = ($alignment); # First, make the first line of the file the first entry in @seq_ID_entries.
	my %out_best; # An array to store the results for _outBest
	my %out_repeat;# An array to store the results for _outRepeat
	# Iterate through the high quality blat results, searching for reads that mapped to multiple
	# places in the genome
	while (my $next_alignment = <$data_in>) {

		
		chomp($next_alignment);
		# Obtain the sequence ID of the next entry. Make $newsequence_ID a junk ID if we're at the end of the file.
		my ($newsequence_ID) = $next_alignment =~ m/\w/ ? (split /\t/, $next_alignment)[9] : "<<EOF>>";
		# If the next sequence ID is the same as the first one, append it to @seq_ID_entries and move on to the next comparison
 		if ($sequence_ID eq $newsequence_ID) {
 			push @seq_ID_entries, $next_alignment;
 		}else{  # If the next sequence ID is different from the first, determine wich, if any, entries
 				# are considered multi-mapped.
 			my $tag_count = @seq_ID_entries;
 			if($tag_count == 1) { # If there is only one instance of a particular sequence ID, print it to _outBest
				# The label 'only_one' will be used to identify sequence IDs that only APPEARED once in the alignment results, as opposed
				# to sequences that aligned more than once but there scores still count them as single-mappers.
 				$alignment .= '#only_one#';					 
 				$out_best{$sequence_ID} = $alignment;
 			}else { # If there are more than 1 entry, check for multi-mapping

 				my $first_entry = shift @seq_ID_entries; # The top scoring entry is the first one to be compared
 				my @multi_mapped = ($first_entry); # Create an array for multi-mapped reads.
 				my ($first_entry_matches) = (split /\t/, $first_entry)[0];
 				
 				while (my $next_entry = shift @seq_ID_entries) { # Check each potentially multi-mapped read in order of descending match score
 					
 					my ($next_entry_matches) = (split /\t/, $next_entry)[0];
 					# Reads are considered multi-mapped if the match score difference of at least one other entry
 					# and the top entry core are within 7 matches of one another. 
 					# If a sequence is multi-mapped, add it to @multi_mapped then check if there are more remaining to be checked
 					if ($first_entry_matches - $next_entry_matches < 7) {
 						push @multi_mapped, $next_entry;
 						# If that was the last element to be checked, add the entries to @out_repeat. If not, move on to the next element
 						if (@seq_ID_entries == 0) {
 							$out_repeat{$sequence_ID} = \@multi_mapped;
 						}
					# If the next entry isn't within 7 matches of the first, the multi-map check is complete
 					}else {
 						# Count the number of mapped reads
 						my $multi_count = @multi_mapped;
 						# If only one read maps, add it to @out_best and end the loop
 						if ($multi_count == 1) {
 							$out_best{$sequence_ID} = $alignment;
 							@seq_ID_entries = undef;
 						# If more than one read map, add them to @out_repeat and end the loop
 						}else {
 							$out_repeat{$sequence_ID} = \@multi_mapped;
 						}
 					}
 				}
			}
			# Prepare the variables for the next comparison
			$alignment = $next_alignment;
			@seq_ID_entries = ($alignment);
			$sequence_ID = $newsequence_ID;
		}
	}

	return(\%out_best, \%out_repeat);
}

####################################################################################################
# Print the results
sub print_best {

	my ($data, $out_file_name) = @_;
	open my $out, '>', $out_file_name;
	foreach (keys %$data) {
		s/#only_one#//;
		say $out $$data{$_};
	}
	close $out;
}
	
sub print_repeat {
	
	my ($data, $out_file_name) = @_;
	open my $out, '>', $out_file_name;
	foreach my $sequence_ID (keys %$data) {
 		my $repeat_array = $$data{$sequence_ID};
 		my $number_of_alignments = @$repeat_array;
  		foreach (@$repeat_array) {
  			say $out "$number_of_alignments\t$_";
  		}
	}
	return;
}

####################################################################################################
# In some instances, one read will map only once while the other maps multiple times. Change the multi-mapper
# to a single-mapper if the other read appears in the alignment data only once.
sub crosscheck_for_singles {

	my ($best_read1, $repeat_read1, $best_read2, $repeat_read2) = @_;
	
	foreach my $sequence_ID (keys %$repeat_read1) {
	
		if (exists $$best_read2{$sequence_ID} and $$best_read2{$sequence_ID} =~ m/#only_one#$/) {
			$$best_read1{$sequence_ID} = $$repeat_read1{$sequence_ID}[0];
			delete $$repeat_read1{$sequence_ID};
		}
	}
	foreach my $sequence_ID (keys %$repeat_read2) {
	
		if (exists $$best_read1{$sequence_ID} and $$best_read1{$sequence_ID} =~ m/\t#only_one#$/) {
			$$best_read2{$sequence_ID} = $$repeat_read2{$sequence_ID}[0];
			delete $$repeat_read2{$sequence_ID};
		}
	}	

	return ($best_read1, $repeat_read1, $best_read2, $repeat_read2);
}

####################################################################################################
# In some instances, alignments to the same integration site will be considered a single mapper
# if the read is long enough but a multi-mapper if it's short. Change these multis to singles if
# the following criteria are met:
#	1. The integration sites are identical.
#	2. The length of the shortest single-mapper is at least 2x + 20 bases longer than the shortest multi-mapper.
#	3. The match score of the single-mapper must be the highest or tied for highest match score.


sub compare_read_lengths {

	my ($best_reads_by_seqID, $repeat_reads) = @_;
	my $best_reads_by_site = convert_to_by_site($best_reads_by_seqID);
	
	#$DB::single=2;	
	COMPARE: foreach my $sequence_ID (keys %$repeat_reads) {
	
		foreach my $alignment (@$repeat_reads{$sequence_ID}) {
			
			foreach my $line (@$alignment) {
				
				my @line_elements = split /\t/, $line;
				my ($matches, $strand, $chr, $start_base) = @line_elements[0,8,13,15];
 				my $site = $chr.$strand.$start_base;
				
				if (exists $$best_reads_by_site{$site}) {
					my @sorted_match_scores = sort keys @$best_reads_by_site{$site};
					if ($matches >= 30 and $sorted_match_scores[0] >= $matches + 20) {
					
						$$best_reads_by_seqID{$sequence_ID} = $line;
						delete $$repeat_reads{$sequence_ID};
						next COMPARE;
					}
				}
			}	
		}
	}
	
	return ($best_reads_by_seqID, $repeat_reads);
}

####################################################################################################
# This subroutine returns hashes of the single-mappers organized by integration site.

sub convert_to_by_site {
	
 	my $data_by_seqID = shift @_;
 	my %data_by_site;
 	
 	foreach my $sequence_ID (keys %$data_by_seqID) {

	 	my $alignment = $$data_by_seqID{$sequence_ID};
	 	my ($matches, $strand, $chr, $start_base) = (split /\t/, $alignment)[0,8,13,15];
		my $site = $chr.$strand.$start_base;
		$data_by_site{$site} = {} unless exists $data_by_site{$site};
		$data_by_site{$site}{$matches} = '' unless exists $data_by_site{$site}{$matches}
		
	}	
	
  	return \%data_by_site;
}

1;