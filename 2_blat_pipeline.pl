#!/usr/bin/perl

# 2_blat_pipeline.pl
# Dr. Xiaolin Wu & Daria W. Wells

# This script is used after the ISA debarcoding script.  It requires BLAT to
# be installed and hg19 to be downloaded. It reads the demultiplexed and
# split LTR sequences, passes them through various quality filters, performs BLAT alignment,
# then reports the highest quality sequences that lie on the same chromosome within 10kb
# of each other.  The user may choose to filter reads sequence by defining $filter_sequence.
# The sequence can be a regular expression to filter out multiple sequences.

use 5.18.4;
use strict;
use warnings;

use Carp;
use File::Basename;
use Getopt::Long qw(GetOptions);
Getopt::Long::Configure qw(gnu_getopt);

# Assign the path to the directory containing this script to a scalar
my $pipeline;
BEGIN 
{
	$pipeline = dirname(__FILE__);
	unshift @INC, $pipeline;
}

# Import the analysis modules
use ISAbatchFilter;
use ISAblatFilter;
use ISAdataReformat;

##########################################################################################
# Parse command line options
# Check for command line options.
# -v [FILE]					Align to a viral reference [FILE] fasta format
# -f '[SEQUENCE]'		Filter reads that match [SEQUENCE] before alignment.
# -b [PATH] Set the path to the command line blat program
# -g [PATH] Set the path tho the reference chromosome fasta files.

# Die if no data was specified on the command line.
# die "\nNo data files specified on the command line\n" unless @ARGV;

our ($opt_virus, $opt_filter, $opt_blat, $opt_genome);
	
our @log_lines;
GetOptions(
	'v=s' => \$opt_virus,
	'b=s' => \$opt_blat,
	'f=s' => \$opt_filter,
	'g=s' => \$opt_genome,
);


# The paths to blat and the reference genome must be set before the pipeline can be used.
if ($opt_blat or $opt_genome) {
	set_blat_and_genome($pipeline);
	exit;
}

# Import the paths for blat and the reference genome
open my $paths, '<', "$pipeline/ALIGNMENT" or croak "File $pipeline/ALIGNMENT not found. Create it by first running 2_blat_pipeline.pl with the -b and -g options.\n";
my $path2blat = <$paths>; # The path of the BLAT application
my $path2genome = <$paths>; # The path to the human chromosome files
chomp($path2blat, $path2genome);
close $paths;

# Define the list of chromosomes to be analyzed

my @chromosomes;
opendir(my $genome, $path2genome);
while (my $chr = readdir($genome)) {
	push @chromosomes, "$path2genome/$chr" if $chr =~ m/chr\w+?\.fa/;
}
closedir($genome);


# Define more variables
my $data_file = shift @ARGV; # The name of the file to be analyzed. From the command line
my $path2ooc = $pipeline.'/11.ooc'; # Required for blat
my $path2uclust = $pipeline.'/uclust'; # Path to the uclust program
my ($LTR, $name) = ($data_file =~ /([35]LTR)\/(.*)\.txt\.gz/); # Remove .txt.gz from the file name.
my $folder = $LTR."/".$name; # Where all report files will be stored
mkdir "$folder" or die "\nFolder $folder already exists.\nWere you about to overwrite your data?\n"; # Create the folder where the results files will be stored.


# Write to log
my $nowtime = localtime;
open my $log_fh,">>","$folder/blat_log.txt" or die "\nUnable to print to log file in main script: $!";
print $log_fh "Analysis started at $nowtime\n\n$data_file\n\n";
close $log_fh;

print "Processing $data_file...\n";

# Filter any sequences that match the filter sequence, if specified. 
if (defined $opt_filter) {
	ISAbatchFilter::batch_filter_1_internal($opt_filter, $LTR, $data_file, $folder);
	$data_file = "PEread_LTR_NonInternal_linker.txt";
}
# Filter out reads shorter than 20 bases long.
ISAbatchFilter::batch_filter_2_length20($data_file,$LTR,$folder);

# Filter out reads with quality scores less than 20.
ISAbatchFilter::batch_filter_3_q20($LTR,$folder);

# Make fasta files for alignment.
ISAbatchFilter::make_unique_pair_fasta($LTR,$folder);

# Align the fasta files to human and HIV
# Example alignment command:
# /Applications/Blat/blat /Path/to/hg19db/hg19chromFA/chr1.2bit ./PEread_LTR_NonInternal_linker_long_UniquePairWct_F.fa -ooc=11.ooc PEread_LTR_NonInternal_linker_long_UniquePair_F_chr1.psl  -minScore=16 -stepSize=8
# Align the forward and reverse reads. $s stands for strand. 	 
foreach my $read (qw(read1 read2)) {
	my @args;
	my $path2data = "$folder/$read\_fasta_for_alignment.fa";
	foreach my $path2chr (@chromosomes) {
		
		$path2chr =~ m/(chr\w+?)\.fa$/;
		my $chr = $1;
		my $path2output = "$folder/$read\_blat_results_$chr.psl";
		# Create the command for each chromosome.
		@args = ($path2blat, $path2chr, $path2data, "-ooc=$path2ooc", $path2output, "minScore=16 -stepSize=8");
		system("@args");

	}
	if ($opt_virus) {
		my $path2output = "$folder/$read\_blat_results_HIV.psl";
		@args[1,4] = ($opt_virus, $path2output);
		system("@args");
		
		
	}
	# Concatenate the BLAT output .psl files from each chromosome into one file text file
	# then delete all of the .psl files. 
	system("cat $folder/*.psl >>$folder/$read\_blat_results.txt");
	system("rm $folder/*.psl");

}

# Filter out alignments from the BLAT output that start more than 3 bases into the 
# sequence and are shorter than 20 bases.	
print "\nFiltering alignments shorter than 20 bases...\n";   
ISAblatFilter::blat_filter1($_,$folder) foreach qw(read1 read2);

# Filtering alignments to repetitive regions 
print "Filtering alignments to repetitive regions...\n";
ISAblatFilter::blat_filter2($folder);

# Combine the read 1 and read 2 blat results
print "Combining paired-end data...\n";
my $repeat_alignments = ISAblatFilter::blat_filter3($folder);

# Filter out reads where the left and right coordinates aren't on the same chromosome, 
# are on the same strand, or are more than 10kb apart.
print "Filtering mismatched read pair alignments...\n";
ISAblatFilter::blat_filter4($folder,$repeat_alignments);

# Merge integration sites if they meet all of the following criteria:
	#	1.	The breakpoints are within 10 bases of each other.
	#	2.	The MIDs differ by no more than two bases. (Hamming distance)
	#	3.	The site with the larger count is greater than double the smaller count plus 1. (count score)
print "Merging reads with small IS or breakpoint base variation...\n";
ISAblatFilter::blat_filter5($folder);


# Predict which integration sites might be the result of mis-priming.
print "Adding mis-priming and sequence info to the results...\n";
ISAdataReformat::annotate_blat_file($LTR,$name,$path2genome,$folder);

# Sort the fasta file by sequence length for uClust.
print "Sorting the data for Uclust...\n";
ISAdataReformat::sort_fasta_by_length($folder);
# Run uClust.
system("$path2uclust --input $folder/read1_sorted_by_length_for_Uclust.fa --uc $folder/uclust_results.uc --id 0.90");
# Reformat the uClust output

print "Re-sorting the Uclust data...\n";
ISAdataReformat::resort_Uclust($name,$folder);


exit;

# This subroutine initializes or changes the paths used for the blat step.
sub set_blat_and_genome {

	my $pipeline = shift;
	my ($line1, $line2);
	my $alignment_file = "$pipeline/ALIGNMENT";
	
	if (-e $alignment_file) {
		open my $paths, '<', $alignment_file;
		$line1 = <$paths>;
		$line2 = <$paths>;
		chomp($line1, $line2);
		close $paths;
	}
	
	no warnings 'uninitialized';
	
	$_ =~ s/\/$// foreach ($line1, $line2);
	$line1 = defined $opt_blat ? $opt_blat : $line1;
	$line2 = defined $opt_genome ? $opt_genome : $line2;
	
	open my $out, '>', $alignment_file;
	say $out "$line1\n$line2";
	close $out;
	
	return;
}