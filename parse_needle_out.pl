#!/usr/bin/env perl

use warnings;
use Data::Dumper;


# Parser of the NW output for each MRE.

# There are 2 distinct CRISPR events being detected:
# - DELETIONS_IN_SEED (either actual deletion or modification or one or more nts)
# - INSERTIONS IN SEED (target seq not found in sense alignment)
# I keep a @deletion_blocks_array array with the number of deletions/insertions (recognised as a single event)
# in each of the 6 pre-defined blocks:
# 	* DU: distal upstream
# 	* PU: proximal upstream
# 	* SD: seed region
# 	* MRE: MRE region downstream of the seed region
# 	* PD: proximal downstream
# 	* DD: distal downstream
#  [ DU ]   [PU]  [SD]   [MRE]  [PD]       [DD]
# (1-50nt) (15nt) (6nt) (14nt) (15nt) (end-50 to end)
# Total length of all blocks for deletions detection: 100nt

#
my $needle_out1 = $ARGV[0];
my $needle_out2 = $ARGV[1];
my $target_seq = $ARGV[2];
my $path_to_save = $ARGV[3];
my $mre = $ARGV[4];


open(FH1, $needle_out1);
open(FH2, $needle_out2);



my $ALLOWED_MISMATCHES_THRESHOLD = 6;
my $ALLOWED_INSERTIONS_THRESHOLD = 5;	
my $TRIM_SEQS_TO_100_NT = "T";






my $total_counts = 0;
my $lines_with_many_mismatches = 0;
my $lines_with_many_insertions = 0;
my $index = 0;
# variables for _sense_ sequences
$sense_id = "";
$amplic_aligned_seq = "";
$sense_alignment = "";
$TOO_MANY_MISMATCHES = "F";
$TOO_MANY_AMPLIC_INSERTIONS = "F";

# variables for _antisense_ sequences
$anti_sense_id = "";
$anti_amplic_aligned_seq = "";
$anti_sense_alignment = "";
$anti_TOO_MANY_MISMATCHES = "F";
$anti_TOO_MANY_AMPLIC_INSERTIONS = "F";

# ************************************
# Hash to store deletion position info
# ************************************
#
# Each element is a reference to an array of the form:
# (DU, PU, SD, MRE, PD, DD) where
# DU: distal upstream
# PU: proximal upstream
# SD: seed region
# MRE: MRE region downstream of the seed region
# PD: proximal downstream
# DD: distal downstream
my @deletion_blocks_array = ();




my $discarded_reads_cnt = 0;
my $target_seq_not_found_cnt = 0;
my $pure_wt_cnt = 0;
my $wt_seed_intact_cnt = 0;
my $deletions_in_seed_cnt = 0;

# debug
my $inconsistent_pairs_cnt = 0;

# Output Files (Some classes of Sequences & Stats)
open(WT_INTACT_SEED, ">$path_to_save/wt_seed_intact.seqs");
open(DELETIONS_IN_SEED, ">$path_to_save/deletions_in_seed.seqs");
open(BASIC_STATS_FH, ">$path_to_save/basic.stats");
open(SEQ_CLASSIFICATION_FH, ">$path_to_save/seq_classification.stats");


my $target_start = -1;
my $anti_target_start = -1;
my $substr_start = -1;
my $anti_substr_start = -1;	

while(defined($l1 = <FH1>) && defined($l2 =<FH2>)){


	chomp($l1);
	chomp($l2);

	my @insertions = ();
	my @mismatches = ();	

	my @anti_insertions = ();
	my @anti_mismatches = ();



	$index++;
	if($index == 1){
		$sense_id = $l1;
		$anti_sense_id = $l2;
		
	}
	if($index == 2){
	
		$amplic_aligned_seq = $l1;
		$anti_amplic_aligned_seq = $l2;	

if($TRIM_SEQS_TO_100_NT eq "T"){
$target_start = index($amplic_aligned_seq , $target_seq);
$anti_target_start = index($anti_amplic_aligned_seq , $target_seq);

#print "target_start: $target_start\n";
#print "anti_target_start: $anti_target_start\n";


#print "Original amplic_aligned_seq:\n";
#print "$amplic_aligned_seq\n";
#print "Original anti_amplic_aligned_seq:\n";
#print "$anti_amplic_aligned_seq\n";
$substr_start = $target_start - 50;


if($target_start != -1){
	$amplic_aligned_seq = substr($amplic_aligned_seq, $substr_start, 100); 
}

$anti_substr_start = $anti_target_start - 50;
if($anti_target_start != -1){
	$anti_amplic_aligned_seq = substr($anti_amplic_aligned_seq, $anti_substr_start, 100); 
}





#print "Trimmed amplic_aligned_seq:\n";
#print "$amplic_aligned_seq\n";
#print "Trimmed anti_amplic_aligned_seq:\n";
#print "$anti_amplic_aligned_seq\n";
}



#		@insertions = $l1 =~ /\-/g;
#		@anti_insertions = $l2 =~ /\-/g;
		@insertions = $amplic_aligned_seq =~ /\-/g;
		@anti_insertions = $anti_amplic_aligned_seq =~ /\-/g;

		my $insertions_count =  @insertions;
		my $anti_insertions_count =  @anti_insertions;

		if($insertions_count >= $ALLOWED_INSERTIONS_THRESHOLD){
			$TOO_MANY_AMPLIC_INSERTIONS = "T";
			$lines_with_many_insertions++;
		}

		if($anti_insertions_count >= $ALLOWED_INSERTIONS_THRESHOLD){
			$anti_TOO_MANY_AMPLIC_INSERTIONS = "T";
			if($TOO_MANY_AMPLIC_INSERTIONS ne "T"){
				$lines_with_many_insertions++;
			}
                }


	} elsif($index == 3){

		$total_counts++;
		$sense_alignment = $l1;
		$anti_sense_alignment = $l2;


if($TRIM_SEQS_TO_100_NT eq "T"){
#print "$sense_alignment\n";
#print "$anti_sense_alignment\n";
#print "substr_start: $substr_start\n";
#print "anti_substr_start: $anti_substr_start\n";
if($target_start != -1){
	$sense_alignment = substr($sense_alignment, $substr_start, 100);
}
if($anti_target_start != -1){
	$anti_sense_alignment = substr($anti_sense_alignment, $anti_substr_start, 100);
}

#print "$sense_alignment\n";
#print "$anti_sense_alignment\n";
}


		@mismatches = $sense_alignment =~ /[ACGTUXN]/g;
		my $mismatch_count =  @mismatches;
		@anti_mismatches = $anti_sense_alignment =~ /[ACGTUXN]/g;
		my $anti_mismatch_count =  @anti_mismatches;


		if($mismatch_count >= $ALLOWED_MISMATCHES_THRESHOLD){
			$TOO_MANY_MISMATCHES = "T";
			$lines_with_many_mismatches++;
		}
		if($anti_mismatch_count >= $ALLOWED_MISMATCHES_THRESHOLD){
			$anti_TOO_MANY_MISMATCHES = "T";
			if($TOO_MANY_MISMATCHES ne "T"){
				$lines_with_many_mismatches++;
			}
		}

		$index = 0;


		# **********************
		# ***** PROCESSING *****
		# **********************
		if($TOO_MANY_MISMATCHES eq "T" || $TOO_MANY_AMPLIC_INSERTIONS eq "T"){

			# process anti-sense read since the sense read has been discarded already based
			# on the allowed insertions and deletions thresholds.
			if(!($anti_TOO_MANY_MISMATCHES eq "T" || $TOO_MANY_AMPLIC_INSERTIONS eq "T")){
				process_valid_alignment_hit($anti_amplic_aligned_seq, $target_seq, \@anti_insertions, $anti_sense_alignment);			
				$inconsistent_pairs_cnt++;
			} else{
				$discarded_reads_cnt++;
			}

			$TOO_MANY_MISMATCHES = "F";
			$TOO_MANY_AMPLIC_INSERTIONS = "F"; 
			$anti_TOO_MANY_MISMATCHES = "F";
			$anti_TOO_MANY_AMPLIC_INSERTIONS = "F"; 
			next;
		} else{ # process the valid alignment hit for the sense read only.

			process_valid_alignment_hit($amplic_aligned_seq, $target_seq, \@insertions, $sense_alignment);

		}	
	}
}
close(FH1);
close(FH2);








print "========= DEBUGGING INFO =========\n";
print "total_counts: $total_counts\n";
print "discarded_reads_cnt: $discarded_reads_cnt\n";
print "==================================\n";







my $valid_reads_cnt = $total_counts - $discarded_reads_cnt;

my $aln_quality = (1-($lines_with_many_mismatches/$total_counts))*100;
my $no_amplicon_insertion_ratios = (1-($lines_with_many_insertions/$total_counts))*100;
my $discarded_reads_ratio = ($discarded_reads_cnt/$total_counts)*100;
my $inconsistent_pairs_ratio = ($inconsistent_pairs_cnt/$total_counts)*100;
# --- Valid reads section ---
my $target_seq_not_found_ratio = ($target_seq_not_found_cnt/$valid_reads_cnt)*100;
my $pure_wt_ratio = ($pure_wt_cnt/$valid_reads_cnt)*100;
my $wt_seed_intact_ratio = ($wt_seed_intact_cnt/$valid_reads_cnt)*100;
my $deletions_in_seed_ratio = ($deletions_in_seed_cnt/$valid_reads_cnt)*100;


#if($inconsistent_pairs_ratio > 1){
print "-----------------------------------\n";
print "-- Total counts: $total_counts\n";
printf("> Alignment quality: %.2f%%\n", $aln_quality);
printf("> No amplicon heavy insertion ratio: %.2f%%\n", $no_amplicon_insertion_ratios);
printf("> Discarded reads ratio: %.2f%%\n", $discarded_reads_ratio);   # the union of the reads with too many mismatches and too many amplicon insertions
printf("> Inconsistent reads ratio: %.2f%%\n", $inconsistent_pairs_ratio); # the union of the reads with too many mismatches and too many amplicon insertions
print "\n--- Valid reads section ---\n";
print "-- Valid counts: $valid_reads_cnt\n";
printf("> Pure WT ratio: %.2f%%\n", $pure_wt_ratio); 
printf("> (extra) WT with intact seed ratio: %.2f%%\n", $wt_seed_intact_ratio); 
printf("> Counts with deletions in seed ratio (CRISPR-ed): %.2f%%\n", $deletions_in_seed_ratio); 
printf("> Target seq not found ratio (Insertion at the MRE seed region) [extra CRISPR-ed]: %.2f%%\n", $target_seq_not_found_ratio); 

print "inconsistent_pairs_cnt: $inconsistent_pairs_cnt\n";
#}

print BASIC_STATS_FH "mre\ttotal_counts\tvalid_reads_cnt\tdiscarded_reads_cnt\taln_quality\tno_amplicon_insertion_ratios\n";
printf BASIC_STATS_FH "$mre\t$total_counts\t$valid_reads_cnt\t$discarded_reads_cnt\t%.6f\t%.6f\n", $aln_quality, $no_amplicon_insertion_ratios;



print SEQ_CLASSIFICATION_FH "MRE\tWT\tWT_SEED_INTACT\tDELETIONS_IN_SEED_RATIO_CRISPR\tTARGET_SEQ_NOT_FOUND_CRISPR\tWT_CNT\tWT_SEED_INTACT_CNT\tDELETIONS_IN_SEED_CNT_CRISPR\tTARGET_SEQ_NOT_FOUND_CNT_CRISPR\n";
printf SEQ_CLASSIFICATION_FH "$mre\t%.6f\t%.6f\t%.6f\t%.6f\t$pure_wt_cnt\t$wt_seed_intact_cnt\t$deletions_in_seed_cnt\t$target_seq_not_found_cnt\n", $pure_wt_ratio, $wt_seed_intact_ratio, $deletions_in_seed_ratio, $target_seq_not_found_ratio;



close(BASIC_STATS_FH);
close(SEQ_CLASSIFICATION_FH);








print "(((((((((((())))))))))))\n";
print "*** Deletion blocks ***\n";
foreach(@deletion_blocks_array){
	@arr = @$_;
	print "@arr\n";
}






sub index_all{

	my $query_string = shift;
	my $char = shift;

	my @insertion_positions = ();
	$offset = 0;
	$result = index($query_string, $char, $offset);
	push @insertion_positions, $result; 

	while ($result != -1) {

		$offset = $result + 1;
		$result = index($query_string, $char, $offset);
		
		if($result != -1){
			push @insertion_positions, $result; 
		}
	}

	return(\@insertion_positions);
}





sub process_valid_alignment_hit {

	my $amplic_aligned_seq = shift;
	my $target_seq =shift;

	my $insertions_ref = shift;	
	my @insertions = @$insertions_ref;
	# The name 'sense_alignment' has been kept only for legacy reasons.
	# It can very well be the 'anti_sense_alignment' too if that's
	# the argument provided to the function.
	my $sense_alignment = shift;


	my $target_start = index($amplic_aligned_seq , $target_seq);



	# > Deal with the amplicons with insertions first
	if($target_start == -1){ # This indicates the insertions ratio at the seed region of the MRE.
		$target_seq_not_found_cnt++;

#			# - Get insertion positions in the target seq
#			# NO insertions is an array with a singe value of -1.
			$insertion_positions_ref = index_all($amplic_aligned_seq, '-');
			my @insertion_positions = @$insertion_positions_ref;
			
			# check if insertions don't overlap with the target MRE!				
#			print "[Insertions in Seed Region]:\n";			
#			print "@insertion_positions\n";	
#			print "target_start: $target_start\n";
#			print "amplic_aligned_seq:\n$amplic_aligned_seq\n";
#			print "sense_alignment:\n$sense_alignment\n\n";
			

			my $substr_length_to_check_for_insertion = 15;


#foreach(@insertion_positions){	
#$pos = $_;
#$toRemove = '-';
#$tmp_amplic_aligned_seq =~ s/(.{$pos})$toRemove/$1/;

#$zz = $amplic_aligned_seq;
#$zz =~ s/(.{$pos})/$1_test_/;
#substr($zz, $pos, 0) = 'F';  #insert character into string at specific position
#print "new tmp_amplic_aligned_seq: $zz\n";
	
				$tmp_amplic_aligned_seq = $amplic_aligned_seq;
				$tmp_amplic_aligned_seq =~ s/\-//g;
#				print "tmp_amplic_aligned_seq:\n$tmp_amplic_aligned_seq\n";	

				$tmp_target_start = index($tmp_amplic_aligned_seq , $target_seq);
				

				$substr_for_insertions_check = substr($amplic_aligned_seq, $tmp_target_start, $substr_length_to_check_for_insertion);

#				print "substr_for_insertions_check: $substr_for_insertions_check\n";

				$original_substr_for_insertions_check = $substr_for_insertions_check;
				
				$substr_for_insertions_check =~ s/\-//g;
				$substr_seed_segment = substr($substr_for_insertions_check, 0, length($target_seq));

#				print "substr_seed_segment: $substr_seed_segment\n";
				
				if($substr_seed_segment ne $target_seq){
					print "[Error]: substr_seed_segment ne target_seq!\n";
					exit;
				}

		
				# Current approach:
				# Classify all '-' from the substr_for_insertions_check as CRISPR insertions (handled as 'deletion')
				# inside the seed region.
				# Later on, I should detect if any '-' exceed the overall region. The ratio of these cases is not 
				# expected to be significant though.

#(DU, PU, SD, MRE, PD, DD)			
#@deletion_blocks_array

				@SD = $original_substr_for_insertions_check =~ /\-/g;
				$SD = scalar @SD;
#				print "SD: $SD\n";   	

				my @deletions_vector = (0, 0, $SD, 0, 0, 0);				 
				push(@deletion_blocks_array, \@deletions_vector);
#}							
#			print "============================\n\n";
	
	} else{

		my $num_of_amplicon_insertions = @insertions;
		
		my $trimmed_sense_alignment = $sense_alignment;
		$trimmed_sense_alignment =~ s/^(\-)+//;
		$trimmed_sense_alignment =~ s/(\-)+$//;
		my @deletions_in_trimmed_read = $trimmed_sense_alignment =~ /\-/g;
		my $num_of_deletions_in_trimmed_read = @deletions_in_trimmed_read;



		if($num_of_amplicon_insertions == 0 && $num_of_deletions_in_trimmed_read == 0){
			$pure_wt_cnt++;
			#print "[Pure WT]:\n";
			#print "num_of_amplicon_insertions: $num_of_amplicon_insertions\n";
			#print "num_of_deletions_in_trimmed_read: $num_of_deletions_in_trimmed_read\n";
			#print "amplic_aligned_seq:\n$amplic_aligned_seq\n";
			#print "trimmed_sense_alignment:\n$trimmed_sense_alignment\n\n";
			#print "============================\n";

		} else{

			my $seed_aln_segment = substr($sense_alignment, $target_start, length($target_seq));


			if($seed_aln_segment !~ /[\-AGCTU]/){
				$wt_seed_intact_cnt++;

				print WT_INTACT_SEED "$amplic_aligned_seq\n";
				print WT_INTACT_SEED "$sense_alignment\n\n";
#				print "[WT Intact Seed]:\n";
#				print "num_of_amplicon_insertions: $num_of_amplicon_insertions\n";
#				print "num_of_deletions_in_trimmed_read: $num_of_deletions_in_trimmed_read\n";
#				print "amplic_aligned_seq:\n$amplic_aligned_seq\n";
#				print "trimmed_sense_alignment:\n$trimmed_sense_alignment\n\n";
#				print "seed_aln_segment: $seed_aln_segment\n";	
#				print "============================\n\n";

				

			} else{
				$deletions_in_seed_cnt++;

				print DELETIONS_IN_SEED "$amplic_aligned_seq\n";
				print DELETIONS_IN_SEED "$sense_alignment\n\n";
#				print "[Deletions in Seed]:\n";
#				print "num_of_amplicon_insertions: $num_of_amplicon_insertions\n";
#				print "num_of_deletions_in_trimmed_read: $num_of_deletions_in_trimmed_read\n";
#				print "amplic_aligned_seq:\n$amplic_aligned_seq\n";
#				print "trimmed_sense_alignment:\n$trimmed_sense_alignment\n\n";
#				print "seed_aln_segment: $seed_aln_segment\n";	
#				print "============================\n\n";

			}
				
			
			# ***********************************************************	
			# WT_INTACT_SEED & DELETIONS_IN_SEED are handled together for
			# recording the deletions across the are of interest.

			# get the number of deletions in each of the pre-defined deletion blocks
			get_deletions_per_block($amplic_aligned_seq, $trimmed_sense_alignment, $target_start);
	

		}
	}
}

sub get_deletions_per_block{

	$amplic_aligned_seq = shift;
	$trimmed_sense_alignment = shift;
	$target_start = shift;

	print "amplic_aligned_seq:\n$amplic_aligned_seq\n";
	print "trimmed_sense_alignment:\n$trimmed_sense_alignment\n\n";

}
