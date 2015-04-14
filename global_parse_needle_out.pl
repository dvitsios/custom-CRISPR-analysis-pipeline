#!/usr/bin/env perl

use warnings;
use Data::Dumper;


$out_dir = "./needle_output";
#$out_dir = "./needle_output-antons_bins";


$DEBUG = 0;


# *** initialise ***
my %target_seqs_hash = ();
get_target_seqs_hash();
#print Dumper(\%target_seqs_hash);



opendir $dh, $out_dir;
@lib_dirs = grep {-d "$out_dir/$_" && ! /^\.{1,2}$/} readdir($dh);
close $dh;



#debug
if($DEBUG){
	@lib_dirs = ("C");
}


my $term_cnt = 0;
foreach(@lib_dirs){

	$lib = $_;

	print " \n======== NEW LIBRARY ========\n";
	print "lib: $lib,\t";

	opendir $dh, "$out_dir/$lib";
	#@lib_sub_dirs = grep {-d "$out_dir/$lib/$_" && ! /^\.{1,2}$/} readdir($dh);
	@lib_sub_dirs = ("CRISPR_MRE_Sample");
	close $dh;


	foreach(@lib_sub_dirs){

		$lib_sub_dir = $_;

		print "sub_dir: $lib_sub_dir\t";

		opendir $dh, "$out_dir/$lib/$lib_sub_dir";
		@mres = grep {-d "$out_dir/$lib/$lib_sub_dir/$_" && ! /^\.{1,2}$/} readdir($dh);
		close $dh;

#debug
if($DEBUG){
	@mres = ("G5");
}
		foreach(@mres){
			$mre = $_;

			$cur_target_seq = $target_seqs_hash{$mre}; 

			$needle_file1 = "$out_dir/$lib/$lib_sub_dir/$mre/$mre\_1.needle";
			$needle_file2 = "$out_dir/$lib/$lib_sub_dir/$mre/$mre\_2.needle";

			$path_to_save = "$out_dir/$lib/$lib_sub_dir/$mre";	
			
			print "mre: $mre, cur_target_seq: $cur_target_seq\n";
			print "needle_file1: $needle_file1\n";	
			print "needle_file2: $needle_file2\n";	
	
			$parse_output = `./parse_needle_out.pl $needle_file1 $needle_file2 $cur_target_seq $path_to_save $mre`; 	
			print "$parse_output\n";



#if($term_cnt > 2){			
#exit;
#} else{
#$term_cnt++;
#}

		}
	}
}


# debug only
if($DEBUG){
	exit;
}


open(LIB_BASIC_FH, ">$out_dir/All_libs_basic.stats");
print LIB_BASIC_FH "LIBRARY\tMRE\tTOTAL_COUNTS\tVALID_READS\tDISCARDED_READS\tALN_QUALITY_RATIO\tNO_AMPLICON_INSERTIONS_RATIO\n";

open(LIB_CLASSIFICATION_FH, ">$out_dir/All_libs_classification.stats");
print LIB_CLASSIFICATION_FH "LIBRARY\tMRE\tWT\tWT_SEED_INTACT\tDELETIONS_IN_SEED_RATIO_CRISPR\tTARGET_SEQ_NOT_FOUND_CRISPR\tWT_CNT\tWT_SEED_INTACT_CNT\tDELETIONS_IN_SEED_CNT_CRISPR\tTARGET_SEQ_NOT_FOUND_CNT_CRISPR\n";



open(LIB_DELETIONS_COVERAGE_FOR_ALL_MRES, ">$out_dir/All_libs_deletions_coverage.stats");
print LIB_DELETIONS_COVERAGE_FOR_ALL_MRES "LIBRARY\tMRE\tDU\tPU\tSD\tMRE_3p\tPD\tDD\n";


# mine the extracted stats
foreach(@lib_dirs){

        $lib = $_;


	$lib_sub_dir = "CRISPR_MRE_Sample";
	opendir $dh, "$out_dir/$lib/$lib_sub_dir";
	@mres = grep {-d "$out_dir/$lib/$lib_sub_dir/$_" && ! /^\.{1,2}$/} readdir($dh);
	close $dh;


	foreach(@mres){
		$mre = $_;

		print "mre: $mre\n";
		
		$stats_input_dir = "$out_dir/$lib/$lib_sub_dir/$mre";
		$seq_classification_file = "$stats_input_dir/seq_classification.stats";
		open(SEQ_FH, $seq_classification_file);
		$basic_stats_file = "$stats_input_dir/basic.stats";
		open(BASIC_FH, $basic_stats_file);

		$deletions_summary_file = "$stats_input_dir/deletions_summary.stats";
		open(DELETIONS_FH, $deletions_summary_file);

		if(-s $seq_classification_file && -s $basic_stats_file){

			# get info from 'seq_classification.stats'
			@classification_file_lines = <SEQ_FH>;
			$classification_stats = $classification_file_lines[1];

			print LIB_CLASSIFICATION_FH "$lib\t$classification_stats";




	#		@classification_vals = split('\t', $classification_stats);

	#		$wt_ratio = $classification_vals[1] + $classification_vals[2];
	#		$crispr_ratio = $classification_vals[3] + $classification_vals[4]; 
		
	#		print "wt_ratio: $wt_ratio\n";
	#		print "crispr_ratio: $crispr_ratio\n";

		

			# get info from 'basic.stats'
			@basic_stats_lines = <BASIC_FH>;
			$basic_stats = $basic_stats_lines[1];

			print "basic_stats: $basic_stats\n";

			print LIB_BASIC_FH "$lib\t$basic_stats";

	#		@basic_stats_vals = split('\t', $basic_stats);
	#		$total_counts = $basic_stats_vals[1];
	#		$valid_reads = $basic_stats_vals[2];
	#		$discarded_reads = $basic_stats_vals[3];
		}

		if(-s $deletions_summary_file){
			@deletions_file_lines = <DELETIONS_FH>;
			$deletions_summary = $deletions_file_lines[1];

			print LIB_DELETIONS_COVERAGE_FOR_ALL_MRES "$lib\t$deletions_summary";
		}	
		
	}

}








sub get_target_seqs_hash{

	$target_seqs_file = "data/all_target_seqs.csv";
	open(FH, $target_seqs_file);

	while(<FH>){
		chomp;

		my @vals = split(',', $_);

		my $target_id = $vals[0];
		my $target_seq  = $vals[2];

		$target_seqs_hash{$target_id} = $target_seq;
	}
}
