#!/usr/bin/env perl

use strict;
use warnings;
use threads;
use threads::shared;
use Data::Dumper;
use String::Util qw(trim);

my @threads = ();

my @libs = ($ARGV[0]);

# depr: 
# my @libs = ("A", "B", "C", "D");


my $MERGE_WXYZ_WITH_ABCD_LIBS = 1;
my %lib_mappings_hash = ( 'W','A', 'X','B', 'Y','C', 'Z','D'); 

print Dumper(\%lib_mappings_hash);

my $MAX_NUM_OF_THREADS = 1;

my $bin_dir = "out_bins";
my $base_out_dir = "./needle_output";
#my $bin_dir = "antons_bins";
#my $base_out_dir = "./needle_output_antons_bins";


my @needle_threads = ();
foreach(@libs){
		
	my $lib = $_;

	my $lib_amplicons_hash_ref = get_amplicons_for_cur_lib($lib);
	my %lib_amplicons_hash = %$lib_amplicons_hash_ref;



#	print "\n\nLibrary: $lib\n";
#	print Dumper($lib_amplicons_hash_ref);


	my @potential_samples_dirs = ("CRISPR_MRE_Sample", "HART_MRE_Sample", "Caribou_CRISPR");

	foreach(@potential_samples_dirs){

		my $cur_sample_dir = $_;
		my $input_dir = "../../make_blast_bins_factory/$bin_dir/$lib/$cur_sample_dir";

#		print "input_dir: $input_dir\n";

		if(-d $input_dir){
			
				
				print "$input_dir exists! Starting samples mapping...\n";
				my @all_input_mres_basenames = glob("$input_dir/*_1.fasta");


				foreach(@all_input_mres_basenames){
					s/$input_dir\///;
					s/_1\.fasta//;
				}

				if($MERGE_WXYZ_WITH_ABCD_LIBS){
					$lib = $lib_mappings_hash{$lib};
				}
				
				#output
				my $wd = "$base_out_dir/$lib";
				`mkdir -p $wd`;
				$wd = $wd."/$cur_sample_dir";
				`mkdir -p $wd`;



#				my $running_threads = 0;

				my @mre_threads = ();


# debug only:			@all_input_mres_basenames = ("A8");
#				@all_input_mres_basenames = ("G7");



				foreach(@all_input_mres_basenames){

#					if($running_threads > $MAX_NUM_OF_THREADS){
#						foreach(@needle_threads){
#						       $_->join();
#						}

#						$running_threads = 0;
#						@needle_threads = ();
#					} else{
#						$running_threads++;
		
	
#						push @needle_threads, threads->create( sub{


							#my %cur_mre_sequences_hash = ();

							my $cur_mre_basename = $_;
							$cur_mre_basename =~ s/ /_/g;
							$cur_mre_basename =~ s/\(//g;
							$cur_mre_basename =~ s/\)//g;
							
							my $cur_amplicon = $lib_amplicons_hash{$cur_mre_basename};
#								print "*** cur_mre: $cur_mre_basename ***\n"; 
#								print "*** cur_amplicon: $cur_amplicon ***\n";

							`mkdir -p $wd/$cur_mre_basename`;


							my $tmp_amplicon_fasta = "$wd/$cur_mre_basename/$cur_mre_basename\_amplicon.fasta";
							open(AMPL_FH, ">$tmp_amplicon_fasta");
							print AMPL_FH ">$cur_mre_basename\n$cur_amplicon";
							close(AMPL_FH);
	

							my $cur_needle_out_file1 = "$wd/$cur_mre_basename/$cur_mre_basename\_1.needle";
							if(!$MERGE_WXYZ_WITH_ABCD_LIBS){	
								open(NEEDLE_FH1, ">$cur_needle_out_file1");
							} else{
								open(NEEDLE_FH1, ">>$cur_needle_out_file1");
							}


							my $cur_needle_out_file2 = "$wd/$cur_mre_basename/$cur_mre_basename\_2.needle";
							if(!$MERGE_WXYZ_WITH_ABCD_LIBS){     
								open(NEEDLE_FH2, ">$cur_needle_out_file2");
							} else{
								open(NEEDLE_FH2, ">>$cur_needle_out_file2");
							}

						
							my $pair1 = "$input_dir/$cur_mre_basename"."_1.fasta";
							open(FILE1, $pair1);

							my $pair1_fh = "$wd/$cur_mre_basename/$cur_mre_basename"."_1_fh.fasta"; #fh: fixed header	

							if(!$MERGE_WXYZ_WITH_ABCD_LIBS){     
								open(FILE1_FH, ">$pair1_fh");
							} else {
								open(FILE1_FH, ">>$pair1_fh");
							}
	
							my $header1 = "";
							while(my $l1 = <FILE1>){
                                                                chomp($l1);
                                                                if($l1 =~ /^>/){
                                                                        $header1 = $l1;
								        $header1 =~ s/\:/_/g;	
                                                                }else{
                                                                        my $seq1 = $l1;
                                                                        print FILE1_FH "$header1\n$seq1\n";
                                                                }                                                               
                                                        }
                                                        close(FILE1);
                                                        close(FILE1_FH);



							open(PROC,"/homes/aje/src/emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0/emboss/needle -asequence $tmp_amplicon_fasta -bsequence $pair1_fh -gapopen 10 -gapextend=0.5 -outfile stdout -aformat markx2 -brief=N -awidth=2000000|");
							#print "*** Pair-1 ***\n";
							my $PRINT_LINE = "FALSE";
							while(<PROC>){
                 						chomp;
								if($_ =~ /^# 2\:/){
									my $query_seq_id = $_;
									$query_seq_id =~ s/^# 2\:\s//;
									print NEEDLE_FH1 ">$query_seq_id\n";

								} elsif($_ =~ /^(\s)+$cur_mre_basename/){
									# found subject (amplicon) sequence line
									$PRINT_LINE = "TRUE";
								
									my $cur_line = trim($_);
	
									my @amplic_seq_vals = split /\s+/ , $cur_line;
									my $amplic_seq = $amplic_seq_vals[1];

									#$amplic_seq = substr($amplic_seq, rindex($amplic_seq, "\s"));
									print NEEDLE_FH1 "$amplic_seq\n";

								} elsif($PRINT_LINE eq "TRUE"){
									my $subject_seq = trim($_);

									my @subj_seq_vals = split /\s+/ , $subject_seq;
									$subject_seq = $subj_seq_vals[1];

									print NEEDLE_FH1 "$subject_seq\n";

									$PRINT_LINE = "FALSE";
								}
                					}
							close(PROC);	
	

							my $pair2 = "$input_dir/$cur_mre_basename"."_2.fasta";
							open(FILE2, $pair2);

							my $pair2_rc = "$wd/$cur_mre_basename/$cur_mre_basename"."_2_rc.fasta";
							if(!$MERGE_WXYZ_WITH_ABCD_LIBS){     
								open(FILE2_RC, ">$pair2_rc"); # reverse complement fasta of antisense sequences!
							} else{
								open(FILE2_RC, ">>$pair2_rc"); # reverse complement fasta of antisense sequences!
							}


							my $header2 = "";
							while(my $l2 = <FILE2>){
								chomp($l2);
								if($l2 =~ /^>/){
									$header2 = $l2;
								        $header2 =~ s/\:/_/g;	
								}else{
									my $seq2 = revcompl($l2);

									print FILE2_RC "$header2\n$seq2\n";
								}								
							}
							close(FILE2);
							close(FILE2_RC);

							open(PROC,"/homes/aje/src/emboss.open-bio.org/pub/EMBOSS/EMBOSS-6.6.0/emboss/needle -asequence $tmp_amplicon_fasta -bsequence $pair2_rc -gapopen 10 -gapextend=0.5 -outfile stdout -aformat markx2 -brief=N -awidth=2000000|");
							#print "*** Pair-2 ***\n";
							while(<PROC>){
                 						chomp;
								if($_ =~ /^# 2\:/){
									my $query_seq_id = $_;
									$query_seq_id =~ s/^# 2\:\s//;
									print NEEDLE_FH2 ">$query_seq_id\n";

								} elsif($_ =~ /^(\s)+$cur_mre_basename/){
									# found subject (amplicon) sequence line
									$PRINT_LINE = "TRUE";
								
									my $cur_line = trim($_);
	
									my @amplic_seq_vals = split /\s+/ , $cur_line;
									my $amplic_seq = $amplic_seq_vals[1];

									#$amplic_seq = substr($amplic_seq, rindex($amplic_seq, "\s"));
									print NEEDLE_FH2 "$amplic_seq\n";

								} elsif($PRINT_LINE eq "TRUE"){
									my $subject_seq = trim($_);

									my @subj_seq_vals = split /\s+/ , $subject_seq;
									$subject_seq = $subj_seq_vals[1];

									print NEEDLE_FH2 "$subject_seq\n";
									$PRINT_LINE = "FALSE";
								}
								
                					}
							close(PROC);	


							print "Finished alignment for $cur_mre_basename!\n\n";
#						});
#					}
				}
		}

	}
}


#foreach(@needle_threads){
#	$_->join();
#}



sub get_amplicons_for_cur_lib{
	my $lib = shift;

	my %amplicons_hash = ();

	open(FH, "../../make_blast_bins_factory/db/$lib\_all_amplicons.fa");
		
	my $header = "";
	while(<FH>){

		chomp;
		if($_ =~ /^>/){
			$header = $_;
			$header =~ s/>//g;
		} else{
			$amplicons_hash{$header} = $_;
		}
	}


	return(\%amplicons_hash);
}


sub revcompl{

        my $seq=$_[0];
        $seq=~ s/A/X/g;
        $seq=~ s/T/Y/g;
        $seq=~ s/U/Y/g;
        $seq=~ s/G/Z/g;
        $seq=~ s/C/W/g;
        $seq=~ s/X/T/g;
        $seq=~ s/Y/A/g;
        $seq=~ s/Z/C/g;
        $seq=~ s/W/G/g;

        return(reverse($seq));
}


