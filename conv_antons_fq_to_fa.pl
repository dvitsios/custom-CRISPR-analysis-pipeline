#!/usr/bin/env perl

use strict;
use warnings;
use threads;



#my @fastq_files = glob("unzipped_original_files/*.fastq");
my @fastq_files = ($ARGV[0]);
my $path_to_save = $ARGV[1];


print "@fastq_files\n";


my @conv_threads = ();

foreach(@fastq_files){

 push @conv_threads, threads->create( sub{
	open(FQ, $_);

	my $fa_name = $_;

	$fa_name =~ s/\.fq//g;
	$fa_name =~ s/.*\/()/$1/g;
	$fa_name = "$path_to_save/$fa_name.fasta";

	open(FA, '>', $fa_name);

	my $line_cnt = 0;
	while(my $line = <FQ>){

	    $line_cnt++;


	    if($line_cnt == 1){
		$line =~ s/^@/>/;
		$line =~ s/\s/_/;;
		print FA $line;
	    } elsif($line_cnt == 2){
		print FA $line;
	    } elsif($line_cnt == 3){
		next;
	    } elsif($line_cnt == 4){
		$line_cnt = 0;
		next;
	    }
	}

	close(FA);
	close(FQ);
 });
}

foreach(@conv_threads){
	$_->join();
}
