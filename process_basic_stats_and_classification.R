input_dir = "needle_output"
#input_dir = "needle_output-antons_bins"

basic_stats_file = paste(input_dir, "All_libs_basic.stats", sep="/")
classifications_file = paste(input_dir, "All_libs_classification.stats", sep="/")


basic_stats_df = read.table(basic_stats_file, header=TRUE) 
classifications_df = read.table(classifications_file, header=TRUE) 


print(head(basic_stats_df))



basic_stats_plots_file = paste(input_dir, "basic_stats_plots.pdf", sep="/")
pdf(basic_stats_plots_file)



valid_discarded_colors = c("#377eb8", "#fbb4ae")


libs = sort(unique(basic_stats_df$LIBRARY))
for(lib in libs){
	
	print(lib)
	cur_lib_basic_df = basic_stats_df[ basic_stats_df$LIBRARY == lib, ]
	print(head(cur_lib_basic_df))

	cur_lib_basic_ratios = data.frame(cur_lib_basic_df$VALID_READS/cur_lib_basic_df$TOTAL_COUNTS, cur_lib_basic_df$DISCARDED_READS/cur_lib_basic_df$TOTAL_COUNTS)

	
	rownames(cur_lib_basic_ratios) = cur_lib_basic_df$MRE
 	colnames(cur_lib_basic_ratios) = c("VALID_READS_RATIO", "DISCARDED_READS_RATIO")

	print(head(cur_lib_basic_ratios))


	# plot valid/discarded % ratios for each MRE
	barplot(as.matrix(t(cur_lib_basic_ratios)), col=valid_discarded_colors, main=paste("Library ", lib, "\nValid (blue) & Discarded (light red)\nreads ratios", sep=""), las = 2, cex.names=0.4)

	# plot reads depth for each MRE
	cur_lib_total_read_depths = data.frame(cur_lib_basic_df$TOTAL_COUNTS)
	rownames(cur_lib_total_read_depths) = cur_lib_basic_df$MRE
	colnames(cur_lib_total_read_depths) = c("TOTAL_COUNTS")

	barplot(as.matrix(t(cur_lib_total_read_depths)), col="#99d8c9", main=paste("Library ", lib, "\nTotal reads depth per MRE", sep=""), las = 2, cex.names=0.4)
#plot.new()

}

dev.off()


classification_plots_file = paste(input_dir, "classification_plots.pdf", sep="/")
pdf(classification_plots_file)


wt_crispr_colors = c("#1f78b4", "#a6cee3", "#33a02c", "#b2df8a")
libs = sort(unique(basic_stats_df$LIBRARY))
for(lib in libs){

        print(lib)
	cur_lib_classifications_df = classifications_df[ classifications_df$LIBRARY == lib, ]


	cur_lib_classification_ratios = cur_lib_classifications_df[ , 3:6]
	print(head(cur_lib_classification_ratios))

	rownames(cur_lib_classification_ratios) = cur_lib_classifications_df$MRE
	colnames(cur_lib_classification_ratios) = c("WT", "WT_SEED_INTACT", "DELETIONS_IN_SEED_RATIO_CRISPR", "TARGET_SEQ_NOT_FOUND_CRISPR")

	barplot(as.matrix(t(cur_lib_classification_ratios)), col=wt_crispr_colors, main=paste("Library ", lib, "\n\n{WT: dark blue}\n{Intact seed region: light blue}\n{CRISPR seed deletion: dark green}\n{CRISPR seed insertion: light green}", sep=""), las=2, cex.names=0.4, cex.main=0.5)
}

dev.off()


final_mir_scores_df = data.frame(MRE=character(), CRISPR_MIR_ratio=numeric(), WT_MIR_ratio=numeric(), ACTIVE_ANNOT=character(), FOLD_CHANGE=numeric(), stringsAsFactors=FALSE) 


class_counts_df = classifications_df[ , colnames(classifications_df) %in% c("LIBRARY", "MRE", "WT_CNT", "WT_SEED_INTACT_CNT", "DELETIONS_IN_SEED_CNT_CRISPR", "TARGET_SEQ_NOT_FOUND_CNT_CRISPR")]
print(head(class_counts_df))

mres = sort(unique(basic_stats_df$MRE))
print(length(mres))

for(mre in mres){

	print(mre)
	tmp_counts_df = class_counts_df[ class_counts_df$MRE == mre, ]
	print(tmp_counts_df)

	gDNA_WT = tmp_counts_df[ tmp_counts_df$LIBRARY == "A", "WT_CNT"] + tmp_counts_df[ tmp_counts_df$LIBRARY == "A", "WT_SEED_INTACT_CNT"] +
                 tmp_counts_df[ tmp_counts_df$LIBRARY == "B", "WT_CNT"] + tmp_counts_df[ tmp_counts_df$LIBRARY == "B", "WT_SEED_INTACT_CNT"]
        print(gDNA_WT)


	cDNA_WT = tmp_counts_df[ tmp_counts_df$LIBRARY == "C", "WT_CNT"] + tmp_counts_df[ tmp_counts_df$LIBRARY == "C", "WT_SEED_INTACT_CNT"] +
		 tmp_counts_df[ tmp_counts_df$LIBRARY == "D", "WT_CNT"] + tmp_counts_df[ tmp_counts_df$LIBRARY == "D", "WT_SEED_INTACT_CNT"] 
	print(cDNA_WT)



	gDNA_MT = tmp_counts_df[ tmp_counts_df$LIBRARY == "A", "DELETIONS_IN_SEED_CNT_CRISPR"] + tmp_counts_df[ tmp_counts_df$LIBRARY == "A", "TARGET_SEQ_NOT_FOUND_CNT_CRISPR"] +
                 tmp_counts_df[ tmp_counts_df$LIBRARY == "B", "DELETIONS_IN_SEED_CNT_CRISPR"] + tmp_counts_df[ tmp_counts_df$LIBRARY == "B", "TARGET_SEQ_NOT_FOUND_CNT_CRISPR"]
	print(gDNA_MT)


	cDNA_MT = tmp_counts_df[ tmp_counts_df$LIBRARY == "C", "DELETIONS_IN_SEED_CNT_CRISPR"] + tmp_counts_df[ tmp_counts_df$LIBRARY == "C", "TARGET_SEQ_NOT_FOUND_CNT_CRISPR"] +
                 tmp_counts_df[ tmp_counts_df$LIBRARY == "D", "DELETIONS_IN_SEED_CNT_CRISPR"] + tmp_counts_df[ tmp_counts_df$LIBRARY == "D", "TARGET_SEQ_NOT_FOUND_CNT_CRISPR"]
	print(cDNA_MT)


if(FALSE){	
	if(gDNA_MT == 0){
		gDNA_MT = 0.00000000001
	}
	if(gDNA_WT == 0){
		gDNA_WT = 0.00000000001 
	} 
	if(cDNA_MT ==0){
		cDNA_MT = 0.00000000001 
	}
	if(cDNA_WT ==0){
		cDNA_WT = 0.00000000001 
	}
}

	CRISPR_MIR_ratio = cDNA_MT/gDNA_MT


	WT_MIR_ratio = cDNA_WT/gDNA_WT

	active_annotation = ""
	fold_change_val = -1

	if(length(gDNA_WT) == 1 && length(cDNA_WT) == 1){
		if(CRISPR_MIR_ratio == "NaN" || CRISPR_MIR_ratio == "Inf" || WT_MIR_ratio == "NaN" || WT_MIR_ratio== "Inf"){
			active_annotation = "NON_DEFINED"
			
		} else if(CRISPR_MIR_ratio > WT_MIR_ratio){
			active_annotation = "ACTIVE"
			fold_change_val = CRISPR_MIR_ratio/WT_MIR_ratio
		} else{
			active_annotation = "Inactive"
			fold_change_val = CRISPR_MIR_ratio/WT_MIR_ratio
		}

			print(paste("CRISPR_MIR_ratio: ", CRISPR_MIR_ratio, sep=""))
			print(paste("WT_MIR_ratio: ", WT_MIR_ratio, sep=""))


		final_mir_scores_df = rbind(final_mir_scores_df, data.frame(mre, CRISPR_MIR_ratio, WT_MIR_ratio, active_annotation, fold_change_val))
	}

}

final_mir_scores_df = final_mir_scores_df[ with(final_mir_scores_df, order(-fold_change_val)), ]

options(digits=3, "scipen"=4)
print(final_mir_scores_df)

final_mir_scores_file = paste(input_dir, "final_mir_scores.csv", sep="/")
write.csv(final_mir_scores_df, final_mir_scores_file, row.names=FALSE)
