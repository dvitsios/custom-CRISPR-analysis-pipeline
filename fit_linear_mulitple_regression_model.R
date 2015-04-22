library(gplots)
library(RColorBrewer)
library(plotrix)

# Read input
df = read.table(file="needle_output/All_libs_deletions_coverage.stats", stringsAsFactors=FALSE, sep="\t", header=T)
mir_scores_df = read.csv(file="needle_output/final_mir_scores.csv")





mir_scores_df = mir_scores_df[ ,colnames(mir_scores_df) %in% c("mre", "fold_change_val")]
colnames(mir_scores_df)[1] = "MRE"



# => I should examine the effect of the fold change of the mean of deletions ratios
# between the gDNA and cDNA libraries.



# ****************************************************************
# *** get the df with the average deletion ratios from all libs ***
# ****************************************************************

print(head(df))
print(head(mir_scores_df))

ab_df = df[ df$LIBRARY %in% c("A", "B"), ]
#print(ab_df)
cd_df = df[ df$LIBRARY %in% c("C", "D"), ]
#print(cd_df)


print("******************************************")

ab_aggr_avg_df = aggregate(cbind(DU, PU, SD, MRE_3p, PD, DD) ~ MRE, data=ab_df, FUN=mean)
cd_aggr_avg_df = aggregate(cbind(DU, PU, SD, MRE_3p, PD, DD) ~ MRE, data=cd_df, FUN=mean)


print(head(ab_aggr_avg_df))
print(head(cd_aggr_avg_df))



#fold_change_df = matrix(0, ncol=0, nrow=nrow(ab_aggr_avg_df))
#fold_change_df = data.frame(fold_change_df)

fold_change_df = data.frame(MRE=ab_aggr_avg_df$MRE)

for(c in 2:ncol(ab_aggr_avg_df)){
	tmp_col = cd_aggr_avg_df[ , c] / ab_aggr_avg_df[ , c] 
	fold_change_df = cbind(fold_change_df, tmp_col)
}

colnames(fold_change_df) = c("MRE", "DU_fold_change", "PU_fold_change", "SD_fold_change", "MRE_3p_fold_change", "PD_fold_change", "DD_fold_change")
print(head(fold_change_df))



print(head(mir_scores_df))


merged_fold_change_df = merge(fold_change_df, mir_scores_df, by="MRE")
merged_fold_change_df = merged_fold_change_df[ with(merged_fold_change_df, order(-fold_change_val)), ]

print(head(merged_fold_change_df))


fold_change_fit = lm(fold_change_val ~ DU_fold_change + PU_fold_change + SD_fold_change + MRE_3p_fold_change + PD_fold_change + DD_fold_change, data=merged_fold_change_df)
#fold_change_fit = lm(fold_change_val ~ DU_fold_change + PU_fold_change + MRE_3p_fold_change + PD_fold_change + DD_fold_change, data=merged_fold_change_df)


summary(fold_change_fit)
plot(fold_change_fit, scale="r2")

library(leaps)
library(car)

subsets(fold_change_fit, statistic="rsq")


stop()

# ****************************************************************
# *** get the df with the average deletion ratios from all libs ***
# ****************************************************************


aggr_avg_df = aggregate(cbind(DU, PU, SD, MRE_3p, PD, DD) ~ MRE, data=df, FUN=mean)
