library(gplots)
library(RColorBrewer)


mir_scores_df = read.csv("needle_output/final_mir_scores.csv")

sorted_mres = mir_scores_df$mre
print(head(sorted_mres))

libs = c("A", "B", "C", "D")
global_sorted_mres = vector()
for(lib in libs){
	for(mre in sorted_mres){
		new_mre = paste(lib, mre, sep="_")
		global_sorted_mres = c(global_sorted_mres, new_mre) 
	}
}
print(global_sorted_mres)

df = read.table(file="needle_output/All_libs_deletions_coverage.stats", stringsAsFactors=FALSE, sep="\t", header=T)
print(head(df))


get_deletions_coverage_across_all_libs <- function(){

		
	mres = paste(df$LIBRARY, df$MRE, sep="_")
	rownames(df) = mres
	del_colnames = colnames(df)
	del_colnames = del_colnames[3:length(del_colnames)]

	df = df[ match(global_sorted_mres, row.names(df)), ]
#	df = df[ order(row.names(df)), ]
	# delete columns: "LIBRARY" and "MRE"
	df[ , 1] = NULL
	df[ , 1] = NULL
	print(head(df))

	mat = as.matrix(df)


	#heatmap.2(mat,  col=brewer.pal(8,"GnBu"), Colv = del_colnames, trace="none", dendrogram="none", sepwidth=c(0.05, 0.45), lhei = c(1, 5), cexRow=0.1,
	heatmap.2(mat,  col=brewer.pal(8,"GnBu"), Colv = del_colnames, Rowv=FALSE, trace="none", dendrogram="none", sepwidth=c(0.35, 0.45), cexRow=0.1, cexCol=0.4,
		lhei = c(0.05,0.15))
		# colsep=c(1:6),rowsep=(1:nrow(mat)), margins=c(1,8) ) 
}



get_deletions_coverage_for_lib <- function(lib_df, lib){

	mres = lib_df$MRE
	rownames(lib_df) = mres
	del_colnames = colnames(df)
        del_colnames = del_colnames[3:length(del_colnames)]

lib_df = lib_df[ match(sorted_mres, lib_df$MRE), ]
#	lib_df = lib_df[ order(sorted_mres), ]

	# delete columns: "LIBRARY" and "MRE"
        lib_df[ , 1] = NULL
        lib_df[ , 1] = NULL
        print(head(lib_df))


        mat = as.matrix(lib_df)
	
	heatmap.2(mat,  col=brewer.pal(8,"GnBu"), Colv = del_colnames, Rowv=FALSE, trace="none", dendrogram="none", sepwidth=c(0.35, 0.45), cexRow=0.3, cexCol=0.4,
                lhei = c(0.05,0.15), main=paste("Deletions coverage, Library:", lib))	

}



pdf(file="needle_output/deletions_coverage.pdf")


# get deletions coverage for all libs together
get_deletions_coverage_across_all_libs()


# get deletions coverage for each lib separately
libs = sort(unique(df$LIBRARY))
for(lib in libs){
 print(paste("lib:", lib))

    lib_df = subset(df, LIBRARY==lib)
    print(head(lib_df))
   
	get_deletions_coverage_for_lib(lib_df, lib) 
}











dev.off()





