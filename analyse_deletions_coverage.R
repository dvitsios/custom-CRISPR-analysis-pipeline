library(gplots)
library(RColorBrewer)

df = read.table(file="needle_output/All_libs_deletions_coverage.stats", stringsAsFactors=FALSE, sep="\t", header=T)
print(head(df))

libs = unique(df$LIBRARY)
for(lib in libs){
print(paste("lib:", lib))
    sub_df = subset(df, "LIBRARY"="D")
    print(head(sub_df))
    stop()

    mres = rownames(df)
    del_colnames = colnames(df)
    del_colnames = del_colnames[3:length(del_colnames)]
    print(del_colnames)


    df[ , 1] = NULL
    df[ , 1] = NULL
    print(head(df))
    str(df)

#df = as.data.frame(lapply(df,as.numeric))

    mat = as.matrix(df)
#apply(mat, 1, as.numeric)
#print(class(mat))


    print(head(df))


    heatmap.2(mat,  col=brewer.pal(8,"GnBu"), Colv = del_colnames, Rowv = mres, trace="none", dendrogram="none") 
}
