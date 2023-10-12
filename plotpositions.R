readlinkerpos = read.table("positions.tsv",header=F,sep="\t")
pdf("linkerpositions.pdf")
hist(readlinkerpos[[1]], xlim=c(1,50), breaks=1)
dev.off()
