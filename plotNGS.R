# Mon 14 Jan 2013 11:02:28 CET: new v takes names from arguments
options(echo=FALSE)

dat = read.table(file = "sample.wintheta", header = T, stringsAsFactor = F)
chrs = unique(dat[,1])

pdf("plot.pdf")

par(mfrow=c(3,2))
for(i in 1:length(chrs)) {
    datToplot = dat[dat[,1] == chrs[i],]
    for (j in 6:11) {
        plot ( datToplot[,c(2,j)], xlab = paste("chr",chrs[i],", window"), ylab = colnames(dat)[j], main=colnames(dat)[j] )
	lines (lowess(datToplot[,c(2,j)]), col=2) 
    }
}
dev.off()



