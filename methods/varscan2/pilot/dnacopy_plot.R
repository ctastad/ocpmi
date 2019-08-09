library(DNAcopy)

cn <- read.table("pilot_output.copynumber.called",header=T)
cna.object <- CNA(cbind(cn$normal_depth, cn$tumor_depth), cn$chrom, cn$chr_start, data.type="logratio", sampleid=c("pm2","pm15"))

# Smoothing object
cna.object.smooth <- smooth.CNA(cna.object)

#Segmentation at default parameters
cna.object.smooth.segd <- segment(cna.object.smooth, verbose=1)

#Plot whole studies
pdf('rplot.pdf')
plot(cna.object.smooth.segd, plot.type="w")
dev.off()
