library(mclust)
args = commandArgs(TRUE)

#Prefix of current experiment
currexpt = args[1]
#cellfloor = No. reads to be considered a cell
#Either specify a number or "mclust" to calculate automatically
cellfloor = args[2]

report2 = read.table(paste0(currexpt,".report.txt"),header=T)
bkgd.ind = grep("bkgd",report2$Tag)
if(length(bkgd.ind) > 0){
	nobkgdmat = report2[-bkgd.ind,]
} else {
	nobkgdmat = report2
}
cellcall = Mclust(data.frame(log10(nobkgdmat$Total)),G=2)
if(cellfloor == "auto"){
  cellfloor = min(nobkgdmat[which(cellcall$classification == 2 & cellcall$uncertainty < 0.05),2])
}else{
  cellfloor = as.numeric(cellfloor)
}
#cellfloor = 1000

subsetmat = nobkgdmat[which(nobkgdmat$Total >= cellfloor),]
write.table(cbind(as.character(rownames(subsetmat)),as.character(subsetmat[,1])),paste0(currexpt,".readdepth.cells.indextable.txt"),row.names=F,col.names=F,sep="\t",quote=F)
subsamples = levels(nobkgdmat$Tag)

pdf(paste0(currexpt,".results.hists.pdf"),height=12,width=12)
par(mfrow=c(2,2))
for(i in 1:length(subsamples)){
  if(subsamples[i] == "bkgd"){next}
  currind = grep(subsamples[i], nobkgdmat$Tag)
  currsubind = grep(subsamples[i], subsetmat$Tag)
  currsub = subset(nobkgdmat,Tag == subsamples[i])
  currsubcells = which(currsub$Total >= cellfloor)
  hist(log10(currsub$Total),breaks=60,col="mediumseagreen",main=subsamples[i],
       xlab="Number of Reads (log10)",las=1,xlim=c(0,7))
  abline(v=log10(cellfloor),lty="dashed",  lwd=2)
  legend("topright",c(paste0("Total Reads: ",sum(currsub$Total)),
                    paste0("\n Total Reads (cells only): ",sum(currsub[currsubcells,2])),
                    paste0("\n Total Barcodes: ",length(currsub$Total)),
                    paste0("\n Number of Cells: ",length(subsetmat$Total[currsubcells])),
                    paste0("\n Median Reads/Cell: ",median(currsub[currsubcells,2])),
                    paste0("\n Range of Reads/Cell: ",min(currsub[currsubcells,2])," - ",max(currsub[currsubcells,2]))),bty="n")
}

hist(log10(nobkgdmat$Total),breaks=60,col="mediumseagreen",main="Overall",
	xlab="Number of Reads (log10)",las=1,xlim=c(0,7))
abline(v=log10(cellfloor),lty="dashed",  lwd=2)
legend("topright",c(paste0("Total Reads: ",sum(nobkgdmat$Total)),
	paste0("Total Reads (cells only): ",sum(subsetmat$Total)),
	paste0("\n Total Barcodes: ",length(nobkgdmat$Total)),
	paste0("\n Number of Cells: ",length(subsetmat$Total)),
	paste0("\n Median Reads/Cell: ",median(subsetmat[,2])),
	paste0("\n Range of Reads/Cell: ",min(subsetmat[,2])," - ",max(subsetmat[,2]))),bty="n")
dev.off()
