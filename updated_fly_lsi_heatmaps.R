#module load gmp/5.0.2 mpfr/3.1.0 mpc/0.8.2 gcc/4.9.1 R/3.2.1
library(Matrix)
library(limma)
library(RColorBrewer)
library(proxy)
library(irlba,lib.loc="~/R/x86_64-unknown-linux-gnu-library/3.1/")
library(gplots)
library(data.table)
library(Rtsne)
library(densityClust)
library(methods)

#Function to perform LSA
my_lsa <- function (x, dims = 2) 
{
  n <- nrow(x)
  p <- ncol(x)
  SVD = irlba(x, dims, dims)
  if (is.function(dims)) {
    dims = dims(SVD$d)
  }
  if (dims < 2) 
    dims = 2
  if (any(SVD$d <= sqrt(.Machine$double.eps))) {
    warning("[lsa] - there are singular values which are zero.")
  }
  space = NULL
  space$tk = SVD$u[, 1:dims]
  space$dk = SVD$v[, 1:dims]
  space$sk = SVD$d[1:dims]
  rownames(space$tk) = rownames(x)
  rownames(space$dk) = colnames(x)
  class(space) = "LSAspace"
  return(space)
}

gz = FALSE
numclust = 5
numcomponents = 6
cellcutoff=0.1
sitecutoff=20000
inmat = "/net/shendure/vol10/projects/scATAC/nobackup/FlyEmbryogenesis/sc_drone_tha/bowtie/SCatac_tha.bowtie.2to4.fixed.nodups.2kbwindows.dhsmatrix.txt"
outhists = paste0("/net/shendure/vol10/projects/scATAC/nobackup/FlyEmbryogenesis/sc_drone_tha/bowtie/SCatac_tha.bowtie.2to4.fixed.nodups.2kbwindows.",cellcutoff,"cells.",sitecutoff,"sites.",numcomponents,"components.",numclust,"clusters.binaryhists.pdf")
outheatmap = paste0("/net/shendure/vol10/projects/scATAC/nobackup/FlyEmbryogenesis/sc_drone_tha/bowtie/SCatac_tha.bowtie.2to4.fixed.nodups.2kbwindows.",cellcutoff,"cells.",sitecutoff,"sites.",numcomponents,"components.",numclust,"clusters.heatmap.jpg")
sitesout = paste0("/net/shendure/vol10/projects/scATAC/nobackup/FlyEmbryogenesis/sc_drone_tha/bowtie/SCatac_tha.bowtie.2to4.fixed.nodups.2kbwindows.",cellcutoff,"cells.",sitecutoff,"sites.",numcomponents,"components.",numclust,"clusters.siteclusters.txt")
cellsout = paste0("/net/shendure/vol10/projects/scATAC/nobackup/FlyEmbryogenesis/sc_drone_tha/bowtie/SCatac_tha.bowtie.2to4.fixed.nodups.2kbwindows.",cellcutoff,"cells.",sitecutoff,"sites.",numcomponents,"components.",numclust,"clusters.cellclusters.txt")
tsneout = paste0("/net/shendure/vol10/projects/scATAC/nobackup/FlyEmbryogenesis/sc_drone_tha/bowtie/SCatac_tha.bowtie.2to4.fixed.nodups.2kbwindows.",cellcutoff,"cells.",sitecutoff,"sites.",numcomponents,"components.",numclust,"clusters.tsne.pdf")
compplotsout = paste0("/net/shendure/vol10/projects/scATAC/nobackup/FlyEmbryogenesis/sc_drone_tha/bowtie/SCatac_tha.bowtie.2to4.fixed.nodups.2kbwindows.",cellcutoff,"cells.",sitecutoff,"sites.",numcomponents,"components.",numclust,"clusters.componentplots.pdf")
promwindows = read.table("/net/shendure/vol10/projects/scATAC/nobackup/genomes/dm3/dm3.2kb.windows.within500bpoftss.whitelist.bed")
oldcells = read.table("/net/shendure/vol10/projects/scATAC/nobackup/FlyEmbryogenesis/sc_drone_tha/sorting/heatmaps/SCatac_tha.2to4.true.nodups.2kbwindows.0.1cells.20000sites.6components.5clusters.cellclusters.txt")
weirdcolors = c("darkmagenta","cyan2","black","gray","darkorange1")

if(gz){
  outprefix = strsplit(inmat,"[.]dhsmatrix[.]txt[.]gz")[[1]]
}else{
  outprefix = strsplit(inmat,"[.]dhsmatrix[.]txt")[[1]]
}
rdsfile = paste0(outprefix,".dhsmatrix.rds")
sparserdsfile = paste0(outprefix,".dhsmatrix.binary.sparse.rds")

if(file.exists(sparserdsfile)){
  bigmat.bin = readRDS(sparserdsfile)
  namelist = readRDS(paste0(outprefix,".dhsmatrix.binary.namelist.rds"))
  label = readRDS(paste0(outprefix,".dhsmatrix.binary.sitelist.rds"))
  annot = rownames(label)
}else{
  #Read in matrix
  if(file.exists(rdsfile)){
    bigmat = readRDS(rdsfile)
  }else{
    if(gz == TRUE){
      bigmat = fread(paste0("zcat ",inmat),header=T)
    } else {
      bigmat = fread(inmat,header=T)
    }
    bigmat = as.data.frame(bigmat)
    saveRDS(bigmat,rdsfile)
  }
  print("Formatting data...")
  sortmat = as.matrix(bigmat[,-c(1:4)])
  #Make a list of sites
  sites = strsplit2(bigmat[,4],"_")
  label = paste0(sites[,1],":",sites[,2],"-",sites[,3])
  annot = label
  currnonzs = unlist(apply(sortmat,2,function(x) which(x > 0)))
  namers = gsub('[[:digit:]]+', '', names(currnonzs))
  namenums = c(1:length(unique(namers)))[unlist(as.factor(namers))]
  mastermat = cbind(currnonzs,namenums)
  namelist = unique(namers)
  bigmat.bin = sparseMatrix(i=mastermat[,1],j=mastermat[,2],x=rep(1,times=dim(mastermat)[1]),dims=dim(sortmat))
  saveRDS(bigmat.bin,sparserdsfile)
  saveRDS(namelist,paste0(outprefix,".dhsmatrix.binary.namelist.rds"))
  saveRDS(label,paste0(outprefix,".dhsmatrix.binary.sitelist.rds"))
  rm(currnonzs)
  rm(mastermat)
  rm(namers)
  rm(namenums)
  gc()
}

colnames(bigmat.bin) = namelist
rownames(bigmat.bin) = annot

num_cells_ncounted = rowSums(bigmat.bin)
num_sites_ncounted = colSums(bigmat.bin)

pdf(outhists)
hist(log10(num_cells_ncounted),main="Number of Cells Each Site is Observed In",breaks=50)
abline(v=log10(num_cells_ncounted[order(num_cells_ncounted,decreasing=T)[sitecutoff]]),lwd=2,col="indianred")
hist(log10(num_sites_ncounted),main="Number of Sites Each Cell Uses",breaks=50)
abline(v=log10(quantile(num_sites_ncounted,probs=cellcutoff)),lwd=2,col="indianred")
annot.ncounts = annot[which(num_cells_ncounted >= num_cells_ncounted[order(num_cells_ncounted,decreasing=T)[sitecutoff]])]
ncounts = bigmat.bin[which(num_cells_ncounted >= num_cells_ncounted[order(num_cells_ncounted,decreasing=T)[sitecutoff]]),]
new_counts = colSums(ncounts)
ncounts = ncounts[,new_counts >= quantile(new_counts,probs=cellcutoff)]
hist(log10(new_counts),main="Number of Sites Each Cell Uses",breaks=50)
abline(v=log10(quantile(new_counts,probs=cellcutoff)),lwd=2,col="indianred")
dev.off()

annot.ncounts = annot.ncounts[rowSums(ncounts) > 0]
ncounts = ncounts[rowSums(ncounts) > 0,]
ncounts = Matrix(ncounts,sparse=T)

nfreqs <- t(t(ncounts) / colSums(ncounts))
tf_idf_counts <- nfreqs * log(1 + ncol(ncounts) / rowSums(ncounts))
set.seed(0)
LLL <- my_lsa(tf_idf_counts,numcomponents)
sk_diag <- matrix(0, nrow=length(LLL$sk), ncol=length(LLL$sk))
diag(sk_diag) <- LLL$sk
sk_diag[1,1] = 0
LLL_d <- t(sk_diag %*% t(LLL$dk))
hclust_cells <- hclust(proxy::dist(LLL_d, method="cosine"), method="ward.D2")
TTT_d <- t(sk_diag %*% t(LLL$tk))
hclust_genes <- hclust(proxy::dist(TTT_d, method="cosine"), method="ward.D2")

genes_tree_cut <- cutree(hclust_genes, numclust)
gene_pal <- weirdcolors
cells_tree_cut <- cutree(hclust_cells, numclust)
cell_pal <- weirdcolors

LSI_out <-  t(t(sk_diag %*% t(LLL$dk)) %*% t(LLL$tk))
LSI_out <- t(scale(t(LSI_out)))

scale_max <- 1.5
scale_min <- -1.5
LSI_out[LSI_out > scale_max] <- scale_max
LSI_out[LSI_out < scale_min] <- scale_min

hmcols <- colorpanel(100, "steelblue", "white", "tomato")
jpeg(outheatmap, width=4, height=6, units="in", res=600,type="cairo")
heatmap.2(LSI_out, 
          col=hmcols,
          ColSideColors=cell_pal[as.factor(cells_tree_cut)],
          RowSideColors=gene_pal[as.factor(genes_tree_cut)],
          Rowv = as.dendrogram(hclust_genes), Colv = as.dendrogram(hclust_cells),
          labRow=FALSE, labCol=FALSE, trace="none",  scale="none",
          useRaster=TRUE)
dev.off()

pdf(paste0(outheatmap,".newcolors.pdf"), width=4, height=6)
heatmap.2(LSI_out, 
          col=hmcols,
          ColSideColors=cell_pal[as.factor(cells_tree_cut)],
#          RowSideColors=gene_pal[as.factor(genes_tree_cut)],
          Rowv = as.dendrogram(hclust_genes), Colv = as.dendrogram(hclust_cells),
          labRow=FALSE, labCol=FALSE, trace="none",  scale="none",
          useRaster=TRUE)
dev.off()

inder = match(oldcells[,1],colnames(LSI_out))
inderna = which(is.na(inder))
oldcols = rep("red",times=length(colnames(LSI_out)))
oldcols[inder[-inderna]] = brewer.pal(5, "Paired")[unlist(as.factor(oldcells[-inderna,2]))] 

jpeg(paste0(outheatmap,"oldcolors.jpg"), width=4, height=6, units="in", res=600,type="cairo")
heatmap.2(LSI_out, 
          col=hmcols,
          ColSideColors=oldcols,
          #RowSideColors=gene_pal[as.factor(genes_tree_cut)],
          Rowv = as.dendrogram(hclust_genes), Colv = as.dendrogram(hclust_cells),
          labRow=FALSE, labCol=FALSE, trace="none",  scale="none",
          useRaster=TRUE)
dev.off()


pdf(compplotsout)
par(mfrow=c(3,3))
plot(LLL$dk[,1],log10(colSums(ncounts)),pch=20,col="dodgerblue2",xlab="PC1",ylab="Read Depth")
for(i in seq(1,dim(LLL$dk)[2],2)){
  plot(LLL$dk[,i],LLL$dk[,i+1],pch=20,col="mediumseagreen",xlab=i,ylab=i+1)
}
dev.off()

sitebeds = strsplit2(annot.ncounts,"-|:")
site_clusters = cbind(annot.ncounts,genes_tree_cut)
site_clusters = cbind(sitebeds,site_clusters)
colnames(site_clusters) = c("Chrom","Start","End","Site","Cluster")
cell_clusters = cbind(colnames(ncounts),cells_tree_cut)
colnames(cell_clusters) = c("Cell","Cluster")
write.table(site_clusters,sitesout,row.names=F,col.names=F,sep="\t",quote=F)
write.table(cell_clusters,cellsout,row.names=F,col.names=F,sep="\t",quote=F)


####6to8hrs
numclust = 4
inmat = "/net/shendure/vol10/projects/scATAC/nobackup/FlyEmbryogenesis/sc_drone_tha/bowtie/SCatac_tha.bowtie.6to8.fixed.nodups.2kbwindows.dhsmatrix.txt"
outhists = paste0("/net/shendure/vol10/projects/scATAC/nobackup/FlyEmbryogenesis/sc_drone_tha/bowtie/SCatac_tha.bowtie.6to8.fixed.nodups.2kbwindows.",cellcutoff,"cells.",sitecutoff,"sites.",numcomponents,"components.",numclust,"clusters.binaryhists.pdf")
outheatmap = paste0("/net/shendure/vol10/projects/scATAC/nobackup/FlyEmbryogenesis/sc_drone_tha/bowtie/SCatac_tha.bowtie.6to8.fixed.nodups.2kbwindows.",cellcutoff,"cells.",sitecutoff,"sites.",numcomponents,"components.",numclust,"clusters.heatmap.jpg")
sitesout = paste0("/net/shendure/vol10/projects/scATAC/nobackup/FlyEmbryogenesis/sc_drone_tha/bowtie/SCatac_tha.bowtie.6to8.fixed.nodups.2kbwindows.",cellcutoff,"cells.",sitecutoff,"sites.",numcomponents,"components.",numclust,"clusters.siteclusters.txt")
cellsout = paste0("/net/shendure/vol10/projects/scATAC/nobackup/FlyEmbryogenesis/sc_drone_tha/bowtie/SCatac_tha.bowtie.6to8.fixed.nodups.2kbwindows.",cellcutoff,"cells.",sitecutoff,"sites.",numcomponents,"components.",numclust,"clusters.cellclusters.txt")
tsneout = paste0("/net/shendure/vol10/projects/scATAC/nobackup/FlyEmbryogenesis/sc_drone_tha/bowtie/SCatac_tha.bowtie.6to8.fixed.nodups.2kbwindows.",cellcutoff,"cells.",sitecutoff,"sites.",numcomponents,"components.",numclust,"clusters.tsne.pdf")
compplotsout = paste0("/net/shendure/vol10/projects/scATAC/nobackup/FlyEmbryogenesis/sc_drone_tha/bowtie/SCatac_tha.bowtie.6to8.fixed.nodups.2kbwindows.",cellcutoff,"cells.",sitecutoff,"sites.",numcomponents,"components.",numclust,"clusters.componentplots.pdf")
oldcells = read.table("/net/shendure/vol10/projects/scATAC/nobackup/FlyEmbryogenesis/sc_drone_tha/sorting/heatmaps/SCatac_tha.6to8.true.nodups.2kbwindows.0.1cells.20000sites.6components.4clusters.cellclusters.txt")
weirdcolors = c("#1F78B4","#FFD700","#60CC52","#E31A1C")

if(gz){
  outprefix = strsplit(inmat,"[.]dhsmatrix[.]txt[.]gz")[[1]]
}else{
  outprefix = strsplit(inmat,"[.]dhsmatrix[.]txt")[[1]]
}
rdsfile = paste0(outprefix,".dhsmatrix.rds")
sparserdsfile = paste0(outprefix,".dhsmatrix.binary.sparse.rds")

if(file.exists(sparserdsfile)){
  bigmat.bin = readRDS(sparserdsfile)
  namelist = readRDS(paste0(outprefix,".dhsmatrix.binary.namelist.rds"))
  label = readRDS(paste0(outprefix,".dhsmatrix.binary.sitelist.rds"))
}else{
  #Read in matrix
  if(file.exists(rdsfile)){
    bigmat = readRDS(rdsfile)
  }else{
    if(gz == TRUE){
      bigmat = fread(paste0("zcat ",inmat),header=T)
    } else {
      bigmat = fread(inmat,header=T)
    }
    bigmat = as.data.frame(bigmat)
    saveRDS(bigmat,rdsfile)
  }
  print("Formatting data...")
  sortmat = as.matrix(bigmat[,-c(1:4)])
  #Make a list of sites
  sites = strsplit2(bigmat[,4],"_")
  label = paste0(sites[,1],":",sites[,2],"-",sites[,3])
  currnonzs = unlist(apply(sortmat,2,function(x) which(x > 0)))
  namers = gsub('[[:digit:]]+', '', names(currnonzs))
  namenums = c(1:length(unique(namers)))[unlist(as.factor(namers))]
  mastermat = cbind(currnonzs,namenums)
  namelist = unique(namers)
  bigmat.bin = sparseMatrix(i=mastermat[,1],j=mastermat[,2],x=rep(1,times=dim(mastermat)[1]),dims=dim(sortmat))
  saveRDS(bigmat.bin,sparserdsfile)
  saveRDS(namelist,paste0(outprefix,".dhsmatrix.binary.namelist.rds"))
  saveRDS(label,paste0(outprefix,".dhsmatrix.binary.sitelist.rds"))
  rm(currnonzs)
  rm(mastermat)
  rm(namers)
  rm(namenums)
  gc()
}

annot = label
colnames(bigmat.bin) = namelist
rownames(bigmat.bin) = annot

num_cells_ncounted = rowSums(bigmat.bin)
num_sites_ncounted = colSums(bigmat.bin)

pdf(outhists)
hist(log10(num_cells_ncounted),main="Number of Cells Each Site is Observed In",breaks=50)
abline(v=log10(num_cells_ncounted[order(num_cells_ncounted,decreasing=T)[sitecutoff]]),lwd=2,col="indianred")
hist(log10(num_sites_ncounted),main="Number of Sites Each Cell Uses",breaks=50)
abline(v=log10(quantile(num_sites_ncounted,probs=cellcutoff)),lwd=2,col="indianred")
annot.ncounts = annot[which(num_cells_ncounted >= num_cells_ncounted[order(num_cells_ncounted,decreasing=T)[sitecutoff]])]
ncounts = bigmat.bin[which(num_cells_ncounted >= num_cells_ncounted[order(num_cells_ncounted,decreasing=T)[sitecutoff]]),]
new_counts = colSums(ncounts)
ncounts = ncounts[,new_counts >= quantile(new_counts,probs=cellcutoff)]
hist(log10(new_counts),main="Number of Sites Each Cell Uses",breaks=50)
abline(v=log10(quantile(new_counts,probs=cellcutoff)),lwd=2,col="indianred")
dev.off()

annot.ncounts = annot.ncounts[rowSums(ncounts) > 0]
ncounts = ncounts[rowSums(ncounts) > 0,]
ncounts = Matrix(ncounts,sparse=T)

nfreqs <- t(t(ncounts) / colSums(ncounts))
tf_idf_counts <- nfreqs * log(1 + ncol(ncounts) / rowSums(ncounts))
set.seed(0)
LLL <- my_lsa(tf_idf_counts,numcomponents)
sk_diag <- matrix(0, nrow=length(LLL$sk), ncol=length(LLL$sk))
diag(sk_diag) <- LLL$sk
sk_diag[1,1] = 0
LLL_d <- t(sk_diag %*% t(LLL$dk))
hclust_cells <- hclust(proxy::dist(LLL_d, method="cosine"), method="ward.D2")
TTT_d <- t(sk_diag %*% t(LLL$tk))
hclust_genes <- hclust(proxy::dist(TTT_d, method="cosine"), method="ward.D2")

genes_tree_cut <- cutree(hclust_genes, numclust)
gene_pal <- weirdcolors
cells_tree_cut <- cutree(hclust_cells, numclust)
cell_pal <- weirdcolors

LSI_out <-  t(t(sk_diag %*% t(LLL$dk)) %*% t(LLL$tk))
LSI_out <- t(scale(t(LSI_out)))

scale_max <- 1.5
scale_min <- -1.5
LSI_out[LSI_out > scale_max] <- scale_max
LSI_out[LSI_out < scale_min] <- scale_min

hmcols <- colorpanel(100, "steelblue", "white", "tomato")
jpeg(outheatmap, width=4, height=6, units="in", res=600,type="cairo")
heatmap.2(LSI_out, 
          col=hmcols,
          ColSideColors=cell_pal[as.factor(cells_tree_cut)],
          RowSideColors=gene_pal[as.factor(genes_tree_cut)],
          Rowv = as.dendrogram(hclust_genes), Colv = as.dendrogram(hclust_cells),
          labRow=FALSE, labCol=FALSE, trace="none",  scale="none",
          useRaster=TRUE)
dev.off()

pdf(paste0(outheatmap,".newcolors.pdf"), width=4, height=6)
heatmap.2(LSI_out, 
          col=hmcols,
          ColSideColors=cell_pal[as.factor(cells_tree_cut)],
          #          RowSideColors=gene_pal[as.factor(genes_tree_cut)],
          Rowv = as.dendrogram(hclust_genes), Colv = as.dendrogram(hclust_cells),
          labRow=FALSE, labCol=FALSE, trace="none",  scale="none",
          useRaster=TRUE)
dev.off()

inder = match(oldcells[,1],colnames(LSI_out))
inderna = which(is.na(inder))
oldcols = rep("red",times=length(colnames(LSI_out)))
oldcols[inder[-inderna]] = brewer.pal(5, "Paired")[unlist(as.factor(oldcells[-inderna,2]))] 

jpeg(paste0(outheatmap,"oldcolors.jpg"), width=4, height=6, units="in", res=600,type="cairo")
heatmap.2(LSI_out, 
          col=hmcols,
          ColSideColors=oldcols,
          #RowSideColors=gene_pal[as.factor(genes_tree_cut)],
          Rowv = as.dendrogram(hclust_genes), Colv = as.dendrogram(hclust_cells),
          labRow=FALSE, labCol=FALSE, trace="none",  scale="none",
          useRaster=TRUE)
dev.off()


pdf(compplotsout)
par(mfrow=c(3,3))
plot(LLL$dk[,1],log10(colSums(ncounts)),pch=20,col="dodgerblue2",xlab="PC1",ylab="Read Depth")
for(i in seq(1,dim(LLL$dk)[2],2)){
  plot(LLL$dk[,i],LLL$dk[,i+1],pch=20,col="mediumseagreen",xlab=i,ylab=i+1)
}
dev.off()

sitebeds = strsplit2(annot.ncounts,"-|:")
site_clusters = cbind(annot.ncounts,genes_tree_cut)
site_clusters = cbind(sitebeds,site_clusters)
colnames(site_clusters) = c("Chrom","Start","End","Site","Cluster")
cell_clusters = cbind(colnames(ncounts),cells_tree_cut)
colnames(cell_clusters) = c("Cell","Cluster")
write.table(site_clusters,sitesout,row.names=F,col.names=F,sep="\t",quote=F)
write.table(cell_clusters,cellsout,row.names=F,col.names=F,sep="\t",quote=F)


####10to12hrs
numclust = 4
inmat = "/net/shendure/vol10/projects/scATAC/nobackup/FlyEmbryogenesis/sc_drone_tha/bowtie/SCatac_tha.bowtie.10to12.fixed.nodups.2kbwindows.dhsmatrix.txt"
outhists = paste0("/net/shendure/vol10/projects/scATAC/nobackup/FlyEmbryogenesis/sc_drone_tha/bowtie/SCatac_tha.bowtie.10to12.fixed.nodups.2kbwindows.",cellcutoff,"cells.",sitecutoff,"sites.",numcomponents,"components.",numclust,"clusters.binaryhists.pdf")
outheatmap = paste0("/net/shendure/vol10/projects/scATAC/nobackup/FlyEmbryogenesis/sc_drone_tha/bowtie/SCatac_tha.bowtie.10to12.fixed.nodups.2kbwindows.",cellcutoff,"cells.",sitecutoff,"sites.",numcomponents,"components.",numclust,"clusters.heatmap.jpg")
sitesout = paste0("/net/shendure/vol10/projects/scATAC/nobackup/FlyEmbryogenesis/sc_drone_tha/bowtie/SCatac_tha.bowtie.10to12.fixed.nodups.2kbwindows.",cellcutoff,"cells.",sitecutoff,"sites.",numcomponents,"components.",numclust,"clusters.siteclusters.txt")
cellsout = paste0("/net/shendure/vol10/projects/scATAC/nobackup/FlyEmbryogenesis/sc_drone_tha/bowtie/SCatac_tha.bowtie.10to12.fixed.nodups.2kbwindows.",cellcutoff,"cells.",sitecutoff,"sites.",numcomponents,"components.",numclust,"clusters.cellclusters.txt")
tsneout = paste0("/net/shendure/vol10/projects/scATAC/nobackup/FlyEmbryogenesis/sc_drone_tha/bowtie/SCatac_tha.bowtie.10to12.fixed.nodups.2kbwindows.",cellcutoff,"cells.",sitecutoff,"sites.",numcomponents,"components.",numclust,"clusters.tsne.pdf")
compplotsout = paste0("/net/shendure/vol10/projects/scATAC/nobackup/FlyEmbryogenesis/sc_drone_tha/bowtie/SCatac_tha.bowtie.10to12.fixed.nodups.2kbwindows.",cellcutoff,"cells.",sitecutoff,"sites.",numcomponents,"components.",numclust,"clusters.componentplots.pdf")
oldcells = read.table("/net/shendure/vol10/projects/scATAC/nobackup/FlyEmbryogenesis/sc_drone_tha/sorting/heatmaps/SCatac_tha.10to12.true.nodups.2kbwindows.0.1cells.20000sites.6components.4clusters.cellclusters.txt")
weirdcolors = c("#1F78B4","#60CC52","#FFD700","#E31A1C")

if(gz){
  outprefix = strsplit(inmat,"[.]dhsmatrix[.]txt[.]gz")[[1]]
}else{
  outprefix = strsplit(inmat,"[.]dhsmatrix[.]txt")[[1]]
}
rdsfile = paste0(outprefix,".dhsmatrix.rds")
sparserdsfile = paste0(outprefix,".dhsmatrix.binary.sparse.rds")

if(file.exists(sparserdsfile)){
  bigmat.bin = readRDS(sparserdsfile)
  namelist = readRDS(paste0(outprefix,".dhsmatrix.binary.namelist.rds"))
  label = readRDS(paste0(outprefix,".dhsmatrix.binary.sitelist.rds"))
}else{
  #Read in matrix
  if(file.exists(rdsfile)){
    bigmat = readRDS(rdsfile)
  }else{
    if(gz == TRUE){
      bigmat = fread(paste0("zcat ",inmat),header=T)
    } else {
      bigmat = fread(inmat,header=T)
    }
    bigmat = as.data.frame(bigmat)
    saveRDS(bigmat,rdsfile)
  }
  print("Formatting data...")
  sortmat = as.matrix(bigmat[,-c(1:4)])
  #Make a list of sites
  sites = strsplit2(bigmat[,4],"_")
  label = paste0(sites[,1],":",sites[,2],"-",sites[,3])
  currnonzs = unlist(apply(sortmat,2,function(x) which(x > 0)))
  namers = gsub('[[:digit:]]+', '', names(currnonzs))
  namenums = c(1:length(unique(namers)))[unlist(as.factor(namers))]
  mastermat = cbind(currnonzs,namenums)
  namelist = unique(namers)
  bigmat.bin = sparseMatrix(i=mastermat[,1],j=mastermat[,2],x=rep(1,times=dim(mastermat)[1]),dims=dim(sortmat))
  saveRDS(bigmat.bin,sparserdsfile)
  saveRDS(namelist,paste0(outprefix,".dhsmatrix.binary.namelist.rds"))
  saveRDS(label,paste0(outprefix,".dhsmatrix.binary.sitelist.rds"))
  rm(currnonzs)
  rm(mastermat)
  rm(namers)
  rm(namenums)
  gc()
}

annot = label
colnames(bigmat.bin) = namelist
rownames(bigmat.bin) = annot

num_cells_ncounted = rowSums(bigmat.bin)
num_sites_ncounted = colSums(bigmat.bin)

pdf(outhists)
hist(log10(num_cells_ncounted),main="Number of Cells Each Site is Observed In",breaks=50)
abline(v=log10(num_cells_ncounted[order(num_cells_ncounted,decreasing=T)[sitecutoff]]),lwd=2,col="indianred")
hist(log10(num_sites_ncounted),main="Number of Sites Each Cell Uses",breaks=50)
abline(v=log10(quantile(num_sites_ncounted,probs=cellcutoff)),lwd=2,col="indianred")
annot.ncounts = annot[which(num_cells_ncounted >= num_cells_ncounted[order(num_cells_ncounted,decreasing=T)[sitecutoff]])]
ncounts = bigmat.bin[which(num_cells_ncounted >= num_cells_ncounted[order(num_cells_ncounted,decreasing=T)[sitecutoff]]),]
new_counts = colSums(ncounts)
ncounts = ncounts[,new_counts >= quantile(new_counts,probs=cellcutoff)]
hist(log10(new_counts),main="Number of Sites Each Cell Uses",breaks=50)
abline(v=log10(quantile(new_counts,probs=cellcutoff)),lwd=2,col="indianred")
dev.off()

annot.ncounts = annot.ncounts[rowSums(ncounts) > 0]
ncounts = ncounts[rowSums(ncounts) > 0,]
ncounts = Matrix(ncounts,sparse=T)

nfreqs <- t(t(ncounts) / colSums(ncounts))
tf_idf_counts <- nfreqs * log(1 + ncol(ncounts) / rowSums(ncounts))
set.seed(0)
LLL <- my_lsa(tf_idf_counts,numcomponents)
sk_diag <- matrix(0, nrow=length(LLL$sk), ncol=length(LLL$sk))
diag(sk_diag) <- LLL$sk
sk_diag[1,1] = 0
LLL_d <- t(sk_diag %*% t(LLL$dk))
hclust_cells <- hclust(proxy::dist(LLL_d, method="cosine"), method="ward.D2")
TTT_d <- t(sk_diag %*% t(LLL$tk))
hclust_genes <- hclust(proxy::dist(TTT_d, method="cosine"), method="ward.D2")

genes_tree_cut <- cutree(hclust_genes, numclust)
gene_pal <- weirdcolors
cells_tree_cut <- cutree(hclust_cells, numclust)
cell_pal <- weirdcolors

LSI_out <-  t(t(sk_diag %*% t(LLL$dk)) %*% t(LLL$tk))
LSI_out <- t(scale(t(LSI_out)))

scale_max <- 1.5
scale_min <- -1.5
LSI_out[LSI_out > scale_max] <- scale_max
LSI_out[LSI_out < scale_min] <- scale_min

hmcols <- colorpanel(100, "steelblue", "white", "tomato")
jpeg(outheatmap, width=4, height=6, units="in", res=600,type="cairo")
heatmap.2(LSI_out, 
          col=hmcols,
          ColSideColors=cell_pal[as.factor(cells_tree_cut)],
          RowSideColors=gene_pal[as.factor(genes_tree_cut)],
          Rowv = as.dendrogram(hclust_genes), Colv = as.dendrogram(hclust_cells),
          labRow=FALSE, labCol=FALSE, trace="none",  scale="none",
          useRaster=TRUE)
dev.off()

pdf(paste0(outheatmap,".newcolors.pdf"), width=4, height=6)
heatmap.2(LSI_out, 
          col=hmcols,
          ColSideColors=cell_pal[as.factor(cells_tree_cut)],
          #          RowSideColors=gene_pal[as.factor(genes_tree_cut)],
          Rowv = as.dendrogram(hclust_genes), Colv = as.dendrogram(hclust_cells),
          labRow=FALSE, labCol=FALSE, trace="none",  scale="none",
          useRaster=TRUE)
dev.off()

inder = match(oldcells[,1],colnames(LSI_out))
inderna = which(is.na(inder))
oldcols = rep("red",times=length(colnames(LSI_out)))
oldcols[inder[-inderna]] = brewer.pal(5, "Paired")[unlist(as.factor(oldcells[-inderna,2]))] 

jpeg(paste0(outheatmap,"oldcolors.jpg"), width=4, height=6, units="in", res=600,type="cairo")
heatmap.2(LSI_out, 
          col=hmcols,
          ColSideColors=oldcols,
          #RowSideColors=gene_pal[as.factor(genes_tree_cut)],
          Rowv = as.dendrogram(hclust_genes), Colv = as.dendrogram(hclust_cells),
          labRow=FALSE, labCol=FALSE, trace="none",  scale="none",
          useRaster=TRUE)
dev.off()


pdf(compplotsout)
par(mfrow=c(3,3))
plot(LLL$dk[,1],log10(colSums(ncounts)),pch=20,col="dodgerblue2",xlab="PC1",ylab="Read Depth")
for(i in seq(1,dim(LLL$dk)[2],2)){
  plot(LLL$dk[,i],LLL$dk[,i+1],pch=20,col="mediumseagreen",xlab=i,ylab=i+1)
}
dev.off()

sitebeds = strsplit2(annot.ncounts,"-|:")
site_clusters = cbind(annot.ncounts,genes_tree_cut)
site_clusters = cbind(sitebeds,site_clusters)
colnames(site_clusters) = c("Chrom","Start","End","Site","Cluster")
cell_clusters = cbind(colnames(ncounts),cells_tree_cut)
colnames(cell_clusters) = c("Cell","Cluster")
write.table(site_clusters,sitesout,row.names=F,col.names=F,sep="\t",quote=F)
write.table(cell_clusters,cellsout,row.names=F,col.names=F,sep="\t",quote=F)
