#module load gmp/5.0.2 mpfr/3.1.0 mpc/0.8.2 gcc/4.9.1 R/3.2.1
library(Matrix)
library(limma)
library(RColorBrewer)
library(proxy)
library(irlba,lib.loc="~/R/x86_64-unknown-linux-gnu-library/3.1/")
library(gplots)
library(data.table)
library(methods)
library(argparse)

parser = argparse::ArgumentParser(description="Script to generate LSI-based heatmaps from sci-ATAC-seq data.")
parser$add_argument('-I','--inmat', help='Site by cell matrix of accessibility in individual cells.')
parser$add_argument('-T','--outhists', help='Name of output pdf for plotting number of sites per cell and number of cells per site.')
parser$add_argument('-O','--outheatmap', help='Name of output jpeg heatmap.')
parser$add_argument('-S','--outsites', help='Name of output file for site clusters (each site gets assigned to one of `numclusters` clusters.')
parser$add_argument('-C','--outcells', help='Name of output file for site clusters (each cell gets assigned to one of `numclusters` clusters.')
parser$add_argument('-P','--outcomps', help='Name of output pdf for individual component plots.')
parser$add_argument('-Z','--gzip', action="store_true", help='Boolean for whether "inmat" is gzipped.')
parser$add_argument('-N','--numclusters', default='4', help='Number of clusters to cut the dendrogram into.')
parser$add_argument('-M','--numcomps', default='6', help='Number of components to include in LSI (the first component will be skipped).')
parser$add_argument('-CC','--cellcutoff', default='0.1', help='Fraction of cells to filter out as low coverage cells.')
parser$add_argument('-SC','--sitecutoff', default='20000', help='Number of sites to include (all sites tied at the threshold will be included).')
args = parser$parse_args()


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

if(args$gzip){
  outprefix = strsplit(args$inmat,"[.]matrix[.]txt[.]gz")[[1]]
}else{
  outprefix = strsplit(args$inmat,"[.]matrix[.]txt")[[1]]
}
rdsfile = paste0(outprefix,".matrix.rds")
sparserdsfile = paste0(outprefix,".matrix.binary.sparse.rds")

if(file.exists(sparserdsfile)){
  bigmat.bin = readRDS(sparserdsfile)
  namelist = readRDS(paste0(outprefix,".matrix.binary.namelist.rds"))
  label = readRDS(paste0(outprefix,".matrix.binary.sitelist.rds"))
  annot = rownames(label)
}else{
  #Read in matrix
  if(file.exists(rdsfile)){
    bigmat = readRDS(rdsfile)
  }else{
    if(args$gzip == TRUE){
      bigmat = fread(paste0("zcat ",args$inmat),header=T)
    } else {
      bigmat = fread(args$inmat,header=T)
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
  saveRDS(namelist,paste0(outprefix,".matrix.binary.namelist.rds"))
  saveRDS(label,paste0(outprefix,".matrix.binary.sitelist.rds"))
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

pdf(args$outhists)
hist(log10(num_cells_ncounted),main="Number of Cells Each Site is Observed In",breaks=50)
abline(v=log10(num_cells_ncounted[order(num_cells_ncounted,decreasing=T)[args$sitecutoff]]),lwd=2,col="indianred")
hist(log10(num_sites_ncounted),main="Number of Sites Each Cell Uses",breaks=50)
abline(v=log10(quantile(num_sites_ncounted,probs=args$cellcutoff)),lwd=2,col="indianred")
annot.ncounts = annot[which(num_cells_ncounted >= num_cells_ncounted[order(num_cells_ncounted,decreasing=T)[args$sitecutoff]])]
ncounts = bigmat.bin[which(num_cells_ncounted >= num_cells_ncounted[order(num_cells_ncounted,decreasing=T)[args$sitecutoff]]),]
new_counts = colSums(ncounts)
ncounts = ncounts[,new_counts >= quantile(new_counts,probs=args$cellcutoff)]
hist(log10(new_counts),main="Number of Sites Each Cell Uses",breaks=50)
abline(v=log10(quantile(new_counts,probs=args$cellcutoff)),lwd=2,col="indianred")
dev.off()

annot.ncounts = annot.ncounts[rowSums(ncounts) > 0]
ncounts = ncounts[rowSums(ncounts) > 0,]
ncounts = Matrix(ncounts,sparse=T)

nfreqs <- t(t(ncounts) / colSums(ncounts))
tf_idf_counts <- nfreqs * log(1 + ncol(ncounts) / rowSums(ncounts))
set.seed(0)
LLL <- my_lsa(tf_idf_counts,args$numcomps)
sk_diag <- matrix(0, nrow=length(LLL$sk), ncol=length(LLL$sk))
diag(sk_diag) <- LLL$sk
sk_diag[1,1] = 0
LLL_d <- t(sk_diag %*% t(LLL$dk))
hclust_cells <- hclust(proxy::dist(LLL_d, method="cosine"), method="ward.D2")
TTT_d <- t(sk_diag %*% t(LLL$tk))
hclust_genes <- hclust(proxy::dist(TTT_d, method="cosine"), method="ward.D2")

genes_tree_cut <- cutree(hclust_genes, args$numclusters)
gene_pal <- brewer.pal(args$numclusters, "Paired")
cells_tree_cut <- cutree(hclust_cells, args$numclusters)
cell_pal <- brewer.pal(args$numclusters, "Paired")

LSI_out <-  t(t(sk_diag %*% t(LLL$dk)) %*% t(LLL$tk))
LSI_out <- t(scale(t(LSI_out)))

scale_max <- 1.5
scale_min <- -1.5
LSI_out[LSI_out > scale_max] <- scale_max
LSI_out[LSI_out < scale_min] <- scale_min

hmcols <- colorpanel(100, "steelblue", "white", "tomato")
jpeg(args$outheatmap, width=4, height=6, units="in", res=600,type="cairo")
heatmap.2(LSI_out, 
          col=hmcols,
          ColSideColors=cell_pal[as.factor(cells_tree_cut)],
#          RowSideColors=gene_pal[as.factor(genes_tree_cut)],
          Rowv = as.dendrogram(hclust_genes), Colv = as.dendrogram(hclust_cells),
          labRow=FALSE, labCol=FALSE, trace="none",  scale="none",
          useRaster=TRUE)
dev.off()

pdf(args$outcomps)
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
write.table(site_clusters,args$outsites,row.names=F,col.names=F,sep="\t",quote=F)
write.table(cell_clusters,args$outcells,row.names=F,col.names=F,sep="\t",quote=F)
