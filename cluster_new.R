library(SeqArray)
library(gdsfmt)
library(SNPRelate)
library(stringr)
library(dendextend)
library(phylogram)
library(ggplot2)

args <- commandArgs(TRUE)
mod <- read.table(args[3])

if (mod == "chunk_mod") {
print('chunk_mod')
file <- list.files(path = args[1], pattern = ".vcf", full.names=TRUE)
for (w in file) {
print(w)
snpgdsVCF2GDS(w,  paste0(w, ".gds"), method="biallelic.only")
}

fn <- list.files(path = args[1], pattern = ".gds", full.names=TRUE)
snpgdsCombineGeno(fn, args[2])
} else {
print('simple_mod')
vcf.fn <- args[4]
snpgdsVCF2GDS(vcf.fn,  args[2], method="biallelic.only")
}

genofile <- snpgdsOpen(args[2])
print(genofile)

set.seed(1000)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.8)
str(snpset)
snpset.id <- unlist(unname(snpset))
head(snpset.id)
samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
samp.id = str_split(samp.id, "/", simplify = T)[,2]
samp.id = str_split(samp.id, ".fastq", simplify = T)[,1]

pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=40)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

tab <- data.frame(sample.id = samp.id,
    EV1 = pca$eigenvect[,1],    # the first eigenvector
    EV2 = pca$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)
write.table(tab, file=paste('results/cluster_tab.tsv', sep=""), sep="\t")

png("results/cluster.png", width=15,height=6,units="in",res=300)
ggplot(tab, aes(x=EV2, y=EV1)) + geom_point() + geom_text(aes(label=sample.id))


ibs.hc<-snpgdsHCluster(snpgdsIBS(genofile,num.thread=2, autosome.only=FALSE))
rv <- snpgdsCutTree(ibs.hc)
write.dendrogram(rv$dendrogram, file='results/dendrogram.tree', edges = TRUE)
png("results/plotdendogram.png",width=800,height=2500)
par(cex=1,font=20)
rv$dendrogram %>% set("labels", c(samp.id)) %>% set("labels_cex", 0.9) %>% plot(horiz = TRUE)

