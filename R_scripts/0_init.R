### collect all datafiles needed for dLGN_cross_species analysis

library(matrixStats)
library(Matrix)
library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)
library(scrattch.hicat)
library(feather)

datlist=list()
metalist=list()

###macaque data###
load("data/20191206_LGN_Mmul_10_star_exon.Rdata")
load("data/20191206_LGN_Mmul_10_star_intron.Rdata")
mac_dat=exon+intron
load("data/20191206_LGN_Mmul_10_star_samp.dat.comb.Rdata")
mac_dat=sweep(mac_dat,2,colSums(mac_dat),"/")*10^6
mac_dat2=Matrix(as.matrix(mac_dat))
rownames(mac_dat2)=rownames(mac_dat)
colnames(mac_dat2)=colnames(mac_dat)
mac_dat2@x=mac_dat2@x/rep.int(Matrix::colSums(mac_dat2),diff(mac_dat2@p))*10^6
samp.dat$donor=do.call(rbind,strsplit(samp.dat$external_donor_name,"\\."))[,1]
datlist[["macaque"]]=mac_dat2
metalist[["macaque"]]=data.frame(samp.dat)


###human###
human_dat=read_feather("data/human_data.feather")
human_meta=read_feather("data/human_anno.feather")
human_dat2=Matrix(as.matrix(human_dat[,-ncol(human_dat)]))
rownames(human_dat2)=human_dat[,ncol(human_dat)][[1]]
human_dat2@x=human_dat2@x/rep.int(Matrix::colSums(human_dat2), diff(human_dat2@p))*10^6
metadat=as.matrix(human_meta)
rownames(metadat)=metadat[,1]
datlist[["human"]]=t(human_dat2)
metalist[["human"]]=data.frame(metadat)

###mouse###
load("data/20180105_RSC-005-049_mouse_LGN_star_cpm.Rdata")
lgn_cpmR=cpmR
load("data/20180105_RSC-005-049_mouse_LGN_star_samp.dat.Rdata")
load("data/LP_ss2v4_RSC191.RDA")
lp_cpmR=cpmR
rm(cpmR)
common.cols=intersect(colnames(samp.dat),colnames(samp.dat.lp))
metadat=rbind(samp.dat[,common.cols],samp.dat.lp[,common.cols])
metadat$roi=gsub("^LGv","Mouse LGv",metadat$roi)
metadat$roi=gsub("^LP","Mouse LP",metadat$roi)
datlist[["mouse"]]=cbind(lgn_cpmR,lp_cpmR)
metalist[["mouse"]]=data.frame(metadat)

save(datlist,metalist,file="int_data/collected_data_20200818.rda")




