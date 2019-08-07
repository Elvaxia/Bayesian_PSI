
setwd('/mnt/data8/xiaoxia/bayes_bulk/3bias/alpha0.8/')

library(data.table)
library(R2jags)
library(dplyr)
load('/mnt/raid61/Personal_data/tangchao/IR6/result/polyester_Simulation/chr1_depth_mat_and_BinExp_3bias.RData')
#genes with only single isoform
#SingleIsoGenes <- all_exp%>% group_by(gene_id)%>%count(freq)%>%count(gene_id)%>%filter(n<2)%>%pull(gene_id)
#Exp.SIgenes <- all_exp[gene_id%in%SingleIsoGenes & region !='SJ',]
#use genes with only single isoform to estimate the GBC (global bias curve)


#for( g in unique(Exp.SIgenes$tx_id)){
  #those genes may have multiple unoverlappled isoforms,so just use isoform of those genes
 # iso <- Exp.SIgenes[tx_id==g,]
  #setkey(iso,start)
  #if(unique(iso$strand) == '+'){
   # seq <- iso$ExpSum
  #}else{
  #  seq <- rev(iso$ExpSum)
  #}
  #normalize the seq to 0-1
  #seq_norm <- (seq-min(seq))/(max(seq)-min(seq)) 
   # res <- data.table(bin = cut(seq(1,length(seq)), 10, include.lowest=F, labels=seq(1,10)),
    #             reads = seq) %>% group_by(bin) %>%  
    #              summarise(combined = sum(reads))%>% pull(combined)
  #bined.reads[[g]] <- res
  #normed.reads[[g]] <- 2*(res-min(res))/(max(res)-min(res))
#}

Iter.sum <- function(vector){
  if(length(vector)==1){
    vector.new <- vector
  }else{
    vector.new <- c()
    for(i in 2:length(vector)){
      vector.new[1] <- vector[1]
      vector.new[i] <- sum(vector[1:i])
    }
  }
  return(vector.new)
}

normalize.seq <- function(seq){
  if(length(seq)==1){
    normalized <- 1
  }else{
    normalized <- (seq-min(seq))/(max(seq)-min(seq))
  }
  return(normalized)
}

excluded.genes <- unique(all_exp[all_exp$freq>2,]$gene_id) 
expres.low <- all_exp[, sum(ExpSum),by=gene_id] %>% filter(V1 <1000)%>% pull(gene_id)
excluded.genes <- unique(c(excluded.genes, expres.low))
all_exp <- all_exp[!all_exp$gene_id%in%excluded.genes,]
tmp <- all_exp[gene_id=="ENSG00000007933.13_4",]
setkey(tmp, start)
head(tmp)
#use gene body coverge to be the global bias curve 
GBC <- read.table('/mnt/data8/xiaoxia/bayes_bulk/GeneBody/3bias.geneBodyCoverage.txt', row.names = 1)
GBC <- as.numeric(GBC[2,])
normalized.GBC <- (GBC-min(GBC))/(max(GBC)-min(GBC))
#LBC ( loca bias curve) depict gene-specific read distibutions



bayes.mod.exon<-function(){
    for (j in 1:nB) {

      p[j] <- Index[j,1]*theta*L[j]/L1+Index[j,2]*(1-theta)*L[j]/L2
      N[j] ~ dpois(p[j]*W)
     # N[j] ~ dbinom(p[j],W)#number of molecules we can observe
    }

                theta ~ dunif(0,1)#prior on theta

}
bayes.mod.params<-c("theta")
bayes.mod.inits<-function(){
        list(theta=runif(1,0,1))
}
posterior.iso1 <- list()
for(g in unique(all_exp$gene_id)){
  tmp <- all_exp[gene_id==g,]
  tmp <- tmp[tmp$ExpSum!=0,]
  setkey(tmp, start)
  if(!'ExonBin'%in%tmp$region){
    posterior.iso1[[g]] <- data.frame(isoform=ifelse(nrow(tmp)==0,NA,unique(unlist(strsplit(tmp$tx_id,";|:")))[1]),
                                      inferenced.psi=NA,gene.id =g )
  }else{
 # only one isoform expressed, skip the model
    if(length(unique(unlist(strsplit(tmp$tx_id,";|:"))))  ==1 ){
      posterior.iso1[[g]] <- data.frame(isoform=unique(unlist(strsplit(tmp$tx_id,";|:"))),
                                        inferenced.psi=1,gene.id =g )
    }else{
    index <- matrix(rep(0,2*nrow(tmp)),c(nrow(tmp),2))
    iso1 <- unique(unlist(strsplit(tmp$tx_id,";|:")))[1]
    iso2 <- unique(unlist(strsplit(tmp$tx_id,":|;")))[2]
    index[grep(iso1, tmp$tx_id),1] <- 1
    index[grep(iso2, tmp$tx_id),2] <- 1
    index.exon <- index[tmp$region == 'ExonBin',]######
    if(is.null(nrow(index.exon))){
      index.exon <- matrix(index.exon, c(1,1))#sometimes there is only one bin
    }
    L1 <- sum(tmp[index[,1]==1 & region == 'ExonBin',]$width)
    L2 <- sum(tmp[index[,2]==1 & region == 'ExonBin',]$width)
    
    W1 <- sum(index[,1]*tmp$ExpSum)
    W2 <- sum(index[,2]*tmp$ExpSum)
    
    type <- tmp$region
    type<-gsub('SJ','0',type)
    type<-gsub('ExonBin','1',type)
    type <- as.numeric(type)
    N.exon <- tmp[tmp$region =='ExonBin',]$ExpSum
    ln <- tmp[tmp$region =='ExonBin',]$width
#weight index matrix with sequcencing preference of each bin of each isoform
    if(nrow(index.exon) ==1){
      weighted.index.exon <- index.exon
    }else{
      seq1 <- tmp[index[,1]==1&region == 'ExonBin', ]$ExpSum
      seq2 <- tmp[index[,2]==1&region == 'ExonBin', ]$ExpSum
      #normalize the seq to 0-1
      seq1_norm <- normalize.seq(seq1+rnorm(length(seq1)))
      seq2_norm <- normalize.seq(seq2+rnorm(length(seq2)))
      LBC.weighted.index.exon <- index.exon
      LBC.weighted.index.exon[LBC.weighted.index.exon[,1] ==1,1] <- seq1_norm
      LBC.weighted.index.exon[LBC.weighted.index.exon[,2] ==1,2] <- seq2_norm
      
      #GBC weighted.index.exon
      if(unique(tmp$strand=='+')){
        L1.w <- tmp[index[,1]==1&region == 'ExonBin', ]$width
        L2.w <- tmp[index[,2]==1&region == 'ExonBin', ]$width
        L1.w.s <- Iter.sum(L1.w)
        qtiles.1 <- as.numeric(quantile(seq(1,sum(L1.w),1),probs = seq(0.1,1,0.1)))
        qtiles.list.1 <- sapply(L1.w.s, function(x)which.min(abs(qtiles.1-x)))
        GBC.weighted.1 <- sapply(qtiles.list.1, function(x)normalized.GBC[x*10])
        #
        L2.w.s <- Iter.sum(L2.w)
        qtiles.2 <- as.numeric(quantile(seq(1,sum(L2.w),1),probs = seq(0.1,1,0.1)))
        qtiles.list.2 <- sapply(L2.w.s, function(x)which.min(abs(qtiles.2-x)))
        GBC.weighted.2 <- sapply(qtiles.list.2, function(x)normalized.GBC[x*10])
        #
        GBC.weighted.index.exon <- index.exon
        GBC.weighted.index.exon[GBC.weighted.index.exon[,1] ==1,1] <- GBC.weighted.1
        GBC.weighted.index.exon[GBC.weighted.index.exon[,2] ==1,2] <- GBC.weighted.2
      }else{
        L1.w <- rev(tmp[index[,1]==1&region == 'ExonBin', ]$width)
        L2.w <- rev(tmp[index[,2]==1&region == 'ExonBin', ]$width)
        L1.w.s <- Iter.sum(L1.w)
        qtiles.1 <- as.numeric(quantile(seq(1,sum(L1.w),1),probs = seq(0.1,1,0.1)))
        qtiles.list.1 <- sapply(L1.w.s, function(x)which.min(abs(qtiles.1-x)))
        GBC.weighted.1 <- sapply(qtiles.list.1, function(x)normalized.GBC[x*10])
        #
        L2.w.s <- Iter.sum(L2.w)
        qtiles.2 <- as.numeric(quantile(seq(1,sum(L2.w),1),probs = seq(0.1,1,0.1)))
        qtiles.list.2 <- sapply(L2.w.s, function(x)which.min(abs(qtiles.2-x)))
        GBC.weighted.2 <- sapply(qtiles.list.2, function(x)normalized.GBC[x*10])
        #
        GBC.weighted.index.exon <- index.exon
        GBC.weighted.index.exon[GBC.weighted.index.exon[,1] ==1,1] <- rev(GBC.weighted.1)
        GBC.weighted.index.exon[GBC.weighted.index.exon[,2] ==1,2] <- rev(GBC.weighted.2)
      }

      #
      alpha <- 0.8
      weighted.index.exon <- GBC.weighted.index.exon*alpha + LBC.weighted.index.exon *(1-alpha)
    }
    dat.jags.exon <- list(Index=weighted.index.exon, nB=nrow(index.exon),
                          N=N.exon, 
                          W = sum(tmp$ExpSum, na.rm = T),L=ln,
                          L1=L1,L2=L2)
    bayes.mod.fit.exon <- jags(data = dat.jags.exon,
                              inits = bayes.mod.inits,
                              parameters.to.save = bayes.mod.params,
                              n.iter = 3000,
                              n.chains = 1, 
                              n.burnin = 1000,
                              model.file = bayes.mod.exon)
    res <- bayes.mod.fit.exon$BUGSoutput$mean$theta
    posterior.iso1[[g]] <- data.frame(isoform=iso1,inferenced.psi=res,gene.id =g )
    }
  }
}
inferenced.psi.exon <- do.call(rbind,posterior.iso1)
##use only SJ
bayes.mod.SJ<-function(){
    for (j in 1:nB) {
      p[j] <- Index[j,1]*theta + Index[j,2]*(1-theta)
      N[j] ~ dpois(p[j]*W)
     # N[j] ~ dbinom(p[j],W)#number of molecules we can observe

    }
                theta ~ dunif(0,1)#prior on theta

}

bayes.mod.params<-c("theta")
bayes.mod.inits<-function(){
        list(theta=runif(1,0,1))
}
posterior.iso1 <- list()
for(g in unique(all_exp$gene_id)){
  tmp <- all_exp[gene_id==g,]
  tmp <- tmp[tmp$ExpSum!=0,]
  setkey(tmp, start)
  #only one isoform expressed, skip the model
  if(!'SJ'%in%tmp$region){
    posterior.iso1[[g]] <- data.frame(isoform=ifelse(nrow(tmp)==0,NA,unique(unlist(strsplit(tmp$tx_id,";|:")))[1]),
                                      inferenced.psi=NA,gene.id =g )
  }else{
      if(length(unique(unlist(strsplit(tmp$tx_id,";|:"))))  ==1 ){
      posterior.iso1[[g]] <- data.frame(isoform=unique(unlist(strsplit(tmp$tx_id,";|:"))),
                                      inferenced.psi=1,gene.id =g )
    }else{
      index <- matrix(rep(0,2*nrow(tmp)),c(nrow(tmp),2))
      iso1 <- unique(unlist(strsplit(tmp$tx_id,";|:")))[1]
      iso2 <- unique(unlist(strsplit(tmp$tx_id,":|;")))[2]
      index[grep(iso1, tmp$tx_id),1] <- 1
      index[grep(iso2, tmp$tx_id),2] <- 1
      index.SJ <- index[tmp$region == 'SJ',]######
    if(is.null(nrow(index.SJ))){
      index.SJ <- matrix(index.SJ, c(1,1))#sometimes there is only one bin
    }
      tmp.new <- tmp
      tmp.new[tmp.new$region=='SJ',]$width <- 0
      ln <- tmp.new$width
      L1 <- sum(tmp.new[index[,1]==1,]$width)
      L2 <- sum(tmp.new[index[,2]==1,]$width)
      W1 <- mean(index[tmp$region == 'SJ',1]*tmp[tmp$region == 'SJ']$ExpSum)#
      W2 <- mean(index[tmp$region == 'SJ',2]*tmp[tmp$region == 'SJ']$ExpSum)#
      type <- tmp.new$region
      type<-gsub('SJ','0',type)
      type<-gsub('ExonBin','1',type)
      type <- as.numeric(type)
      N.SJ <- tmp[tmp$region =='SJ',]$ExpSum
      
      dat.jags.SJ <- list(Index=index.SJ, nB=nrow(index.SJ),
                          nIso =2, N=N.SJ, 
                          W = sum(W1+W2),L=ln,
                          Le = sum(ln),
                          L1=L2,L2=L2)
      bayes.mod.fit.SJ <- jags(data = dat.jags.SJ,
                              inits = bayes.mod.inits,
                              parameters.to.save = bayes.mod.params,
                              n.iter = 3000,
                              n.chains = 1, 
                              n.burnin = 1000,
                              model.file = bayes.mod.SJ)
      res <- bayes.mod.fit.SJ$BUGSoutput$mean$theta
      posterior.iso1[[g]] <- data.frame(isoform=iso1,inferenced.psi=res,gene.id =g )
    }
  }
}
inferenced.psi.SJ <- do.call(rbind,posterior.iso1)


###use all region
bayes.mod<-function(){
    for (j in 1:nB) {
      p[j] <- ifelse(type[j] == 1, Index[j,1]*theta*L[j]/L1+Index[j,2]*(1-theta)*L[j]/L2,
                         Index[j,1]*theta + Index[j,2]*(1-theta))
      w[j] <- ifelse(type[j] == 1, W, sum(W1+W2))
      N[j] ~ dpois(p[j]*w[j])
     # N[j] ~ dbinom(p[j],W)#number of molecules we can observe

    }
                theta ~ dunif(0,1)#prior on theta

}
bayes.mod.params<-c("theta")
bayes.mod.inits<-function(){
        list(theta=runif(1,0,1))
}
posterior.iso1 <- list()
for(g in unique(all_exp$gene_id)){
  tmp <- all_exp[gene_id==g,]
  tmp <- tmp[tmp$ExpSum!=0,]
  setkey(tmp, start)
  if(nrow(tmp) == 0){
    posterior.iso1[[g]] <- data.frame(isoform=ifelse(nrow(tmp)==0,NA,unique(unlist(strsplit(tmp$tx_id,";|:")))[1]),
                                      inferenced.psi=NA,gene.id =g )
  }else{
  #only one isoform expressed, skip the model
    if(length(unique(unlist(strsplit(tmp$tx_id,";|:"))))  ==1 ){
      posterior.iso1[[g]] <- data.frame(isoform=unique(unlist(strsplit(tmp$tx_id,";|:"))),
                                        inferenced.psi=1,gene.id =g )
    }else{
      iso1 <- unique(unlist(strsplit(tmp$tx_id,";|:")))[1]
      iso2 <- unique(unlist(strsplit(tmp$tx_id,";|:")))[2]
      index <- matrix(rep(0,nrow(tmp)*2),c(nrow(tmp),2))
      index[grep(iso1, tmp$tx_id),1] <- 1 #index is 1 if isoform express this bin,otherwise 0
      index[grep(iso2, tmp$tx_id),2] <- 1
      tmp.new <- tmp
      tmp.new[tmp.new$region=='SJ',]$width <- 0
      ln <- tmp.new$width
      L1 <- sum(tmp.new[index[,1]==1,]$width)
      L2 <- sum(tmp.new[index[,2]==1,]$width)
      W1 <- mean(index[tmp$region == 'SJ',1]*tmp[tmp$region == 'SJ']$ExpSum)#
      W2 <- mean(index[tmp$region == 'SJ',2]*tmp[tmp$region == 'SJ']$ExpSum)#
      type <- tmp.new$region
      type<-gsub('SJ','0',type)
      type<-gsub('ExonBin','1',type)
      type <- as.numeric(type)
      #weight index matrix with sequcencing preference of each bin of each isoform
      if(nrow(index) ==1){
        weighted.index.exon <- index
      }else{
        seq1 <- tmp[index[,1]==1&region == 'ExonBin', ]$ExpSum
        seq2 <- tmp[index[,2]==1&region == 'ExonBin', ]$ExpSum
        #normalize the seq to 0-1
        seq1_norm <- normalize.seq(seq1+rnorm(length(seq1)))
        seq2_norm <- normalize.seq(seq2+rnorm(length(seq2)))
        LBC.weighted.index.exon <- index
        LBC.weighted.index.exon[type==1 & index[,1] ==1,1] <- seq1_norm
        LBC.weighted.index.exon[type==1 & index[,2]==1,2] <- seq2_norm
        
        #GBC weighted.index.exon
        if(unique(tmp$strand=='+')){
          L1.w <- tmp[index[,1]==1&region == 'ExonBin', ]$width
          L2.w <- tmp[index[,2]==1&region == 'ExonBin', ]$width
          L1.w.s <- Iter.sum(L1.w)
          qtiles.1 <- as.numeric(quantile(seq(1,sum(L1.w),1),probs = seq(0.1,1,0.1)))
          qtiles.list.1 <- sapply(L1.w.s, function(x)which.min(abs(qtiles.1-x)))
          GBC.weighted.1 <- sapply(qtiles.list.1, function(x)normalized.GBC[x*10])
          #
          L2.w.s <- Iter.sum(L2.w)
          qtiles.2 <- as.numeric(quantile(seq(1,sum(L2.w),1),probs = seq(0.1,1,0.1)))
          qtiles.list.2 <- sapply(L2.w.s, function(x)which.min(abs(qtiles.2-x)))
          GBC.weighted.2 <- sapply(qtiles.list.2, function(x)normalized.GBC[x*10])
          #
          GBC.weighted.index.exon <- index
          GBC.weighted.index.exon[type==1 & index[,1] ==1,1] <- GBC.weighted.1
          GBC.weighted.index.exon[type==1 & index[,2] ==1,2] <- GBC.weighted.2
        }else{
          L1.w <- rev(tmp[index[,1]==1&region == 'ExonBin', ]$width)
          L2.w <- rev(tmp[index[,2]==1&region == 'ExonBin', ]$width)
          L1.w.s <- Iter.sum(L1.w)
          qtiles.1 <- as.numeric(quantile(seq(1,sum(L1.w),1),probs = seq(0.1,1,0.1)))
          qtiles.list.1 <- sapply(L1.w.s, function(x)which.min(abs(qtiles.1-x)))
          GBC.weighted.1 <- sapply(qtiles.list.1, function(x)normalized.GBC[x*10])
          #
          L2.w.s <- Iter.sum(L2.w)
          qtiles.2 <- as.numeric(quantile(seq(1,sum(L2.w),1),probs = seq(0.1,1,0.1)))
          qtiles.list.2 <- sapply(L2.w.s, function(x)which.min(abs(qtiles.2-x)))
          GBC.weighted.2 <- sapply(qtiles.list.2, function(x)normalized.GBC[x*10])
          #
          GBC.weighted.index.exon <- index
          GBC.weighted.index.exon[type==1 & index[,1] ==1,1] <- rev(GBC.weighted.1)
          GBC.weighted.index.exon[type==1 & index[,2] ==1,2] <- rev(GBC.weighted.2)
        }
        
        #
        alpha <- 0.8
        weighted.index.exon <- GBC.weighted.index.exon*alpha + LBC.weighted.index.exon *(1-alpha)
      }
      dat.jags <- list(Index=weighted.index.exon, nB=nrow(tmp),
                       N=tmp.new$ExpSum, 
                       W = sum(tmp.new$ExpSum),
                       L=ln,type=type,
                       L1=L1,L2=L2,W1=W1,W2=W2)
    bayes.mod.fit <- jags(data = dat.jags,
                          inits = bayes.mod.inits,
                          parameters.to.save = bayes.mod.params,
                          n.iter = 3000,
                          n.chains = 1, 
                          n.burnin = 1000,
                          model.file = bayes.mod)
    res <- bayes.mod.fit$BUGSoutput$mean$theta
    posterior.iso1[[g]] <- data.frame(isoform=iso1,inferenced.psi=res,gene.id =g )
    }
  }
}
inferenced.psi.both <- do.call(rbind,posterior.iso1)










identical(inferenced.psi.SJ$gene.id,inferenced.psi.exon$gene.id)
identical(inferenced.psi.SJ$gene.id,inferenced.psi.both$gene.id)
inferenced.all <- data.frame(SJ = inferenced.psi.SJ$inferenced.psi,
                             Exon = inferenced.psi.exon$inferenced.psi,
                             both = inferenced.psi.both$inferenced.psi,
                             gene.id = inferenced.psi.both$gene.id,
                             isoform = inferenced.psi.both$isoform)
inferenced.all <- inferenced.all[!is.na(inferenced.all$isoform),]
saveRDS(inferenced.all, file = 'inferenced.all.rds')

library(corrplot)
test.data<- readRDS('inferenced.all.rds')
true.psi <- as.data.frame(depth_mat[tx_id%in%as.character(test.data$isoform),])
rownames(true.psi) <- true.psi$tx_id
true.psi <- true.psi[as.character(test.data$isoform), ]
identical(as.character(true.psi$tx_id),as.character(test.data$isoform))

plot(true.psi$PSI,test.data$both)

#calculate observed psi
observed.psi.SJ <- c()
for(g in test.data$gene.id){
  iso <- test.data[test.data$gene.id == g, ]$isoform
  tmp <- all_exp[all_exp$gene_id == g& region == 'SJ'&freq == 1,]
  observed.psi.SJ[g] <- mean(tmp[tmp$tx_id == iso,]$ExpSum)/(mean(tmp[tmp$tx_id == iso,]$ExpSum)+mean(tmp[tmp$tx_id != iso,]$ExpSum))
}

observed.psi.exon <- c()
for(g in test.data$gene.id){
  iso <- test.data[test.data$gene.id == g, ]$isoform
  tmp <- all_exp[all_exp$gene_id == g& region == 'ExonBin'&freq == 1,]
  observed.psi.exon[g] <- mean(tmp[tmp$tx_id == iso,]$ExpSum/tmp[tmp$tx_id == iso,]$width)/(mean(tmp[tmp$tx_id == iso,]$ExpSum/tmp[tmp$tx_id == iso,]$width)+mean(tmp[tmp$tx_id != iso,]$ExpSum/tmp[tmp$tx_id != iso,]$width))
}

all.data <- test.data
all.data$true.psi <- true.psi$PSI
all.data$observed.psi.SJ <- observed.psi.SJ
all.data$observed.psi.exon <- observed.psi.exon
saveRDS(all.data, file = 'res.all.data.rds')





colnames(all.data)[1:3] <- c("inferenced.SJ","inferenced.exon","inferenced.both")
all.data$observed.both.psi <- ifelse(!is.na(all.data$observed.psi.SJ),
                                 all.data$observed.psi.SJ, all.data$observed.psi.exon)


MSE <- function(x,y){
  df <- data.frame(x=x,y=y)
  df <- na.omit(df)
  mse <- round(mean((df$x-df$y)^2), 4)
  return(mse)
}

COR <- function(x,y){
  df <- data.frame(x=x,y=y)
  df <- na.omit(df)
  cor <- round(cor(df$x,df$y),4)
  return(cor)
}
Plot <- function(x,y,label){
  plot(x,y,main = paste0('MSE ', MSE(x,y),
                         '\n', 'Number of genes ', length(x) -sum(is.na(y)),
                         '\n', 'cor ', COR(x,y)),
       xlab = "True.psi", ylab = label)
}
pdf('test.all.3bias.pdf')
op <- par(mfrow = c(2, 3), pty = "s") 
Plot(all.data$true.psi, all.data$inferenced.both, "Inferenced.Both")
Plot(all.data$true.psi, all.data$inferenced.exon, "Inferenced.Exon")
Plot(all.data$true.psi, all.data$inferenced.SJ, "Inferenced.SJ")
Plot(all.data$true.psi, all.data$observed.both.psi, "Observed.Both")
Plot(all.data$true.psi, all.data$observed.psi.exon, "Observed.Exon")
Plot(all.data$true.psi, all.data$observed.psi.SJ, "Observed.SJ")
par(op)
dev.off()

pdf('test.SJ.3bias.pdf')
df <- all.data[!is.na(all.data$observed.psi.SJ),]
op <- par(mfrow = c(2, 3), pty = "s") 
Plot(df$true.psi, df$inferenced.both, "Inferenced.Both")
Plot(df$true.psi, df$inferenced.exon, "Inferenced.Exon")
Plot(df$true.psi, df$inferenced.SJ, "Inferenced.SJ")
Plot(df$true.psi, df$observed.both.psi, "Observed.Both")
Plot(df$true.psi, df$observed.psi.exon, "Observed.Exon")
Plot(df$true.psi, df$observed.psi.SJ, "Observed.SJ")
par(op)
dev.off()


