#compare RSEM and Bayes
#RSEM results
nobias_rsem <- read.table('/mnt/raid61/APA/Project/RSEM_SIMULATION/nobias.isoforms.results',
                           header = T)[,c(1,8)]
nobias.bayes.true <- readRDS('/mnt/data8/xiaoxia/bayes_bulk/Nobias/res.all.data.rds')

nobias <- merge(nobias_rsem, nobias.bayes.true, by.x = 'transcript_id',by.y = 'isoform',
                 all.x = F, all.y = T)
nobias$IsoPct <- nobias$IsoPct*0.01
bias_rsem <- read.table('/mnt/raid61/APA/Project/RSEM_SIMULATION/bias.isoforms.results',
                        header = T)[,c(1,8)]


bias.bayes.true <- readRDS('/mnt/data8/xiaoxia/bayes_bulk/3bias/res.all.data.rds')
bias <- merge(bias_rsem, bias.bayes.true, by.x = 'transcript_id',by.y = 'isoform',
                 all.x = F, all.y = T)
bias$IsoPct <- bias$IsoPct*0.01

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
Plot <- function(x,y,labelx,labely){
  plot(x,y,main = paste0('MSE ', MSE(x,y),
                         '\n', 'Number of genes ', length(x) -sum(is.na(y)),
                         '\n', 'cor ', COR(x,y)),
       xlab = labelx, ylab = labely)
}
pdf('biasCompare.pdf')
Plot(bias$true.psi, bias$IsoPct,"True.psi","RSEM")
Plot(bias$true.psi, bias$both,"True.psi","Bayes")
Plot(bias$IsoPct, bias$both, "RSEM","Bayes")

dev.off()

pdf('nobiasCompare.pdf')
Plot(nobias$true.psi, nobias$IsoPct,"True.psi","RSEM")
Plot(nobias$true.psi, nobias$both,"True.psi","Bayes")
Plot(nobias$IsoPct, nobias$both, "RSEM","Bayes")
dev.off()





