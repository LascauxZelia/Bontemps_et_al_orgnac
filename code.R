library(phyloseq)
library(minpack.lm)
library(Hmisc)
library(vegan)
library(stats)
library(ade4)
library(picante)
library(iCAMP)

# NNMDS

table = OTU[,5:ncol(OTU)]
m_com = as.matrix(table)
nmds = metaMDS(m_com,distance = "bray")
nmds

data.scores = as.data.frame(scores(nmds))
data.scores$Sample = OTU$Sample
data.scores$Month = OTU$Replicat
data.scores$Zones = OTU$Zone
data.scores$Inter = OTU$Interaction

#Neutral Community Model (NCM)
  require(minpack.lm)
  require(Hmisc)
  require(stats4)
  
  #Number of individuals per community
  Nu <- mean(apply(echt[,-1], 1, sum))
  
  #Calculate the average relative abundance of each taxa across communities
  if(is.null(meta)){
    p.m <- apply(echt[,-1], 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/Nu
  } else {
    p.m <- apply(meta[,-1], 2, mean)
    p.m <- p.m[p.m != 0]
    p <- p.m/Nu
  }
  
  #Calculate the occurrence frequency of each taxa
  echt.bi <- 1*(echt[,-1]>0)
  freq <- apply(echt.bi, 2, mean)
  freq <- freq[freq != 0]
  
  C <- merge(p, freq, by=0)
  C <- C[order(C[,2]),]
  C <- as.data.frame(C)
  C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
  freq <- C.0[,3]
  names(p) <- C.0[,1]
  names(freq) <- C.0[,1]
  d = 1/N
  
  ##Fit model parameter m
  m.fit <- nlsLM(freq ~ pbeta(d, Nu*m*p, Nu*m*(1-p), lower.tail=FALSE), start=list(m=0.1))
  m.ci <- confint(m.fit, 'm', level=0.95)
  
  ##Fit neutral model parameter m 
  sncm.LL <- function(m, sigma){
    R = freq - pbeta(d, Nu*m*p, N*m*(1-p), lower.tail=FALSE)
    R = dnorm(R, 0, sigma)
    -sum(log(R))
  }
  m.mle <- mle(sncm.LL, start=list(m=0.1, sigma=0.1), nobs=length(p))
  
  ##Akaike's Information Criterion (AIC)
  aic.fit <- AIC(m.mle, k=2)
  bic.fit <- BIC(m.mle)
  
  ##R-squared and Root Mean Squared Error
  freq.pred <- pbeta(d, Nu*coef(m.fit)*p, Nu*coef(m.fit)*(1-p), lower.tail=FALSE)
  Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE <- sqrt(sum((freq-freq.pred)^2)/(length(freq)-1))
  
  pred.ci <- binconf(freq.pred*nrow(echt[,-1]), nrow(echt[,-1]), alpha=0.05, method="wilson", return.df=TRUE)
  
  ##binomial model
  bino.pred <- pbinom(d, Nu, p, lower.tail=FALSE)
  Rsqr.bino <- 1 - (sum((freq - bino.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE.bino <- sqrt(sum((freq - bino.pred)^2)/(length(freq) - 1))
  
  bino.pred.ci <- binconf(bino.pred*nrow(echt[,-1]), nrow(echt[,-1]), alpha=0.05, method="wilson", return.df=TRUE)
  
    
  ##Poisson model
  pois.pred <- ppois(d, Nu*p, lower.tail=FALSE)
  Rsqr.pois <- 1 - (sum((freq - pois.pred)^2))/(sum((freq - mean(freq))^2))
  RMSE.pois <- sqrt(sum((freq - pois.pred)^2)/(length(freq) - 1))
  
  pois.pred.ci <- binconf(pois.pred*nrow(echt[,-1]), nrow(echt[,-1]), alpha=0.05, method="wilson", return.df=TRUE)
  
      ##AIC binomial model
  bino.LL <- function(mu, sigma){
    R = freq - pbinom(d, Nu, p, lower.tail=FALSE)
    R = dnorm(R, mu, sigma)
    -sum(log(R))
  }
  bino.mle <- mle(bino.LL, start=list(mu=0.1, sigma=0.1), nobs=length(p))
  
  aic.bino <- AIC(bino.mle, k=2)
  bic.bino <- BIC(bino.mle)
  
  ##AIC Poisson model
  pois.LL <- function(mu, sigma){
    R = freq - ppois(d, N*p, lower.tail=FALSE)
    R = dnorm(R, mu, sigma)
    -sum(log(R))
  }
  pois.mle <- mle(pois.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
  
  aic.pois <- AIC(pois.mle, k=2)
  bic.pois <- BIC(pois.mle)

  
  ##Results
    fitstats <- data.frame(m=numeric(),m.ci=numeric(), m.mle=numeric(), maxLL=numeric(), binoLL=numeric(), poisLL=numeric(), Rsqr=numeric(), Rsqr.bino=numeric(), Rsqr.pois=numeric(), RMSE=numeric(), RMSE.bino=numeric(), RMSE.pois=numeric(), AIC=numeric(), BIC=numeric(), AIC.bino=numeric(), BIC.bino=numeric(), AIC.pois=numeric(), BIC.pois=numeric(), N=numeric(), Samples=numeric(), Richness=numeric(), Detect=numeric())
    fitstats[1,] <- c(coef(m.fit),coef(m.fit)-m.ci[1], m.mle@coef['m'], m.mle@details$value, bino.mle@details$value, pois.mle@details$value, Rsqr, Rsqr.bino, Rsqr.pois, RMSE, RMSE.bino, RMSE.pois, aic.fit, bic.fit, aic.bino, bic.bino, aic.pois, bic.pois, Nu, nrow(echt), length(p), d)
    A <- cbind(p, freq, freq.pred, pred.ci[,2:3], bino.pred, bino.pred.ci[,2:3])
    A <- as.data.frame(A)
    colnames(A) <- c('p', 'freq', 'freq.pred', 'pred.lower', 'pred.uper', 'binomial.pred', 'binomial.lower', 'binomial.uper')

#BNTI (Stegen et al. 2013, https://doi.org/10.1038/ismej.2013.93).
require(picante)
    
setwd("~/")
bnti <- read_excel("~/")
setwd("~/")
phylo = read.tree("XX.txt");

match.phylo.otu = match.phylo.data(phylo, t(bnti));
str(match.phylo.otu);

beta.mntd.weighted = as.matrix(comdistnt(t(match.phylo.otu$data),cophenetic(match.phylo.otu$phy),abundance.weighted=T));
dim(beta.mntd.weighted);
beta.mntd.weighted[1:5,1:5];
beta.reps = 999;

rand.weighted.bMNTD.comp = array(c(-999),dim=c(ncol(match.phylo.otu$data),ncol(match.phylo.otu$data),beta.reps));
dim(rand.weighted.bMNTD.comp);

for (rep in 1:beta.reps) {
  
  rand.weighted.bMNTD.comp[,,rep] = as.matrix(comdistnt(t(match.phylo.otu$data),taxaShuffle(cophenetic(match.phylo.otu$phy)),abundance.weighted=T,exclude.conspecifics = F));
  
  print(c(date(),rep));
  
}

weighted.bNTI = matrix(c(NA),nrow=ncol(match.phylo.otu$data),ncol=ncol(match.phylo.otu$data));
dim(weighted.bNTI);

for (columns in 1:(ncol(match.phylo.otu$data)-1)) {
  for (rows in (columns+1):ncol(match.phylo.otu$data)) {
    
    rand.vals = rand.weighted.bMNTD.comp[rows,columns,];
    weighted.bNTI[rows,columns] = (beta.mntd.weighted[rows,columns] - mean(rand.vals)) / sd(rand.vals);
    rm("rand.vals");
    
  };
};

rownames(weighted.bNTI) = colnames(match.phylo.otu$data);
colnames(weighted.bNTI) = colnames(match.phylo.otu$data);
weighted.bNTI;
setwd("~/")
write.csv(weighted.bNTI,"XX_weighted_bNTI.csv",quote=F, row.names = FALSE)

#RCbray
require(iCAMP)
bnti <- read_excel("~/")

x <- RC.pc(bnti, rand = 1000, na.zero = TRUE, nworker = 4,
      memory.G = 50, weighted = TRUE, unit.sum = NULL,
      meta.ab = NULL,
      detail.null=FALSE,output.bray=FALSE,silent=FALSE)
write.table(x,"XX.xls",quote = F, row.names = FALSE)
