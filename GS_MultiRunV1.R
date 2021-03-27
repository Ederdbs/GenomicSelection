### Simulation Breeding Selection Strategy, Index of Selection and GS Models ---
## Packagers -------------------------------------------------------------------
rm(list=ls())
if(TRUE){
  #source('Founder_pop.R')
  load('Founders.RData')
}else{
  library('AlphaSimR')
  set.seed(123)
  Ne <- 106  # Tsuda 2015 BMC
  segSites  <-  round(20000/20,0)  # Genome sequence of the palaeopolyploid soybean
  founderPop <- quickHaplo(nInd=200,
                           nChr=20,
                           segSites=segSites,
                           genLen = 1.15,
                           ploidy = 2L,
                           inbred = TRUE)
  save('Founders.RData')
}
Xpackagers <- c('AlphaSimR','bWGR','parallel','foreach','doParallel',
                'reshape','ggplot2','gridExtra','lubridate','plyr',
                'ranger','Rcpp','keras','verification','rrBLUP')
XXX <- lapply(Xpackagers, function(x){suppressMessages(require(x,quietly = TRUE, character.only = TRUE))})
#Rcpp::sourceCpp('dnn_chol_norm.cpp')
# Simulation Parameters  -------------------------------------------------------
Number_of_runs = 2
Number_of_generations = 200
Intensity <- c(2.5,5,7.5,10)/100#c(0.005,0.01,0.025,0.05)
NF1 <-  300#200 ## NUmber of new population to be created 
NF2 <- 50#200 ### Number of individual by population in F2 and F_infinity
GS_Model <- list(
  #'XGB' = function(y,gen){require(xgboost); X = as(data.matrix(gen), "dgCMatrix"); fit0 = xgboost(data=X,label=y,params=list(subsample=0.25),nrounds=20,objective="reg:squarederror"); return(list(hat=predict(fit0,X)))},
  #'DNN' = function(y,gen){FNN(as.matrix(y),gen)},
  'RF' = function(y,gen,...){list(hat=ranger::ranger(y~.,data.frame(y=y,gen),verbose = FALSE, save.memory = TRUE,write.forest = FALSE,...)$predictions)},
  #'RKHS' = function(y,gen){ K = GAU(gen); diag(K)=diag(K)+0.00001; E = eigen(K,symmetric = T); fit = emML(y,E$vectors,E$values); return(fit)},
  #'BayesCpi'=BayesCpi,
  #'BayesDpi'=BayesDpi,
  'GBLUP'=emML,
  'RR'=emRR,
  'BayesA'=emBA,
  'BayesB'=emBB,
  'BayesC'=emBC,
  'BayesL'=emBL,
  'FLM' = emDE,
  'Randon' = function(y,gen){return(list(hat=sample(y)))},
  'Gv'= function(y,gen){return(list(hat=NA))},
  'Pheno' = function(y,gen){return(list(hat=NA))}
)
Sel_Strategy <- as.character(c('WBF', # Within the best family
                               'WIF', # Within family
                               'ACF'  # Across Family
))
# Genetics parameters ----------------------------------------------------------
POPBestIntensity <- 0.3
NCgs <- 3    
## Trait parameter -------------------------------------------------------------
#GxE_corr = 0.4912 ;#By_Loc_H2 = 0.0971 ; #Across_Loc_H2 = 0.5769;#h2= 0.12; #Acroos location 77; # GxE = 77; #Env = 120
mean = 60 # From SoyNAN 
var = 77  # From SoyNAN 
varGxE = 77 # From SoyNAN 
varEnv= 200 # From SoyNAN 
corA = matrix(1,nrow=1)
corGxE = matrix(1,nrow=1)
nReps=1
POSSI <- expand.grid(Sel_Strategy=as.character(Sel_Strategy),
                     run=1:Number_of_runs,
                     Intensity= Intensity,
                     Dmodel=1:length(GS_Model))
cat(paste('Number of Simulation: ',Number_of_runs*length(Intensity)*length(GS_Model)*length(Sel_Strategy),'\n'))
## Founder POP and Simulation Parameters ----------------------------------------
nSnpPerChr <-  round(6000/20,0) # Illumina 20K same soyNAN http://journals.atlas-publishing.org/index.php/PGGB/article/view/154
nQtlPerChr <- segSites * 0.7
SP <- SimParam$new(founderPop)
SP$addSnpChip(nSnpPerChr=nSnpPerChr)
SP$addTraitAEG(nQtlPerChr=nQtlPerChr, 
               mean = mean, 
               var = var, 
               varEnv=0,
               relAA = 0.5,
               varGxE = varGxE,
               corA =corA,
               corGxE = corGxE)
SP$setVarE(varE=varEnv) ##0.281 from soyNAN ##0.281 from soyNAN
# function ---------------------------------------------------------------
Mysimu <- function(j,...){
  run <- POSSI[j,'run']
  Dmodel <- POSSI[j,'Dmodel']
  Int <- POSSI[j,'Intensity']
  Strategy <- as.character(POSSI[j,'Sel_Strategy'])
  #cat(paste('Simu: ',sprintf("%04i",j),'time: ',Sys.time(),'\n'))
  pop = newPop(founderPop, simParam=SP)
  genMean = c();genVar = c();H2 = c();Accuracy = c();nIndSel = c();
  AccuracyF = c();AccuracyGvPhe = c();AccuracyPheEbv = c();CRPS=c();
  geno_MC <- data.frame(NULL);
  pheno_MC <- data.frame(NULL);
  Npred <- c();
  useM <- ifelse(names(GS_Model)[Dmodel] == 'Pheno','pheno',
                 ifelse(names(GS_Model)[Dmodel] == 'Gv','gv',
                        'ebv'))
  pop = randCross(pop, nCrosses=20, nProgeny = NF1, simParam=SP)
  pop = makeDH(pop, nDH=1, simParam=SP)
  for(i in 1:Number_of_generations){
    Time_started <- ymd_hms(Sys.time())
    if(i==1){
      NIndSel <- pop@nInd
      pop = selectCross(pop, nInd=pop@nInd, nCrosses=NF1, use="pheno",nProgeny=1, simParam=SP)
    }else{
      if(Strategy == 'WBF'){ #WithinBestFamily POPBestIntensity (0.3)
        pop = selectFam(pop,nFam=round(NF1*POPBestIntensity),use=useM,simParam=SP)
        pop = selectWithinFam(pop,nInd=round((NF2*Int)/POPBestIntensity),use=useM,simParam=SP)
        nIndSel = pop@nInd
        pop = selectCross(pop, nInd=pop@nInd, nCrosses=NF1,use=useM, nProgeny=1, simParam=SP)
      }
      if(Strategy == 'WIF'){ #WithinFamily
        pop = selectWithinFam(pop,nInd=round(NF2*Int),use=useM,simParam=SP)
        nIndSel = pop@nInd
        pop = selectCross(pop, nInd=pop@nInd, nCrosses=NF1,use=useM, nProgeny=1, simParam=SP)
      } 
      if(Strategy == 'ACF'){ #AcrossFamily
        nIndSel = round(NF2*NF1*Int)
        pop = selectCross(pop,nInd=round(NF2*NF1*Int),use=useM,nCrosses=NF1,nProgeny=1,simParam=SP)  #F1
      }
    }
    pop = self(pop, nProgeny=NF2, simParam=SP) #F2
    pop = self(pop, nProgeny=1, simParam=SP)
    pop = self(pop, nProgeny=1, simParam=SP)
    gen = pullSnpGeno(pop, simParam = SP) ## Genotype
    ### multiple cycles prediction 
    geno_MC <- rbind(data.frame(i=i,gen),geno_MC)
    pheno_MC <- rbind(data.frame(i=i,pop@pheno[,1]-mean(pop@pheno[,1])),pheno_MC)
    cycles <- i:(i-NCgs+1)
    cycles <- cycles[cycles>0]
    geno_MC <- geno_MC[geno_MC$i %in% cycles,]
    pheno_MC <- pheno_MC[pheno_MC$i %in% cycles,]
    if(length(unique(gen%*%rnorm(ncol(gen))))<=5){
      fit = list(hat=rnorm(length(fit$hat)))
    }else{
      fit = GS_Model[[Dmodel]](as.matrix(pheno_MC[,-1]),as.matrix(geno_MC[,-1]))
      if(anyNA(fit$hat)){
        fit$hat=rnorm(length(fit$hat))
      } 
    }
    pop@ebv <- as.matrix(fit$hat[pheno_MC$i == i],ncol=1)
    genMean = c(genMean, meanG(pop))
    genVar = c(genVar, varG(pop))
    H2 = c(H2,varG(pop)/varP(pop))
    Accuracy = c(Accuracy,cor(pop@gv,pop@ebv))
    AccuracyGvPhe = c(AccuracyGvPhe,cor(pop@gv,pop@pheno))
    AccuracyPheEbv = c(AccuracyPheEbv,cor(pop@pheno,pop@ebv))
    CRPS = c(CRPS,crps(as.numeric(pop@ebv),c(mean(pop@gv),sd(pop@gv)))$CRPS)
    NIndSel = c(NIndSel, nIndSel)
    Npred = c(Npred,nrow(pheno_MC))
    # Ppairs <- GGally::ggpairs(data.frame(Pheno=pop@pheno,gv=pop@gv,ebv=pop@ebv),#mapping=ggplot2::aes(colour = 'blue'),
    #                           upper = list(continuous = GGally::wrap("cor", method= "spearman")),
    #                           lower=list(continuous=function(data, mapping, ...){ 
    #                             ggplot(data = data, mapping = mapping, ...) + 
    #                               geom_point(...,size = 1) +
    #                               geom_abline(...)+
    #   scale_y_continuous(limits = c(50,150)) +
    #   scale_x_continuous(limits = c(50,150))}))+
    #   labs(title=paste(names(GS_Model)[Dmodel], as.character(Strategy), 'Generation: ',i))
    # print(Ppairs)
    cat(paste('Simu: ',sprintf("%04i",j),'Generation: ',sprintf("%03i",i),'time: ',Sys.time(),'Processing: ',
              sprintf("%03i",round(as.numeric(as.duration(interval(Time_started,ymd_hms(Sys.time()))),"minutes"))
              ),'min model: ',names(GS_Model)[Dmodel],'\n')) ## Remove
    
  }
  # Store run 
  RES <- list()
  RES[[paste0(names(GS_Model)[Dmodel],'_',Strategy,'_',Int,'_',run)]] <- t(unlist(list(GS_Model=names(GS_Model)[Dmodel],
                                                                                       Strategy=as.character(Strategy),
                                                                                       Intensity=Int,
                                                                                       Mu=genMean,
                                                                                       GV=genVar,
                                                                                       He=H2,
                                                                                       Accuracy = Accuracy,
                                                                                       AccuracyGvPhe = AccuracyGvPhe,
                                                                                       AccuracyPheEbv = AccuracyPheEbv,
                                                                                       CRPS = CRPS,
                                                                                       IndSel=NIndSel,
                                                                                       Npred = Npred)))
  return(RES)
  rm(pop)
}
# Replicate loop ---------------------------------------------------------------
sysinf <- Sys.info()
os <- sysinf['sysname']
if(os == 'Windows'){
  lP <- list()
  pdf('Coorelation.pdf')
  for(j in 1:nrow(POSSI)){
    lP[[j]] <-  try(Mysimu(j),TRUE)
  }
  dev.off()
}else{
  registerDoSEQ()
  hosts <- as.vector(unique(unlist(strsplit(as.character(Sys.getenv("LSB_HOSTS"))," "))))
  nh <- length(hosts);nc <- length(unlist(strsplit(as.character(Sys.getenv("LSB_HOSTS"))," ")))
  cl <- parallel::makePSOCKcluster(names=rep(hosts , each = floor(nc/nh)),
                                   outfile = "debug.txt",
                                   master=nsl(Sys.info()['nodename']),
                                   revtunnel = TRUE,
                                   outfile='',
                                   useXDR = TRUE)
  #clusterSetRNGStream(cl,7777777)
  doParallel::registerDoParallel(cl=cl)
  lP <- foreach(j = 1:nrow(POSSI),
                .packages = Xpackagers,
                .verbose=FALSE,
                #.export = ls(globalenv()),
                .inorder=FALSE) %dopar% {
                  tryCatch(Mysimu(j),error=function(e){return('try-error')})
                }
  parallel::stopCluster(cl)
}
RES <- plyr::ldply(lapply(lP[sapply(lP,function(x){return(!(x=='try-error'))})],plyr::ldply))
JOBID <- Sys.getenv("LSB_JOBID")
write.csv(RES,paste0(JOBID,'_o_Results.csv'),row.names = F)
save.image(file=paste0(JOBID,'_All.RData'))
#RES <- read.csv('203825_o_Results.csv')
#load('203825_All.RData')
### Get results and print ------------------------------------------------------
RES1 <- melt(RES[,-1],id=1:3)
RES1$par <- gsub('[0-9]','',RES1$variable)
RES1$sim <- gsub('[A-z]','',RES1$variable)
RES1 <- transform(RES1,
                  value=as.numeric(as.character(value)),
                  sim=as.numeric(as.numeric(sim)),
                  par=as.character(par))

### compare Models--------------------------------------------------------------
plots <- list()
for(i in unique(RES1$par)){
  plots[[i]] <- ggplot(RES1[RES1$par %in% i,],aes(x=sim, y=value,color=GS_Model)) +
    geom_point(cex=0.5)+
    geom_smooth(method = "loess",lwd=1.5,se=FALSE) +
    xlab('Generation')+
    ylab(i)+
    facet_grid(Strategy~Intensity)+
    ggtitle(paste("NF1:",NF1,"|NF2:",NF2,"|segSites:",segSites,"|nQtlPerChr:" ,nQtlPerChr, '|nSnpPerChr:',nSnpPerChr))
}
pdf(paste0(JOBID,'_ResultsModels.pdf'),w=20,h=65)
eval(parse(text=eval(paste('grid.arrange(',paste0(paste0('plots[["',names(plots),'"]]'),collapse=','),',nrow=',length(names(plots)),')'))))
dev.off()

### compare Strategy -----------------------------------------------------------
plots <- list()
for(i in unique(RES1$par)){
  plots[[i]] <- ggplot(RES1[RES1$par %in% i,],aes(x=sim, y=value,color=Strategy)) +
    geom_point(cex=0.5)+
    geom_smooth(method = "loess",lwd=1.5,se=FALSE) +
    xlab('Generation')+
    ylab(i)+
    facet_grid(GS_Model ~ Intensity)+
    ggtitle(paste("NF1:",NF1,"|NF2:",NF2,"|segSites:",segSites,"|nQtlPerChr:" ,nQtlPerChr, '|nSnpPerChr:',nSnpPerChr))
}
pdf(paste0(JOBID,'_ResultsStrategy.pdf'),w=20,h=65)
eval(parse(text=eval(paste('grid.arrange(',paste0(paste0('plots[["',names(plots),'"]]'),collapse=','),',nrow=',length(names(plots)),')'))))
dev.off()

### compare Intensity-----------------------------------------------------------
plots <- list()
for(i in unique(RES1$par)){
  plots[[i]] <- ggplot(RES1[RES1$par %in% i,],aes(x=sim, y=value,color=as.factor(Intensity))) +
    geom_point(cex=0.5)+
    geom_smooth(method = "loess",lwd=1.5,se=FALSE) +
    xlab('Generation')+
    ylab(i)+
    facet_grid( Strategy ~ GS_Model)+
    ggtitle(paste("NF1:",NF1,"|NF2:",NF2,"|segSites:",segSites,"|nQtlPerChr:" ,nQtlPerChr, '|nSnpPerChr:',nSnpPerChr))
}
pdf(paste0(JOBID,'_ResultsIntensity.pdf'),w=20,h=65)
eval(parse(text=eval(paste('grid.arrange(',paste0(paste0('plots[["',names(plots),'"]]'),collapse=','),',nrow=',length(names(plots)),')'))))
dev.off()

### Analysis of genetic mean --------------------------------------------------
RES_sel <- RES1[RES1$par %in% c('Mu'),]
sModels <- unique(RES_sel$GS_Model)
sModels <- sModels[sModels != 'Gv']
resoutC <- data.frame(NULL)
for(iGS_Model in sModels){
  for(iStrategy in unique(RES_sel$Strategy)){
    for(iIntensity in unique(RES_sel$Intensity)){
      bsel <- dplyr::filter(RES_sel,
                            GS_Model %in% c(iGS_Model,'Gv') &
                              Strategy == iStrategy &
                              Intensity == iIntensity)
      spf_c <- splinefun(x=bsel[bsel$GS_Model == iGS_Model,'sim'],y=bsel[bsel$GS_Model == iGS_Model,'value'])
      spf_b <- splinefun(x=bsel[bsel$GS_Model == 'Gv','sim'],y=bsel[bsel$GS_Model == 'Gv','value'])
      Ncycles = max(bsel$sim)
      Int_b <- Int_c <- c()
      for(i in 1:Ncycles){
        Int_b <- c(Int_b,integrate(spf_b,1,i)$value)
        Int_c <- c(Int_c,integrate(spf_c,1,i)$value)
      }
      pred_b = spf_b(1:Ncycles)
      pred_c = spf_c(1:Ncycles)
      diff = pred_b-pred_c
      rela = pred_c/pred_b
      gg = c(c(pred_c,pred_c[Ncycles]) - c(pred_c[1],pred_c))[1:Ncycles]
      gg_per = c(c(c(pred_c,pred_c[Ncycles]) - c(pred_c[1],pred_c))[1:Ncycles]/pred_c)*100
      accum = Int_c/Int_b
      accum[1] <- 1
      resout <- data.frame(
        pred_b = pred_b,
        pred_c = pred_c,
        diff = diff,
        rela = rela,
        gg = gg,
        gg_per = gg_per,
        Int_b=Int_b,
        Int_c=Int_c,
        accum = accum)
      resout$iGS_Model <- iGS_Model
      resout$iStrategy <- iStrategy
      resout$iIntensity <- iIntensity
      resout$sim <- 1:Ncycles
      resoutC <- rbind(resoutC,resout)
    }
  }
}
write.csv(resoutC,paste0(JOBID,'_o_Results_mean.csv'),row.names = F)
RES_x <- melt(resoutC,id=c('iGS_Model','iStrategy','iIntensity','sim'))
RES_x <- transform(RES_x,
                   iGS_Model=as.factor(iGS_Model),
                   iStrategy=as.factor(iStrategy))
write.csv(RES_x,paste0(JOBID,'_o_Results_meanMelt.csv'),row.names = F)
### by mean GS MOdels-----------------------------------------------------------
plots <- list()
for(i in unique(RES_x$variable)){
  plots[[i]] <- ggplot(RES_x[RES_x$variable %in% i,],aes(x=sim, y=value,color=iGS_Model)) +
    geom_smooth(method = "loess",lwd=1.5,se=FALSE) +
    xlab('Generation')+
    ylab(i)+
    facet_grid(iStrategy~iIntensity)+
    ggtitle(paste("NF1:",NF1,"|NF2:",NF2,"|segSites:",segSites,"|nQtlPerChr:" ,nQtlPerChr, '|nSnpPerChr:',nSnpPerChr))
}
pdf(paste0(JOBID,'_ResultsGmeanGSModel.pdf'),w=20,h=65)
eval(parse(text=eval(paste('grid.arrange(',paste0(paste0('plots[["',names(plots),'"]]'),collapse=','),',nrow=',length(names(plots)),')'))))
dev.off()
### by mean Strategy -----------------------------------------------------------
plots <- list()
for(i in unique(RES_x$variable)){
  plots[[i]] <- ggplot(RES_x[RES_x$variable %in% i,],aes(x=sim, y=value,color=iStrategy)) +
    geom_smooth(method = "loess",lwd=1.5,se=FALSE) +
    xlab('Generation')+
    ylab(i)+
    facet_grid(iIntensity ~ iGS_Model)+
    ggtitle(paste("NF1:",NF1,"|NF2:",NF2,"|segSites:",segSites,"|nQtlPerChr:" ,nQtlPerChr, '|nSnpPerChr:',nSnpPerChr))
}
pdf(paste0(JOBID,'_ResultsGmeanStrategy.pdf'),w=20,h=65)
eval(parse(text=eval(paste('grid.arrange(',paste0(paste0('plots[["',names(plots),'"]]'),collapse=','),',nrow=',length(names(plots)),')'))))
dev.off()
### by mean GS MOdels-----------------------------------------------------------
plots <- list()
for(i in unique(RES_x$variable)){
  plots[[i]] <- ggplot(RES_x[RES_x$variable %in% i,],aes(x=sim, y=value,color=as.factor(iIntensity))) +
    geom_smooth(method = "loess",lwd=1.5,se=FALSE) +
    xlab('Generation')+
    ylab(i)+
    facet_grid(iStrategy ~ iGS_Model)+
    ggtitle(paste("NF1:",NF1,"|NF2:",NF2,"|segSites:",segSites,"|nQtlPerChr:" ,nQtlPerChr, '|nSnpPerChr:',nSnpPerChr))
}
pdf(paste0(JOBID,'_ResultsGmeanIntensity.pdf'),w=20,h=65)
eval(parse(text=eval(paste('grid.arrange(',paste0(paste0('plots[["',names(plots),'"]]'),collapse=','),',nrow=',length(names(plots)),')'))))
dev.off()

##### classic all parameters ---------------------------------------------------
pMu <- ggplot(RES1[RES1$par %in% c('Mu'),],aes(x=sim, y=value,color=GS_Model)) +
  geom_smooth(method = "loess",lwd=1.5) +
  xlab('Generation')+ylab(expression(mu))+
  facet_grid(Strategy~Intensity)+
  ggtitle(paste("NF1:",NF1,"|NF2:",NF2,"|segSites:",segSites,"|nQtlPerChr:" ,nQtlPerChr, '|nSnpPerChr:',nSnpPerChr))

pGv <- ggplot(RES1[RES1$par %in% c('GV'),],aes(x=sim, y=value,color=GS_Model)) +
  geom_smooth(method = "loess",lwd=1.5) +
  xlab('Generation')+ylab(expression(sigma[g]^2))+
  facet_grid(Strategy~Intensity)+
  ggtitle(paste("NF1:",NF1,"|NF2:",NF2,"|segSites:",segSites,"|nQtlPerChr:" ,nQtlPerChr, '|nSnpPerChr:',nSnpPerChr)) 

pHe <- ggplot(RES1[RES1$par %in% c('He'),],aes(x=sim, y=value,color=GS_Model)) +
  geom_smooth(method = "loess",lwd=1.5) +
  xlab('Generation')+ylab(expression(H^2))+
  facet_grid(Strategy~Intensity)+
  ggtitle(paste("NF1:",NF1,"|NF2:",NF2,"|segSites:",segSites,"|nQtlPerChr:" ,nQtlPerChr, '|nSnpPerChr:',nSnpPerChr))

pAc <- ggplot(RES1[RES1$par %in% c('Accuracy'),],aes(x=sim, y=value,color=GS_Model)) +
  geom_smooth(method = "loess",lwd=1.5) +
  xlab('Generation')+ylab('Accuracy')+
  facet_grid(Strategy~Intensity)+
  ggtitle(paste("NF1:",NF1,"|NF2:",NF2,"|segSites:",segSites,"|nQtlPerChr:" ,nQtlPerChr, '|nSnpPerChr:',nSnpPerChr))

pAcGvPhe <- ggplot(RES1[RES1$par %in% c('AccuracyGvPhe'),],aes(x=sim, y=value,color=GS_Model)) +
  geom_smooth(method = "loess",lwd=1.5) +
  xlab('Generation')+ylab('Accuracy gv x  pheno')+
  facet_grid(Strategy~Intensity)+
  ggtitle(paste("NF1:",NF1,"|NF2:",NF2,"|segSites:",segSites,"|nQtlPerChr:" ,nQtlPerChr, '|nSnpPerChr:',nSnpPerChr))

pAcPheEbv <- ggplot(RES1[RES1$par %in% c('AccuracyPheEbv'),],aes(x=sim, y=value,color=GS_Model)) +
  geom_smooth(method = "loess",lwd=1.5) +
  xlab('Generation')+ylab('Accuracy pheno x ebv')+
  facet_grid(Strategy~Intensity)+
  ggtitle(paste("NF1:",NF1,"|NF2:",NF2,"|segSites:",segSites,"|nQtlPerChr:" ,nQtlPerChr, '|nSnpPerChr:',nSnpPerChr))

pCRPS <- ggplot(RES1[RES1$par %in% c('CRPS'),],aes(x=sim, y=value,color=GS_Model)) +
  geom_smooth(method = "loess",lwd=1.5) +
  xlab('Generation')+ylab('CRPS')+
  facet_grid(Strategy~Intensity)+
  ggtitle(paste("NF1:",NF1,"|NF2:",NF2,"|segSites:",segSites,"|nQtlPerChr:" ,nQtlPerChr, '|nSnpPerChr:',nSnpPerChr))

pNp <- ggplot(RES1[RES1$par %in% c('Npred'),],aes(x=sim, y=value,color=GS_Model)) +
  geom_point()+
  xlab('Generation')+ylab('Npred')+
  facet_grid(Strategy~Intensity)+
  ggtitle(paste("NF1:",NF1,"|NF2:",NF2,"|segSites:",segSites,"|nQtlPerChr:" ,nQtlPerChr, '|nSnpPerChr:',nSnpPerChr))

pNs <- ggplot(RES1[RES1$par %in% c('IndSel'),],aes(x=sim, y=value,color=GS_Model)) +
  geom_point()+
  xlab('Generation')+ylab('N')+
  facet_grid(Strategy~Intensity)+
  ggtitle(paste("NF1:",NF1,"|NF2:",NF2,"|segSites:",segSites,"|nQtlPerChr:" ,nQtlPerChr, '|nSnpPerChr:',nSnpPerChr))

pdf(paste0(JOBID,'_ResultsAllLines.pdf'),w=20,h=65)
grid.arrange(pMu,pGv,pHe,pAc,pAcGvPhe,pAcPheEbv,pCRPS,pNp,pNs,nrow = 9)
dev.off()

### Points ---------------------------
pMu <- ggplot(RES1[RES1$par %in% c('Mu'),],aes(x=sim, y=value,color=GS_Model)) +
  geom_point(cex=0.5)+
  geom_smooth(method = "loess",lwd=1.5) +
  xlab('Generation')+ylab(expression(mu))+
  facet_grid(Strategy~Intensity)+
  ggtitle(paste("NF1:",NF1,"|NF2:",NF2,"|segSites:",segSites,"|nQtlPerChr:" ,nQtlPerChr, '|nSnpPerChr:',nSnpPerChr))

pGv <- ggplot(RES1[RES1$par %in% c('GV'),],aes(x=sim, y=value,color=GS_Model)) +
  geom_point(cex=0.5)+
  geom_smooth(method = "loess",lwd=1.5) +
  xlab('Generation')+ylab(expression(sigma[g]^2))+
  facet_grid(Strategy~Intensity)+
  ggtitle(paste("NF1:",NF1,"|NF2:",NF2,"|segSites:",segSites,"|nQtlPerChr:" ,nQtlPerChr, '|nSnpPerChr:',nSnpPerChr)) 

pHe <- ggplot(RES1[RES1$par %in% c('He'),],aes(x=sim, y=value,color=GS_Model)) +
  geom_point(cex=0.5)+
  geom_smooth(method = "loess",lwd=1.5) +
  xlab('Generation')+ylab(expression(H^2))+
  facet_grid(Strategy~Intensity)+
  ggtitle(paste("NF1:",NF1,"|NF2:",NF2,"|segSites:",segSites,"|nQtlPerChr:" ,nQtlPerChr, '|nSnpPerChr:',nSnpPerChr))

pAc <- ggplot(RES1[RES1$par %in% c('Accuracy'),],aes(x=sim, y=value,color=GS_Model)) +
  geom_point(cex=0.5)+
  geom_smooth(method = "loess",lwd=1.5) +
  xlab('Generation')+ylab('Accuracy')+
  facet_grid(Strategy~Intensity)+
  ggtitle(paste("NF1:",NF1,"|NF2:",NF2,"|segSites:",segSites,"|nQtlPerChr:" ,nQtlPerChr, '|nSnpPerChr:',nSnpPerChr))

pAcF <- ggplot(RES1[RES1$par %in% c('AccuracyF'),],aes(x=sim, y=value,color=GS_Model)) +
  geom_point(cex=0.5)+
  geom_smooth(method = "loess",lwd=1.5) +
  xlab('Generation')+ylab('Accuracy Pheno')+
  facet_grid(Strategy~Intensity)+
  ggtitle(paste("NF1:",NF1,"|NF2:",NF2,"|segSites:",segSites,"|nQtlPerChr:" ,nQtlPerChr, '|nSnpPerChr:',nSnpPerChr))

pAcGvPhe <- ggplot(RES1[RES1$par %in% c('AccuracyGvPhe'),],aes(x=sim, y=value,color=GS_Model)) +
  geom_point(cex=0.5)+
  geom_smooth(method = "loess",lwd=1.5) +
  xlab('Generation')+ylab('Accuracy gv x  pheno')+
  facet_grid(Strategy~Intensity)+
  ggtitle(paste("NF1:",NF1,"|NF2:",NF2,"|segSites:",segSites,"|nQtlPerChr:" ,nQtlPerChr, '|nSnpPerChr:',nSnpPerChr))

pAcPheEbv <- ggplot(RES1[RES1$par %in% c('AccuracyPheEbv'),],aes(x=sim, y=value,color=GS_Model)) +
  geom_point(cex=0.5)+
  geom_smooth(method = "loess",lwd=1.5) +
  xlab('Generation')+ylab('Accuracy pheno x ebv')+
  facet_grid(Strategy~Intensity)+
  ggtitle(paste("NF1:",NF1,"|NF2:",NF2,"|segSites:",segSites,"|nQtlPerChr:" ,nQtlPerChr, '|nSnpPerChr:',nSnpPerChr))

pCRPS <- ggplot(RES1[RES1$par %in% c('CRPS'),],aes(x=sim, y=value,color=GS_Model)) +
  geom_point(cex=0.5)+
  geom_smooth(method = "loess",lwd=1.5) +
  xlab('Generation')+ylab('CRPS')+
  facet_grid(Strategy~Intensity)+
  ggtitle(paste("NF1:",NF1,"|NF2:",NF2,"|segSites:",segSites,"|nQtlPerChr:" ,nQtlPerChr, '|nSnpPerChr:',nSnpPerChr))

pNp <- ggplot(RES1[RES1$par %in% c('Npred'),],aes(x=sim, y=value,color=GS_Model)) +
  geom_point()+
  xlab('Generation')+ylab('Npred')+
  facet_grid(Strategy~Intensity)+
  ggtitle(paste("NF1:",NF1,"|NF2:",NF2,"|segSites:",segSites,"|nQtlPerChr:" ,nQtlPerChr, '|nSnpPerChr:',nSnpPerChr))

pNs <- ggplot(RES1[RES1$par %in% c('IndSel'),],aes(x=sim, y=value,color=GS_Model)) +
  geom_point()+
  xlab('Generation')+ylab('N')+
  facet_grid(Strategy~Intensity)+
  ggtitle(paste("NF1:",NF1,"|NF2:",NF2,"|segSites:",segSites,"|nQtlPerChr:" ,nQtlPerChr, '|nSnpPerChr:',nSnpPerChr))

pdf(paste0(JOBID,'_ResultsAllPoints.pdf'),w=20,h=65)
grid.arrange(pMu,pGv,pHe,pAc,pAcGvPhe,pAcPheEbv,pCRPS,pNp,pNs,nrow = 9)
dev.off()

