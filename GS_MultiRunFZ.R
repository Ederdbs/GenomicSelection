### Simulation Breeding Selection Strategy, Index of Selection and GS Models ---
## Packagers -------------------------------------------------------------------
rm(list=ls())
if(TRUE){
  load('Founders.RData')
}else{
  library('AlphaSimR')
  set.seed(123)
  Ne <- 106  
  segSites  <-  round(20000/20,0)  # Genome sequence of the palaeopolyploid soybean
  founderPop <- quickHaplo(nInd=200,
                           nChr=20,
                           segSites=segSites,
                           genLen = 1.15,
                           ploidy = 2L,
                           inbred = TRUE)
}
Xpackagers <- c('AlphaSimR','bWGR','parallel','foreach','doParallel',
                'reshape','ggplot2','gridExtra','lubridate','plyr',
                'ranger','Rcpp','keras','verification','rrBLUP')
XXX <- lapply(Xpackagers, function(x){suppressMessages(require(x,quietly = TRUE, character.only = TRUE))})
# Simulation Parameters  -------------------------------------------------------
Number_of_runs = 1
Number_of_generations = 100
Intensity <- c(2.5,5,7.5,10)/100
FS <- list('1'=c(300,50), #F1 and F2 size
           '2'=c(250,60), #F1 and F2 size
           '3'=c(200,75), #F1 and F2 size
           '4'=c(150,100), #F1 and F2 size
           '5'=c(100,150)) #F1 and F2 size

GS_Model <- list(
  'GBLUP'=emML,
  'BayesA'=emBA,
  'BayesB'=emBB,
  'FLM' = emDE,
  'Gv'= function(y,gen){return(list(hat=NA))},
  'Pheno' = function(y,gen){return(list(hat=NA))}
)
Sel_Strategy <- as.character(c('WBF', # Within the best family
                               'WIF', # Within family
                               'ACF'  # Across Family
))
POSSI <- expand.grid(Sel_Strategy=as.character(Sel_Strategy),
                     run=1:Number_of_runs,
                     Intensity= Intensity,
                     Dmodel=1:length(GS_Model),
                     FS=names(FS))
FS <- ldply(FS)
colnames(FS) <- c('FS','F1','F2')
POSSI <- merge(POSSI,FS,all.x=TRUE)
cat(paste('Number of Simulation: ',nrow(POSSI),'\n'))
# Genetics parameters ----------------------------------------------------------
POPBestIntensity <- 0.3
NCgs <- 3    
## Trait parameter -------------------------------------------------------------
mean = 60 # From SoyNAN 
var = 77  # From SoyNAN 
varGxE = 77 # From SoyNAN 
varEnv= 200 # From SoyNAN 
corA = matrix(1,nrow=1)
corGxE = matrix(1,nrow=1)
nReps=1
## Founder POP and Simulation Parameters ----------------------------------------
nSnpPerChr <-  round(6000/20,0) 
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
  NF1 <- POSSI[j,'F1']
  NF2 <- POSSI[j,'F2']
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
    cat(paste('Simu: ',sprintf("%04i",j),'Generation: ',sprintf("%03i",i),'time: ',Sys.time(),'Processing: ',
              sprintf("%03i",round(as.numeric(as.duration(interval(Time_started,ymd_hms(Sys.time()))),"minutes"))
              ),'min model: ',names(GS_Model)[Dmodel],'\n')) ## Remove
    
  }
  # Store run 
  RES <- list()
  RES[[paste0(names(GS_Model)[Dmodel],'_',Strategy,'_',Int,'_',run)]] <- t(unlist(list(GS_Model=names(GS_Model)[Dmodel],
                                                                                       Strategy=as.character(Strategy),
                                                                                       Intensity=Int,
                                                                                       NF1=NF1,
                                                                                       NF2=NF2,
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
  cl <- makeCluster(6,revtunnel = TRUE)
  JOBID <- 0
}else{
  registerDoSEQ()
  hosts <- as.vector(unique(unlist(strsplit(as.character(Sys.getenv("LSB_HOSTS"))," "))))
  nh <- length(hosts)
  nc <- length(unlist(strsplit(as.character(Sys.getenv("LSB_HOSTS"))," ")))-1
  
  
  cl <- parallel::makePSOCKcluster(names=rep(hosts , each = floor(nc/nh)),
                                   outfile = "debug.txt",
                                   master=nsl(Sys.info()['nodename']),
                                   revtunnel = TRUE,
                                   outfile='',
                                   useXDR = TRUE)
  JOBID <- Sys.getenv("LSB_JOBID")
}
print(  rep(hosts , each = floor(nc/nh)))
print(hosts)
print(nh)
print(nc)


doParallel::registerDoParallel(cl=cl)
lP <- foreach(j = 1:nrow(POSSI),
              .packages = Xpackagers,
              .verbose=FALSE,
              .inorder=FALSE) %dopar% {
                tryCatch(Mysimu(j),error=function(e){return('try-error')})
              }
doParallel::stopImplicitCluster()

RES <- plyr::ldply(lapply(lP[sapply(lP,function(x){return(!(x=='try-error'))})],plyr::ldply))
write.csv(RES,paste0(JOBID,'_o_Results.csv'),row.names = F)
save.image(file=paste0(JOBID,'_All.RData'))
#load('245180_All.RData')
### Get results and print ------------------------------------------------------
RES1 <- melt(RES[,-1],id=1:5)
RES1$par <- gsub('[0-9]','',RES1$variable)
RES1$sim <- gsub('[A-z]','',RES1$variable)
RES1 <- transform(RES1,
                  value=as.numeric(as.character(value)),
                  sim=as.numeric(as.numeric(sim)),
                  par=as.character(par),
                  FS=paste0(NF1,'_',NF2))

### Points ---------------------------
pMu <- ggplot(RES1[RES1$par %in% c('Mu'),],aes(x=sim, y=value,color=FS)) +
  geom_point(cex=0.5)+
  geom_smooth(method = "loess",lwd=1.5) +
  xlab('Generation')+ylab(expression(mu))+
  facet_grid(GS_Model+Strategy~Intensity)

pGv <- ggplot(RES1[RES1$par %in% c('GV'),],aes(x=sim, y=value,color=FS)) +
  geom_point(cex=0.5)+
  geom_smooth(method = "loess",lwd=1.5) +
  xlab('Generation')+ylab(expression(sigma[g]^2))+
  facet_grid(GS_Model+Strategy~Intensity)

pHe <- ggplot(RES1[RES1$par %in% c('He'),],aes(x=sim, y=value,color=FS)) +
  geom_point(cex=0.5)+
  geom_smooth(method = "loess",lwd=1.5) +
  xlab('Generation')+ylab(expression(H^2))+
  facet_grid(GS_Model+Strategy~Intensity)

pAc <- ggplot(RES1[RES1$par %in% c('Accuracy'),],aes(x=sim, y=value,color=FS)) +
  geom_point(cex=0.5)+
  geom_smooth(method = "loess",lwd=1.5) +
  xlab('Generation')+ylab('Accuracy')+
  facet_grid(GS_Model+Strategy~Intensity)

pdf(paste0(JOBID,'_ResultsAllPoints.pdf'),w=20,h=120)
grid.arrange(pMu,pGv,pHe,pAc,nrow = 4)
dev.off()

