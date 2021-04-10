### Functions for Simulation Breeding Selection Strategy, Index of Selection and GS Models ###
### Generate Founders populations ----------------------------------------------
#set.seed(123)
founderPop <- quickHaplo(nInd=nInd,
                         nChr=20,
                         segSites=segSites,
                         genLen = 1.15,
                         ploidy = 2L,
                         inbred = TRUE)
## Founder POP and Simulation Parameters ----------------------------------------
FS <- ldply(FS)
colnames(FS) <- c('FS','F1','F2')
POSSI <- expand.grid(Sel_Strategy=as.character(Sel_Strategy),
                     run=1:Number_of_runs,
                     Intensity= Intensity,
                     Dmodel=1:length(GS_Model),
                     FS=FS$FS)
POSSI <- merge(POSSI,FS,all.x=TRUE)
cat(paste('Number of Simulation: ',nrow(POSSI),'\n'))
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
# function for main simulation -------------------------------------------------
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
      if(Strategy == 'WPSF'){ #WithinBestFamily POPBestIntensity (0.3)
        pop = selectFam(pop,nFam=round(NF1*POPBestIntensity),use=useM,simParam=SP)
        pop = selectWithinFam(pop,nInd=round((NF2*Int)/POPBestIntensity),use=useM,simParam=SP)
        nIndSel = pop@nInd
        pop = selectCross(pop, nInd=pop@nInd, nCrosses=NF1,use=useM, nProgeny=1, simParam=SP)
      }
      if(Strategy == 'WF'){ #WithinFamily
        pop = selectWithinFam(pop,nInd=round(NF2*Int),use=useM,simParam=SP)
        nIndSel = pop@nInd
        pop = selectCross(pop, nInd=pop@nInd, nCrosses=NF1,use=useM, nProgeny=1, simParam=SP)
      } 
      if(Strategy == 'AF'){ #AcrossFamily
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
