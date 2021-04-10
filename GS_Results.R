### Analysis of simulated data -------------------------------------------------
MYC <- c("blue",'yellow',"darkorange","blueviolet","green","red",'aquamarine',"darkgreen")
Generations <- paste("Mu",BGenerations,sep='')
# load data   ------------------------------------------------------------------
JOBID <- 'ALL'
RES <- ldply(lapply(list.files(pattern='Results.csv$'),function(x){read.csv(x)}))

### Change Accuracy of pheno to accuracy GvPhe
RES[RES$GS_Model == 'Pheno',grep('Accuracy[0-9]',colnames(RES))] <- RES[RES$GS_Model == 'Pheno',grep('AccuracyGvPhe[0-9]',colnames(RES))]
RES[RES$GS_Model == 'Gv',grep('Accuracy[0-9]',colnames(RES))] <- 1

RES[is.na(RES)] <- 0
RES <- cbind(S=1:nrow(RES),RES[,-1])
RES <- RES[RES$GS_Model %in% c("BayesA","BayesB","FLM","GBLUP","Gv","Pheno","Random","RF"),]
RES$GS_Model <- as.character(RES$GS_Model)
RES[RES$GS_Model=='Gv','GS_Model'] <- 'TBV'
RES$Intensity <- RES$Intensity*100

RES[,grep('^Mu',colnames(RES))] <- RES[,grep('^Mu',colnames(RES))]*60/1000
write.csv(RES,paste0(JOBID,'_o_ResultsALL.csv'),row.names = F)

RES$Strategy <- as.character(RES$Strategy)
RES$Strategy <- factor(RES$Strategy,levels = c('AF','WPSF','WF'))

RES$Selection_Intensity <- paste(sprintf('%.1f',RES$Intensity),'%')
RES$Selection_Intensity <- factor(RES$Selection_Intensity,levels=c("2.5 %","5.0 %","7.5 %","10.0 %"))

RES1 <- melt(RES,id=c('S','GS_Model','Strategy','Intensity','Selection_Intensity'))
RES1$par <- gsub('[0-9]','',RES1$variable)
RES1$sim <- gsub('[A-z]','',RES1$variable)
RES1 <- transform(RES1,
                  value=as.numeric(as.character(value)),
                  sim=as.numeric(as.numeric(sim)),
                  par=as.character(par))
write.csv(RES1,paste0(JOBID,'_o_ResultsALLMelt.csv'),row.names = F)
save.image('ALLData.RData')

### Gain versus selection intensity --------------------------------------------
pdf('ALL_Gain_x_SI_1.pdf',w=15,h=10)
res <- list();con <- 1
for(t in Generations){
  a <- ggplot(RES,aes(x=Intensity, y=RES[,t],color=GS_Model)) +
    geom_point(cex=0.5)+
    geom_smooth(method='lm',se=FALSE,lwd=2)+
    xlab('Intensity Index')+
    ylab(t)+
    facet_grid(~Strategy, labeller = label_both)+
    theme_bw()+
    scale_colour_manual(values = MYC )+
    guides(color=guide_legend("Models:",nrow=1,byrow=TRUE))+
    theme(legend.position = "bottom", legend.box = "horizontal")
  print(a)
  for(i in unique(RES$GS_Model)){
    for(j in unique(RES$Strategy)){
      dataS <- RES[RES$GS_Model==i & RES$Strategy == j,c(t,'Intensity')]
      m0 <- lm(dataS[,1]~dataS[,2])
      res[[con]] <- c(t,i,j,round(summary(m0)$coef[2,],3),R2=round(summary(m0)$r.squared,2))
      con <- con+1
    }
  }
  
}
dev.off()
rest <- ldply(res)
rest$V1 <- as.numeric(as.character(gsub('Mu','',rest$V1)))
rest$coeff <- paste0(sprintf('%.2f',as.numeric(rest$Estimate)),' (sd:',
                     sprintf('%.2f',as.numeric(rest$`Std. Error`)),' p-value:',
                     sprintf('%.2f',as.numeric(rest$`Pr(>|t|)`)),' R²:',
                     sprintf('%.2f',as.numeric(rest$R2)),')')
rest1 <- melt(rest[,c("V1","V2","V3","coeff")], id.vars = c("V1","V2","V3"))

res <- list();con <- 1
for(t in Generations){
  dataS1 <- RES[,c(t,'Intensity','GS_Model','Strategy')]
  dataS1$GSS <- as.factor(paste0(dataS1$GS_Model,'_',dataS1$Strategy))
  dataS1$Int <- as.factor(dataS1$Intensity)
  
  out <- ldply(lapply(1:length(levels(dataS1$Int)),
                      function(x){summary(SK.nest(aov(dataS1[,t] ~ Int*GSS,data=dataS1),
                                                  which='Int:GSS',
                                                  fl1=x))}))
  out <- cbind(ldply(strsplit(as.character(out$Levels),'/')),out)
  out1 <- ldply(lapply(1:length(levels(dataS1$GSS)),
                       function(x){summary(SK.nest(aov(dataS1[,t] ~ Int*GSS,data=dataS1),
                                                   which='GSS:Int',
                                                   fl1=x))}))
  out1 <- cbind(ldply(strsplit(as.character(out1$Levels),'/')),out1)
  out1$`SK(5%)` <- toupper(out1$`SK(5%)`)
  out$key <- paste0(out$V2,out$V1)
  out1$key <- paste0(out1$V1,out1$V2)
  outf <- merge(out,out1,by.x='key',by.y='key')
  outf$value <- paste0(sprintf('%.1f',outf$Means.x)," ",outf$`SK(5%).x`,outf$`SK(5%).y`)
  res[[con]] <- cbind(t,outf[,c('key','V1.x','V2.x','value')])
  con <- con+1
  
}
restSK <- ldply(res)
restSK$t <- as.numeric(as.character(gsub('Mu','',restSK$t)))
restSK <- cbind(restSK,ldply(strsplit(restSK$V2.x,'_')))
rest1SK <- melt(restSK[,c('t','V1','V2','V1.x','value')], id.vars = c('t',"V1","V2","V1.x"))
rest1$V1.x <- 'XX'
rest1SK <- rest1SK[c('t','V1','V2','variable','value','V1.x')]
colnames(rest1SK) <- colnames(rest1)
res_all <- rbind(rest1,rest1SK)
rest_allx <- dcast(res_all, V2+V3 ~ V1 + V1.x + variable)
write.csv(rest_allx,'Gain_x_SI.csv',row.names = FALSE)

## Estimable to lose 80% genetic variance --------------------------------------
RESGV <- RES[,c(1:4,grep('Selection_Intensity',colnames(RES)),grep('GV',colnames(RES)))]
RESGV$GVres80 <- NA
for(i in 1:nrow(RESGV)){
  RESGV$GVres80[i] <- sum(RESGV[i,(7:ncol(RESGV)-1)] >= (RESGV[i,7]*0.2))
}
RESGV <- RESGV[,c('S','GS_Model','Strategy','Selection_Intensity','GVres80')]
a <- ggplot(RESGV,aes(x=GS_Model, y=GVres80)) +
  geom_boxplot()+
  xlab('GS Models')+
  ylab('Generation to lose 80% Genetic variance')+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  facet_grid(Strategy~Selection_Intensity, labeller = label_both)+
  theme_bw()+
  scale_colour_manual(values = MYC )+
  guides(color=guide_legend("Models:",nrow=1,byrow=TRUE))+
  theme(legend.position = "bottom", legend.box = "horizontal")

jpeg('ALL_Residual80GV1.jpg',width = 250, height = 180 ,units='mm',quality=95,res=900)
print(a)
dev.off()

a <- ggplot(RESGV[RESGV$GS_Model %in% c('BayesA','GBLUP','Pheno'),],aes(x=GS_Model, y=GVres80,colour=GS_Model)) +
  geom_boxplot()+
  xlab('GS Models')+
  ylab('Generation to lose 80% Genetic variance')+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  facet_grid(Strategy~Selection_Intensity, labeller = label_both)+
  theme_bw()+
  scale_colour_manual(values = MYC )+
  guides(color=guide_legend("Models:",nrow=1,byrow=TRUE))+
  theme(legend.position = "bottom", legend.box = "horizontal")

jpeg('ALL_Residual80GV2.jpg',width = 250, height = 180 ,units='mm',quality=95,res=900)
print(a)
dev.off()

Residual80GV <- cast(RESGV,GS_Model+Strategy ~ Selection_Intensity, 
                     fun.aggregate =function(x){return(paste(sprintf('%.2f',mean(x)),
                                                             '(',
                                                             sprintf('%.2f',sd(x)),')',sep=''))},
                     value='GVres80')
write.csv(Residual80GV,'Residual80GV.csv',row.names = FALSE)

RESGV$Intensity <- as.factor(RESGV$Selection_Intensity)
RESGV$GSS <- as.factor(paste0(RESGV$GS_Model,'_',RESGV$Strategy))
out <- ldply(lapply(1:length(unique(RESGV$Intensity)),
                    function(x){summary(SK.nest(aov(RESGV$GVres80 ~ Intensity*GSS,data=RESGV),
                                                which='Intensity:GSS',
                                                fl1=x))}))
out <- cbind(ldply(strsplit(as.character(out$Levels),'/')),out)
out1 <- ldply(lapply(1:length(unique(RESGV$GSS)),
                     function(x){summary(SK.nest(aov(RESGV$GVres80 ~ Intensity*GSS,data=RESGV),
                                                 which='GSS:Intensity',
                                                 fl1=x))}))
out1 <- cbind(ldply(strsplit(as.character(out1$Levels),'/')),out1)
out1$`SK(5%)` <- toupper(out1$`SK(5%)`)
out$key <- paste0(out$V2,out$V1)
out1$key <- paste0(out1$V1,out1$V2)
outf <- merge(out,out1,by.x='key',by.y='key')
outf$value <- paste0(sprintf('%.0f',outf$Means.x)," ",outf$`SK(5%).x`,outf$`SK(5%).y`)
outf <- outf[,c('key','V1.x','V2.x','value')]
restSK <- cbind(outf,ldply(strsplit(outf$V2.x,'_')))
rest1SK <- melt(restSK[,c('V1.x','V1','V2','value')], id.vars = c("V1","V2","V1.x"))
rest_allx <- dcast(rest1SK, V1+V2 ~ V1.x + variable)
write.csv(rest_allx ,'Residual80GVSK.csv',row.names = FALSE)

## Estimable to lose 80% Acurracy -------------------------------------
RESAC <- RES[,c(1:4,grep('Selection_Intensity',colnames(RES)),grep('Accuracy[0-9]',colnames(RES)))]
RESAC$ACres80 <- NA
for(i in 1:nrow(RESAC)){
  RESAC$ACres80[i] <- sum(RESAC[i,(7:ncol(RESAC))-1] >= (RESAC[i,7]*0.2))
}
RESAC <- RESAC[,c('S','GS_Model','Strategy','Selection_Intensity','ACres80')]
a <- ggplot(RESAC,aes(x=GS_Model, y=ACres80)) +
  geom_boxplot()+
  xlab('GS Models')+
  ylab('Generation to lose 80% Acurracy')+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  facet_grid(Strategy~Selection_Intensity, labeller = label_both)+
  theme_bw()+
  scale_colour_manual(values = MYC )+
  guides(color=guide_legend("Models:",nrow=1,byrow=TRUE))+
  theme(legend.position = "bottom", legend.box = "horizontal")

jpeg('ALL_Residual80AC1.jpg',width = 250, height = 180 ,units='mm',quality=95,res=900)
print(a)
dev.off()

a <- ggplot(RESAC[RESAC$GS_Model %in% c('BayesA','GBLUP','Pheno'),],aes(x=GS_Model, y=ACres80,colour=GS_Model)) +
  geom_boxplot()+
  xlab('GS Models')+
  ylab('Generation to lose 80% Acurracy')+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  facet_grid(Strategy~Selection_Intensity, labeller = label_both)+
  theme_bw()+
  scale_colour_manual(values = MYC )+
  guides(color=guide_legend("Models:",nrow=1,byrow=TRUE))+
  theme(legend.position = "bottom", legend.box = "horizontal")

jpeg('ALL_Residual80AC2.jpg',width = 250, height = 180 ,units='mm',quality=95,res=900)
print(a)
dev.off()

Residual80AC <- cast(RESAC,GS_Model+Strategy ~ Selection_Intensity, 
                     fun.aggregate =function(x){return(paste(sprintf('%.2f',mean(x)),
                                                             '(',
                                                             sprintf('%.2f',sd(x)),')',sep=''))},
                     value='ACres80')
write.csv(Residual80AC,'Residual80AC.csv',row.names = FALSE)

RESAC$Intensity <- as.factor(RESAC$Selection_Intensity)
RESAC$GSS <- as.factor(paste0(RESAC$GS_Model,'_',RESAC$Strategy))
out <- ldply(lapply(1:length(unique(RESAC$Intensity)),
                    function(x){summary(SK.nest(aov(RESAC$ACres80 ~ Intensity*GSS,data=RESAC),
                                                which='Intensity:GSS',
                                                fl1=x))}))
out <- cbind(ldply(strsplit(as.character(out$Levels),'/')),out)
out1 <- ldply(lapply(1:length(unique(RESAC$GSS)),
                     function(x){summary(SK.nest(aov(RESAC$ACres80 ~ Intensity*GSS,data=RESAC),
                                                 which='GSS:Intensity',
                                                 fl1=x))}))
out1 <- cbind(ldply(strsplit(as.character(out1$Levels),'/')),out1)
out1$`SK(5%)` <- toupper(out1$`SK(5%)`)
out$key <- paste0(out$V2,out$V1)
out1$key <- paste0(out1$V1,out1$V2)
outf <- merge(out,out1,by.x='key',by.y='key')
outf$value <- paste0(sprintf('%.0f',outf$Means.x)," ",outf$`SK(5%).x`,outf$`SK(5%).y`)
outf <- outf[,c('key','V1.x','V2.x','value')]
restSK <- cbind(outf,ldply(strsplit(outf$V2.x,'_')))
rest1SK <- melt(restSK[,c('V1.x','V1','V2','value')], id.vars = c("V1","V2","V1.x"))
rest_allx <- dcast(rest1SK, V1+V2 ~ V1.x + variable)
write.csv(rest_allx ,'Residual80ACSK.csv',row.names = FALSE)

### Genetic drift -------------------------------------------------------------
GenDR <- RES1[RES1$par == 'GV' &
                RES1$GS_Model == 'Random',]
a <- ggplot(GenDR,aes(x=sim, y=value)) +
  geom_point(cex=0.5)+
  geom_smooth(method='lm')+
  xlab('Intensity Index')+
  facet_grid(Strategy~Intensity, labeller = label_both)+
  theme_bw()+
  scale_colour_manual(values = MYC )+
  guides(color=guide_legend("Models:",nrow=1,byrow=TRUE))+
  theme(legend.position = "bottom", legend.box = "horizontal")

jpeg('ALL_GeneticDrifft.jpg',width = 250, height = 180 ,units='mm',quality=95,res=900)
print(a)
dev.off()

res <- list();con <- 1  
for(i in unique(RES$Intensity)){
  for(j in unique(RES$Strategy)){
    dataS <- GenDR[GenDR$Intensity==i & GenDR$Strategy == j,]
    m0 <- lm(value~sim,data=dataS)
    res[[con]] <- c(i,j,round(summary(m0)$coef[2,],3),round(summary(m0)$r.squared,3))
    con <- con+1
  }
}
rest <- ldply(res)
rest <- rest[,!(colnames(rest) %in% 't value')]
for(i in 3:ncol(rest)){rest[,i] <- as.numeric(as.character(rest[,i]))}
rest1 <- melt(rest, id.vars = c("V1","V2"))
rest1 <- dcast(rest1, V2 + V1 ~ variable)
colnames(rest1) <- c('Strategy','SI','Estimate','StdError','Pvalue','R^2')
write.csv(rest1,'GeneticDrifft.csv',row.names = FALSE)


rest1$SI <- factor(as.character(rest1$SI),levels = c( '2.5',  '5',  '7.5', '10'))

jpeg(paste0('Final_VarianceDrifft_Results.jpg'),width = 150, height = 100 ,units='mm',quality=95,res=900)
a <- ggplot(rest1,aes(x=Strategy, y=Estimate,fill=SI)) +
  geom_bar(stat="identity", position=position_dodge())+
  geom_errorbar(aes(ymin=Estimate- StdError, ymax=Estimate+ StdError), width=.2,
                position=position_dodge(.9))+
  xlab('Breeding Strategy')+
  ylab(expression("%"~sigma[g]^{2}~"loss per cycle"))+
  theme_bw()+
  scale_colour_manual(values = MYC )+
  guides(color=guide_legend("Models:",nrow=1,byrow=TRUE))
print(a)
dev.off()

### Accuracy and genetic variance ----------------------------------------------
RES_AGV <- cbind(RES[,c('S','GS_Model','Strategy','Intensity','Selection_Intensity')],
                 ACbyGV=RES[,grep('Accuracy[0-9]$',colnames(RES))]/sqrt(RES[,grep('GV[0-9]$',colnames(RES))])
)

RES_AGVM <- melt(RES_AGV,id=c('S','GS_Model','Strategy','Intensity','Selection_Intensity'))
RES_AGVM$sim <- gsub('[A-z.]','',RES_AGVM$variable)
RES_AGVM <- transform(RES_AGVM,
                      value=as.numeric(as.character(value)),
                      sim=as.numeric(as.numeric(sim)))

a <- ggplot(RES_AGVM,aes(x=sim, y=value,color=GS_Model)) +
  geom_smooth(method = "loess",lwd=1.25,se=FALSE,span = .4) +
  xlab('Breeding cycles')+
  ylab(' Accuracy by Genetic variance')+
  #ylim(c(0,0.25))+
  facet_grid(Strategy~Selection_Intensity, labeller = label_both)+
  theme_bw()+
  scale_colour_manual(values = MYC )+
  guides(color=guide_legend("Models:",nrow=1,byrow=TRUE))+
  theme(legend.position = "bottom", legend.box = "horizontal")


jpeg('Final_Accuracy_by_GeneticVariance.jpg',width = 250, height = 180 ,units='mm',quality=95,res=900)
print(a)
dev.off()

jpeg(paste0('Final_mu_Results.jpg'),width = 250, height = 180 ,units='mm',quality=95,res=900)
a <- ggplot(RES1[RES1$par %in% 'Mu',],aes(x=sim, y=value,color=GS_Model)) +
  geom_smooth(method = "loess",lwd=0.8,se=FALSE,span = 1)+
  xlab('Breeding cycles')+
  ylab(expression('Population mean'~(t.ha^{-1})))+
  #ylim(c(0,53))+
  facet_grid(Strategy~Selection_Intensity, labeller = label_both)+
    theme_bw()+
  scale_colour_manual(values = MYC )+
  guides(color=guide_legend("Models:",nrow=1,byrow=TRUE))+
  theme(legend.position = "bottom", legend.box = "horizontal")
print(a)
dev.off()

jpeg(paste0('Final_variance_Results.jpg'),width = 250, height = 180 ,units='mm',quality=95,res=900)
a <- ggplot(RES1[RES1$par %in% 'GV',],aes(x=sim, y=value,color=GS_Model)) +
  geom_smooth(method = "loess",lwd=0.8,se=FALSE,span = 0.4) +
  xlab('Breeding cycles')+
  ylab(expression(sigma[g]^{2}))+
  ylim(c(0,120))+
  facet_grid(Strategy~Selection_Intensity, labeller = label_both)+
  theme_bw()+
  scale_colour_manual(values = MYC )+
  guides(color=guide_legend("Models:",nrow=1,byrow=TRUE))+
  theme(legend.position = "bottom", legend.box = "horizontal")
print(a)
dev.off()

jpeg(paste0('Final_Accuracy_Results.jpg'),width = 250, height = 180 ,units='mm',quality=95,res=900)
a <- ggplot(RES1[RES1$par %in% 'Accuracy',],aes(x=sim, y=value,color=GS_Model)) +
  geom_smooth(method = "loess",lwd=0.8,se=FALSE,span = 1) +
  xlab('Breeding cycles')+
  ylab(expression(Accuracy))+
  ylim(c(0,1.05))+
  facet_grid(Strategy~Selection_Intensity, labeller = label_both)+
  theme_bw()+
  scale_colour_manual(values = MYC )+
  guides(color=guide_legend("Models:",nrow=1,byrow=TRUE))+
  theme(legend.position = "bottom", legend.box = "horizontal")
print(a)
dev.off()



### Analysis of mean by regression ---------------------------------------------
RES_sel <- RES1[RES1$par %in% c('Mu') & !(RES1$GS_Model %in% c('Random')),]
RES_sel$GS_Model <- as.factor(RES_sel$GS_Model)
res <- list()
dat_out <- data.frame()
con <- 1
for(GS in unique(RES_sel$GS_Model)){
  for(ST in unique(RES_sel$Strategy)){
    bsel <- dplyr::filter(RES_sel,GS_Model == GS,
                          Strategy == ST)
    n1 <- nls(formula=value~A*(1-exp(-V*sim))+D*Intensity, data=bsel,
              start=list(A=500, V=1,D=1), trace=FALSE)
    n1out <- summary(n1)
    n2 <- nls(formula=value~A*(1-exp(sim*log(1-0.8)/V))+D*Intensity, data=bsel,
              start=list(A=500, V=1,D=1), trace=FALSE)
    n2out <- summary(n2)
    ### R2 Modelo conjunto
    R2 <- 1-deviance(n1)/deviance(lm(value~1, data=bsel))
    res[[con]] <- cbind(GS,
                        ST,
                        PAR=c('A','V','D','Vstar'),
                        rbind(n1out$parameters,n2out$parameters['V',]),                    
                        R2)
    con <- con+1
    dat <- expand.grid(GS_Model=GS,Strategy=ST,sim=1:200,Intensity=seq(2.5,10,l=30))
    dat$pred <- predict(n1, list(sim = dat$sim,Intensity=dat$Intensity))
    dat_out <- rbind(dat_out,dat) 
  }
}
res_out <- ldply(res)
for(i in 4:ncol(res_out)){
  res_out[,i] <- as.numeric(as.character(res_out[,i]))
}

res_out$value <- paste0(sprintf("%.3f",res_out$Estimate),
                       '(',
                       sprintf("%.3f",res_out$`Std. Error`),
                       ',',
                       sprintf("%.3f",res_out$`Pr(>|t|)`),
                       ')')
head(res_out)
res_out1 <- dcast(res_out, GS + ST + R2 ~ PAR,value='value')
write.csv(res_out1,paste0(JOBID,'_ResultsRegression.csv'),row.names = F)
mycol <- c("red","yellow","green","darkgreen")
jpeg(paste0('Final_regression_Model.jpg'),width = 250, height = 150 ,units='mm',quality=95,res=900)
ggplot(data = dat_out) +
  theme_bw() +
  geom_tile(aes(x = sim, y = Intensity, fill = pred)) +
  scale_fill_gradientn(colours = mycol,limits=c(2,max(dat_out$pred)))+
  facet_grid(Strategy~GS_Model, labeller = label_both)+
  xlab('Breeding cycles')+
  ylab('Selection intensity (%)')
dev.off()
