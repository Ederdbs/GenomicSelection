### Simulation Breeding Selection Strategy, Index of Selection and GS Models ---
### April 2021
### coded by Eder Silva and Alencar Xavier
### Load packagers -------------------------------------------------------------
rm(list=ls())
Xpackagers <- c('AlphaSimR','bWGR','parallel','foreach','doParallel',
                'reshape','ggplot2','gridExtra','lubridate','plyr',
                'ranger','Rcpp','keras','verification','rrBLUP',
                'reshape2','ScottKnott','viridis')
XXX <- lapply(Xpackagers, function(x){suppressMessages(require(x,quietly = TRUE, character.only = TRUE))})

### Simulation Parameters  -----------------------------------------------------
Number_of_runs = 2 #60 number of replicated scenario 
Number_of_generations = 5 #200 number of breeding cycles  

### Simulation settings --------------------------------------------------------
Intensity <- c(2.5, #just select 2 intensity for run more fast
               #5,
               #7.5,
               10)/100

FS <- list('1'=c(30,20), #F1 and F2 size #just select 2 families size for run more fast
          # '2'=c(250,60), #F1 and F2 size
          # '3'=c(200,75), #F1 and F2 size
         #  '4'=c(150,100), #F1 and F2 size
           '5'=c(10,60)) #F1 and F2 size

GS_Model <- list( #list of genomic selection models
  #'RF' = function(y,gen,...){list(hat=ranger::ranger(y~.,data.frame(y=y,gen),verbose = FALSE, save.memory = TRUE,write.forest = FALSE,...)$predictions)},
  'GBLUP'=emML,
  #'RR'=emRR,
   'BayesA'=emBA,
  #'BayesB'=emBB,
  #'BayesC'=emBC,
  #'BayesL'=emBL,
  #'FLM' = emDE,
  'Random' = function(y,gen){return(list(hat=sample(y)))},
  #'Gv'= function(y,gen){return(list(hat=NA))},
  'Pheno' = function(y,gen){return(list(hat=NA))}
)
Sel_Strategy <- as.character(c('WPSF', # Within pre-selected families
                               'WF', # Within family
                               'AF'  # Across Family
))

### Genetics parameters --------------------------------------------------------
POPBestIntensity = 0.3 ## percentage of selection the best population 
NCgs = 3    # Breeding cycles used in estimation set
Ne = 106   # Population effective size for Mac
segSites  = 1000 #segregate site by chromosome 
nInd = 200 # number of individual available to star simulation
nSnpPerChr <-  300 # SNP per chromosome  
nQtlPerChr <- segSites * 0.7 # Number of Qtl per chromosome
### Yield Trait parameter ------------------------------------------------------
mean = 60 # From SoyNAN 
var = 77  # From SoyNAN 
varGxE = 77 # From SoyNAN 
varEnv= 200 # From SoyNAN 
corA = matrix(1,nrow=1)
corGxE = matrix(1,nrow=1)
nReps=1

### Run Simulation ------------------------------------------------------------- 
source('GS_Function.R')

### Run Analysis of results ---------------------------------------------------- 
BGenerations <- c(2,3,5) ### Breeding cycle to breakout
source('GS_Results.R')
# see output files in folder
### the end.....