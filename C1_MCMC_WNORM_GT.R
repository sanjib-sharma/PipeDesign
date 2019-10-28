
#Initialize Data
rm(list=ls())
main_path="/gpfs/group/aim127/default/Sanjib/IDF/ACIS/data/stclg/IDF/StormSurge/"


infile<-paste('/gpfs/group/aim127/default/Sanjib/IDF/ACIS/data/stclg/IDF/AMS_stclg.txt',sep='')
big<-scan(infile,skip=0,list(date=0,depth=0))
depth<-big$depth # year 1896-2018
datset<-depth*25.4/24

infile<-paste('/gpfs/group/aim127/default/Sanjib/IDF/ACIS/data/stclg/IDF/Globaltemp.txt',sep='')
big<-scan(infile,skip=21,list(date=0,temp=0))
year<-big$date
tempset<-big$temp
###########################################################################################
# (1) Stationary
source("/gpfs/group/aim127/default/Sanjib/IDF/ACIS/data/stclg/IDF/StormSurge/SourceCode/Prior2SourceStat.R")
# Initial Conditions
start<-c(rep(1,3)) # Start Value

errvect<-c(0.0625, 0.03125, 0.06255) # Proposal Distribution SD
uvect<-c(rep(0,3)) # Prior Distribution Parameter 1: Does not matter for Uniform Priors
sdvect<-c(rep(100,3))  # Prior Distribution Parameter 2: Does not matter for Uniform Priors
#MCMC
run1<-MCMC_breaks_resume(50000,2500,res=FALSE,resdat=NA,datset,tempset,
                         start,errvect,uvect,sdvect,filetitle="Test1") # First Run
run2<-MCMC_breaks_resume(50000,2500,res=TRUE,resdat=run1,datset,tempset,
                         start,errvect,uvect,sdvect,filetitle="Test2")
# Diagnostics
# Run Diagnostics Again
BM_MCMC(run2,int=500,burn=0,bline=6000)
Mean_MCMC(run2,int=500,burn=0,bline=6000)
burnin=10000 # Burnin
MULTnsgevplots(run2,rm.burn=TRUE,burn=burnin) #Trace Plots Without Burn-In
MULTncrej(run2,burn=burnin) #RejectionRate
cred.table<-CredIntervalsGEV(run2,burn=burnin) # Credible Intervals
#save(run2,file="Prior2_Grinsted_Stat_Output.RData") # Save Data
#load("Prior2_Grinsted_Stat_Output.RData")
save(run2,file=paste(main_path,"/GTRESULT/stat_widenorm_run2_GT.RData",sep=""))
cred.table
save(cred.table,file=paste(main_path,"/GTRESULT/stat_widenorm_CI_GT.RData",sep=""))
#100-yr Return level
mu=run2$finmat[,1]
sigma=exp(run2$finmat[,2])
xi=run2$finmat[,3]
mu_chain <-mu[(length(xi)-40000+1):length(xi)]
sigma_chain <- sigma[(length(sigma)-40000+1):length(sigma)]
xi_chain <- xi[(length(xi)-40000+1):length(xi)]
retint<-benreturn(100,mu_chain,sigma_chain,xi_chain)
save(mu_chain,xi_chain,sigma_chain,file=paste(main_path,"/GTRESULT/stat_widenorm_param_GT.RData",sep=""))
save(retint,file=paste(main_path,"/GTRESULT/stat_widenorm_rtnlevel.RData",sep=""))

###########################################################################################
# (2) Mu Non-Stationary
source("/gpfs/group/aim127/default/Sanjib/IDF/ACIS/data/stclg/IDF/StormSurge/SourceCode/Prior2SourceMu.R")
# Initial Conditions
start<-c(rep(1,4)) # Start Value
errvect<-c(0.0625 ,0.03125, 0.0625 ,0.0625) # SD for Proposal 
uvect<-c(rep(0,4)) # Prior Distribution Parameter 1: Does not matter for Uniform Priors
sdvect<-c(rep(100,4))  # Prior Distribution Parameter 2: Does not matter for Uniform Priors
#MCMC
run1<-MCMC_breaks_resume(50000,2500,res=FALSE,resdat=NA,datset,tempset,
                         start,errvect,uvect,sdvect,filetitle="Test1") # First Run
run2<-MCMC_breaks_resume(50000,2500,res=TRUE,resdat=run1,datset,tempset,
                         start,errvect,uvect,sdvect,filetitle="Test2") # Resume 
# Diagnostics
BM_MCMC(run2,int=500,burn=0,bline=6000)
Mean_MCMC(run2,int=500,burn=0,bline=6000)
burnin=1000 # Burnin
MULTnsgevplots(run2,rm.burn=TRUE,burn=burnin) #Trace Plots Without Burn-In
MULTncrej(run2,burn=burnin) #RejectionRate
cred.table<-CredIntervalsGEV(run2,burn=burnin) # Credible Intervals
save(run2,file=paste(main_path,"/GTRESULT/Mu_nonstat_widenorm_run2_GT.RData",sep=""))
#save(run2,file="Prior2_Grinsted_Mu_Output.RData") # Save Data
#load("Prior2_Grinsted_Mu_Output.RData")
#source("/gpfs/group/aim127/default/Sanjib/IDF/ACIS/data/stclg/IDF/StormSurge/SourceCode/Prior2SourceMu.R")
cred.table
save(cred.table,file=paste(main_path,"/GTRESULT/Mu_nonstat_widenorm_CI_GT.RData",sep=""))
#100-yr Return level
mu=run2$finmat[,1]
sigma=exp(run2$finmat[,2])
xi=run2$finmat[,3]
mu_chain <-mu[(length(xi)-40000+1):length(xi)]
sigma_chain <- sigma[(length(sigma)-40000+1):length(sigma)]
xi_chain <- xi[(length(xi)-40000+1):length(xi)]
retint<-benreturn(100,mu_chain,sigma_chain,xi_chain)
save(mu_chain,xi_chain,sigma_chain,file=paste(main_path,"/GTRESULT/Mu_nonstat_widenorm_param_GT.RData",sep=""))
save(retint,file=paste(main_path,"/GTRESULT/Mu_nonstat_widenorm_rtnlevel.RData",sep=""))
###########################################################################################
# (3) Mu Sigma Non-Stationary
source("/gpfs/group/aim127/default/Sanjib/IDF/ACIS/data/stclg/IDF/StormSurge/SourceCode/Prior2SourceMuSigma.R")
# Initial Conditions
start<-c(rep(1,5)) # Start Value
errvect<-c(0.0625, 0.03125, 0.0625,0.0625, 0.25) # SD for Proposal
uvect<-c(rep(0,5)) # Prior Distribution Parameter 1: Does not matter for Uniform Priors
sdvect<-c(rep(100,5))  # Prior Distribution Parameter 2: Does not matter for Uniform Priors

##############################################################################################
#Testing for resuming MCMC
run1<-MCMC_breaks_resume(50000,2500,res=FALSE,resdat=NA,datset,tempset,
                         start,errvect,uvect,sdvect,filetitle="Test1") # First Run
run2<-MCMC_breaks_resume(50000,2500,res=TRUE,resdat=run1,datset,tempset,
                         start,errvect,uvect,sdvect,filetitle="Test2") # Resume 
# Run Diagnostics Again
BM_MCMC(run2,int=500,burn=0,bline=6000)
Mean_MCMC(run2,int=500,burn=0,bline=6000)
burnin=6000 # Burnin
MULTnsgevplots(run2,rm.burn=TRUE,burn=burnin) #Trace Plots Without Burn-In
MULTncrej(run2,burn=burnin) #RejectionRate
cred.table<-CredIntervalsGEV(run2,burn=burnin) # Credible Intervals

#load("Prior2_Grinsted_MuSigma_Output.RData")
burnin=6000 # Burnin
cred.table<-CredIntervalsGEV(run2,burn=burnin) # Credible Intervals
cred.table
save(run2,file=paste(main_path,"/GTRESULT/MuSigma_nonstat_widenorm_run2_GT.RData",sep=""))
save(cred.table,file=paste(main_path,"/GTRESULT/MuSigma_nonstat_widenorm_CI_GT.RData",sep=""))
#100-yr Return level
mu=run2$finmat[,1]
sigma=exp(run2$finmat[,2])
xi=run2$finmat[,3]
mu_chain <-mu[(length(xi)-40000+1):length(xi)]
sigma_chain <- sigma[(length(sigma)-40000+1):length(sigma)]
xi_chain <- xi[(length(xi)-40000+1):length(xi)]
retint<-benreturn(100,mu_chain,sigma_chain,xi_chain)
save(mu_chain,xi_chain,sigma_chain,file=paste(main_path,"/GTRESULT/MuSigma_nonstat_widenorm_param_GT.RData",sep=""))
save(retint,file=paste(main_path,"/GTRESULT/MuSigma_nonstat_widenorm_rtnlevel.RData",sep=""))

###########################################################################################
# (4) Fully Non-Stationary
source("/gpfs/group/aim127/default/Sanjib/IDF/ACIS/data/stclg/IDF/StormSurge/SourceCode/Prior2Source.R") # Load Source File
# Initial Conditions
start<-c(rep(1,6)) # Start Value
errvect<-c(0.0625, 0.03125, 0.0625,0.0625, 0.25 ,0.25) # Proposal SD
uvect<-c(rep(0,6)) # Prior Distribution Parameter 1: Does not matter for Uniform Priors
sdvect<-c(rep(100,6))  # Prior Distribution Parameter 2: Does not matter for Uniform Priors

#Testing for resuming MCMC

run1<-MCMC_breaks_resume(50000,2500,res=FALSE,resdat=NA,datset,tempset,
                         start,errvect,uvect,sdvect,filetitle="Test1") # First Run
run2<-MCMC_breaks_resume(50000,2500,res=TRUE,resdat=run1,datset,tempset,
                         start,errvect,uvect,sdvect,filetitle="Test2") # Resume 
# Run Diagnostics 
burnin=6000 # Burnin
MULTnsgevplots(run2,rm.burn=TRUE,burn=burnin) #Trace Plots Without Burn-In
MULTncrej(run2,burn=burnin) #RejectionRate
cred.table<-CredIntervalsGEV(run2,burn=burnin) # Credible Intervals
save(run2,file=paste(main_path,"/GTRESULT/Fully_nonstat_widenorm_run2_GT.RData",sep=""))
save(cred.table,file=paste(main_path,"/GTRESULT/Fully_nonstat_widenorm_CI_GT.RData",sep=""))
#100-yr Return level
mu=run2$finmat[,1]
sigma=exp(run2$finmat[,2])
xi=run2$finmat[,3]
mu_chain <-mu[(length(xi)-40000+1):length(xi)]
sigma_chain <- sigma[(length(sigma)-40000+1):length(sigma)]
xi_chain <- xi[(length(xi)-40000+1):length(xi)]
retint<-benreturn(100,mu_chain,sigma_chain,xi_chain)
save(mu_chain,xi_chain,sigma_chain,file=paste(main_path,"/GTRESULT/Fully_nonstat_widenorm_param_GT.RData",sep=""))
save(retint,file=paste(main_path,"/GTRESULT/Fully_nonstat_widenorm_rtnlevel.RData",sep=""))

###########################################################################################
