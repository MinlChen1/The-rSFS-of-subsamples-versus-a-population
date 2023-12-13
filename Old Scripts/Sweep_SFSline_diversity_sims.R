#functions to estimate diversity using the SFS/line data (option -o 92 in mstatspop) in relation to the total population

#structure of the input data: infile:	./chr16.DUROC_100.RAW.ID.tfa.gz	scaffold_name:	16	start_window:	250212	end_window:	322075	missing:	0	iteration:	0	npermutations:	0	seed:	123456	Length:	71864.00	Lengtht:	71864	mh:	0	Ratio_S/V:	2.082	Ratio_Missing:	NA	Variants:	1279	popsize:	1	nsam[0]:	100	Eff_length1_pop[0]:	71864.00	Eff_length2_pop[0]:	71864.00	Eff_length3_pop[0]:	71864.00	
#linexfreqy_pop[0]: .. HERE there is a vector divided in line (rows) x freqs (columns): in this case is 100 lines and 50 diff folded frequencies (no fixed included)
#field 16 (nsam[0]:) and field 20 (linexfreqy_pop[0]:) are necessary. to estimate diversity per position the field 18 (Eff_length2_pop[0]:) is also necessary.

#Functions to calculate Theta Watterson (1979) Theta Tajima (1983) and Theta FuLi (1993) from folded SFS using the Achaz approach (2009)

deltak <- function(i,j) {
  if(i==j) return(1)
  return(0)
}

weight.watt.unfolded <-function(nsam) {
  w <- array(0,dim=c(floor(nsam)))
  for(i in 1:length(w)) {
    w[i] <- 1/i
  }
  w
}

weight.taj.unfolded <-function(nsam) {
  w <- array(0,dim=c(floor(nsam)))
  for(i in 1:length(w)) {
    w[i] <- nsam-i
  }
  w
}

weight.fuli.unfolded <-function(nsam) {
  w <- array(0,dim=c(floor(nsam)))
  w[1] <- 1
  w
}

weight.fw.unfolded <- function(nsam) {
  w <- array(0,dim=c(floor(nsam)))
  for(i in 1:length(w)) {
    w[i] <- i
  }
  w
}

weight.watt.folded <-function(nsam) {
  w <- array(0,dim=c(floor(nsam/2)))
  for(i in 1:length(w)) {
    w[i] <- nsam/(i*(nsam-i)*(1+deltak(i,nsam-i)))
  }
  w
}

weight.taj.folded <-function(nsam) {
  w <- array(0,dim=c(floor(nsam/2)))
  for(i in 1:length(w)) {
    w[i] <- nsam/(1+deltak(i,nsam-i))
  }
  w
}

weight.fuli.folded <-function(nsam) {
  w <- array(0,dim=c(floor(nsam/2)))
  w[1] <- nsam
  w
}

phi.i <- function(nsam) {
  w <- array(0,dim=c(floor(nsam/2)))
  for(i in 1:length(w)) { 
    w[i] <- (nsam/(i*(nsam-i))) / (1+deltak(i,nsam-i))
  }
  w
}

psi.ij <- function(nsam,subnsam) {
  w <- array(0,dim=c(floor(nsam)))
  for(i in 1:length(w)) { 
    w[i] <- (1-dhyper(x=0,k=subnsam,m=i,n=nsam-i))
  }
  w
}

psi.fold.ij <- function(nsam,subnsam) {
  w <- array(0,dim=c(floor(nsam/2)))
  for(i in 1:length(w)) { 
    w[i] <- (1-dhyper(x=0,k=subnsam,m=i,n=nsam-i))
  }
  w
}

Calc.Theta.unfolded <- function(sfs,w,psi) {
  th <- 0
  for(i in 1:length(sfs)) {
    th <- th + w[i] * i * sfs[i] * 1/psi[i]
  }
  th <- th/(sum(w))
}

Calc.Theta.folded <- function(sfs,w,phi,psi) {
  th <- 0
  for(i in 1:length(sfs)) {
    th <- th + w[i] * sfs[i] * 1/phi[i] * 1/psi[i]
  }
  th <- th/(sum(w))
}

if (!require("ggplot2")) install.packages("ggplot2") 
library(ggplot2)
if (!require("akima")) install.packages("akima") 
library(akima)
if (!require("reshape2")) install.packages("reshape2") 
library(reshape2)

######################################################################
args = commandArgs(trailingOnly=TRUE)
if(length(args)<4) {
  popsize <- 2
  nsam <- 100
  niter <- 10
  path_file <- "~/Desktop/YuliDUROC-rSFS/testing_simulations/Sweeps/Slim-INDIVIDUALS"
} else {
  popsize <- as.numeric(args[1])
  nsam <- as.numeric(args[2])
  niter <- as.numeric(args[3])
  path_file <- args[4]
}
nind <- nsam/popsize
npops <- nsam/popsize

######################################################################
######################################################################
#
# UNFOLDED
#
#ESTIMATE RELATIVE NUCLEOTIDE DIVERSITY PER POP (ind) (that is, IN RELATION TO TOTAL, option -o 92)

model <- c("SNM","Complete_Sel_Sweep","Incomplete_Sel_Sweep","Standing_Sel_Sweep")
if(popsize==2) data.names  <- c("slim.snm.output_file.ms_mstatspopo92.txt","slim.selsweep.output_file.ms_mstatspopo92.txt","slim.incomplete_selsweep.output_file.ms_mstatspopo92.txt","slim.standing_selsweep.output_file.ms_mstatspopo92.txt")
if(popsize==1) data.names  <- c("slim.snm.output_file.ms_hap_mstatspopo92.txt","slim.selsweep.output_file.ms_hap_mstatspopo92.txt","slim.incomplete_selsweep.output_file.ms_hap_mstatspopo92.txt","slim.standing_selsweep.output_file.ms_hap_mstatspopo92.txt")
data.names1 <- c("slim.snm.output_file.ms_mstatspopo1.txt","slim.selsweep.output_file.ms_mstatspopo1.txt","slim.incomplete_selsweep.output_file.ms_mstatspopo1.txt","slim.standing_selsweep.output_file.ms_mstatspopo1.txt")

pdf(sprintf("%s/Plot_sweeps_rSFS_unfolded_popsize%.0f.pdf",path_file,popsize))
fil <- 3
for(fil in c(1:4)) {
  #define data frames for sending results
  columns.fuli <- c("NAME","nsam","Length",sprintf("ThetaRFuLi_%0.f",c(1:npops)))
  data.results.Fuli <- data.frame(matrix(nrow=0,ncol=length(columns.fuli)))
  colnames(data.results.Fuli) <- columns.fuli
  
  columns.watt <- c("NAME","nsam","Length",sprintf("ThetaRWatt_%0.f",c(1:npops)))
  data.results.Watt <- data.frame(matrix(nrow=0,ncol=length(columns.watt)))
  colnames(data.results.Watt) <- columns.watt
  
  columns.taj  <- c("NAME","nsam","Length",sprintf("ThetaRTaj_%0.f",c(1:npops)))
  data.results.Taj  <- data.frame(matrix(nrow=0,ncol=length(columns.taj)))
  colnames(data.results.Taj) <- columns.taj
  
  columns.fw  <- c("NAME","nsam","Length",sprintf("ThetaRFW_%0.f",c(1:npops)))
  data.results.fw  <- data.frame(matrix(nrow=0,ncol=length(columns.fw)))
  colnames(data.results.fw) <- columns.fw
  
  #data rSFS
  data.freqs <- read.table(file=sprintf("%s/%s",path_file,data.names[fil]))
  
  w.fuli <- weight.fuli.unfolded(nsam)
  w.watt <- weight.watt.unfolded(nsam)
  w.taj  <- weight.taj.unfolded(nsam)
  w.fayw <- weight.fw.unfolded(nsam)
  
  w.psi <- psi.ij(nsam,popsize)
  
  #for each gene/iteration extract the sfs per lineage and calculate their theta statistics
  init.mat <- grep("rSFS",data.freqs[1,])
  num.gene <- 1
  for(num.gene in 1:dim(data.freqs)[1]) {
    p.len.seq <- grep("Length:",data.freqs[num.gene,]) + 1
    sfs.line <- matrix(as.numeric(data.freqs[num.gene,((init.mat+1):(init.mat+(npops*floor(nsam-1))))]),byrow=T,nrow=npops,ncol=(floor(nsam-1)))
    
    #plot(x=c(1:(nsam-1)),y=sfs.line[1,],type="l")
    #for(f in 1:npops) lines(x=c(1:39),y=sfs.line[f,],type="l")
    
    Theta.FuLi.line <- round(apply(sfs.line,1,Calc.Theta.unfolded,w=w.fuli,psi=w.psi)/data.freqs[num.gene,p.len.seq],5)
    Theta.Watt.line <- round(apply(sfs.line,1,Calc.Theta.unfolded,w=w.watt,psi=w.psi)/data.freqs[num.gene,p.len.seq],5)
    Theta.Taji.line <- round(apply(sfs.line,1,Calc.Theta.unfolded,w=w.taj ,psi=w.psi)/data.freqs[num.gene,p.len.seq],5)
    Theta.FayW.line <- round(apply(sfs.line,1,Calc.Theta.unfolded,w=w.fayw,psi=w.psi)/data.freqs[num.gene,p.len.seq],5)
    
    data.results.Fuli <- rbind(data.results.Fuli,data.frame(num.gene,nsam,data.freqs[num.gene,p.len.seq],t(Theta.FuLi.line)))#,(data.freqs[num.gene,c(8)]),data.freqs[num.gene,p.len.seq],t(Theta.FuLi.line)))
    data.results.Watt <- rbind(data.results.Watt,data.frame(num.gene,nsam,data.freqs[num.gene,p.len.seq],t(Theta.Watt.line)))#,(data.freqs[num.gene,c(8)]),data.freqs[num.gene,p.len.seq],t(Theta.Watt.line)))
    data.results.Taj  <- rbind(data.results.Taj, data.frame(num.gene,nsam,data.freqs[num.gene,p.len.seq],t(Theta.Taji.line)))#,(data.freqs[num.gene,c(8)]),data.freqs[num.gene,p.len.seq],t(Theta.Taji.line)))
    data.results.fw   <- rbind(data.results.fw,  data.frame(num.gene,nsam,data.freqs[num.gene,p.len.seq],t(Theta.FayW.line)))#,(data.freqs[num.gene,c(8)]),data.freqs[num.gene,p.len.seq],t(Theta.FayW.line)))
  }
  colnames(data.results.Fuli) <- columns.fuli
  colnames(data.results.Watt) <- columns.watt
  colnames(data.results.Taj)  <- columns.taj
  colnames(data.results.fw)   <- columns.fw
  
  data.results.Fuli <- data.results.Fuli[!duplicated(data.results.Fuli),]
  data.results.Watt <- data.results.Watt[!duplicated(data.results.Watt),]
  data.results.Taj  <- data.results.Taj [!duplicated(data.results.Taj) ,]
  data.results.fw   <- data.results.fw  [!duplicated(data.results.fw)  ,]
  
  data.results.Fuli <- as.matrix(data.results.Fuli)
  data.results.Watt <- as.matrix(data.results.Watt)
  data.results.Taj  <- as.matrix(data.results.Taj)
  data.results.fw   <- as.matrix(data.results.fw)
  
  write.table(file=sprintf("%s/Results_RDiversity_%s_RThetaFuli_unfold_popsize%.0f.txt",path_file,data.names[fil],popsize),x=data.results.Fuli,quote=F,row.names=F)
  write.table(file=sprintf("%s/Results_RDiversity_%s_RThetaWatt_unfold_popsize%.0f.txt",path_file,data.names[fil],popsize),x=data.results.Watt,quote=F,row.names=F)
  write.table(file=sprintf("%s/Results_RDiversity_%s_RThetaTaji_unfold_popsize%.0f.txt",path_file,data.names[fil],popsize),x=data.results.Taj, quote=F,row.names=F)
  write.table(file=sprintf("%s/Results_RDiversity_%s_RThetaFayW_unfold_popsize%.0f.txt",path_file,data.names[fil],popsize),x=data.results.fw , quote=F,row.names=F)
  
  #ESTIMATE REAL NUCLEOTIDE DIVERSITY PER POP (ind) (that is, NO RELATIVE TO TOTAL, option -o 1)
  #define data frames for sending results
  columns.A  <- c("NAME","nsam","Length","ThetaAFuLi","ThetaAWatt","ThetaATaj","ThetaAFayWu")#sprintf("ThetaATaj_%0.f",c(1:npops)))
  data.results.A  <- data.frame(matrix(nrow=0,ncol=length(columns.A)))
  colnames(data.results.A) <- columns.A
  
  data.freqs <- read.table(file=sprintf("%s/%s.columns.txt",path_file,data.names1[fil]),header=T)
  
  #for each gene extract the sfs per lineage and calculate their theta statistics
  init.pi.values <- grep("Theta.Taj.",colnames(data.freqs))
  init.wat.values <- grep("Theta.Wat.",colnames(data.freqs))
  init.fl.values <- grep("Theta.Fu.Li.",colnames(data.freqs))
  init.fw.values <- grep("Theta.Fay.Wu.",colnames(data.freqs))
  
  p.len.seq <- grep("Eff_length2_pop",colnames(data.freqs))[1]#take only one. All are equal under option -u 0
  num.gene <- 1
  for(num.gene in 1:dim(data.freqs)[1]) {
    Theta.A.line <- c(
      data.freqs[num.gene,init.fl.values],
      data.freqs[num.gene,init.wat.values],
      data.freqs[num.gene,init.pi.values],
      data.freqs[num.gene,init.fw.values]
    )/data.freqs[num.gene,p.len.seq]
    data.results.A  <- rbind(data.results.A,data.frame(num.gene,nsam,data.freqs[num.gene,p.len.seq],Theta.A.line[1],Theta.A.line[2],Theta.A.line[3],Theta.A.line[4])) #,data.freqs[num.gene,p.len.seq],Theta.A.line[1],Theta.A.line[2],Theta.A.line[3],Theta.A.line[4]))
  }
  colnames(data.results.A)  <- columns.A
  data.results.A  <- data.results.A[!duplicated(data.results.A),]
  write.table(file=sprintf("%s/Results_ADiversity_unfold.txt",path_file),x=data.results.A, quote=F,row.names=F)
  
  # Difference between relative and absolute pi within ind
  write.table(file=sprintf("%s/Results_proportions_%s_unfold_popsize%.0f.txt",path_file,data.names[fil],popsize),sprintf("npops=%.0f. Times ThetaFuLi is larger than ThetaRel.FuLi: %.3f",npops,mean(data.results.A[,4])/mean(apply(data.results.Fuli[,-c(1:3)],1,mean))),quote=F,row.names=F, append = F)
  write.table(file=sprintf("%s/Results_proportions_%s_unfold_popsize%.0f.txt",path_file,data.names[fil],popsize),sprintf("npops=%.0f. Times ThetaWat is larger than ThetaRel.Watt: %.3f" ,npops,mean(data.results.A[,5])/mean(apply(data.results.Watt[,-c(1:3)],1,mean))),quote=F,row.names=F,append = T)
  write.table(file=sprintf("%s/Results_proportions_%s_unfold_popsize%.0f.txt",path_file,data.names[fil],popsize),sprintf("npops=%.0f. Times ThetaTaj is larger than ThetaRel.Taj: %.3f"  ,npops,mean(data.results.A[,6])/mean(apply(data.results.Taj[,-c(1:3)],1,mean))),quote=F,row.names=F,append = T)
  write.table(file=sprintf("%s/Results_proportions_%s_unfold_popsize%.0f.txt",path_file,data.names[fil],popsize),sprintf("npops=%.0f. Times ThetaFayW is larger than ThetaRel.Fayw: %.3f",npops,mean(data.results.A[,7])/mean(apply(data.results.fw[,-c(1:3)],1,mean))),quote=F,row.names=F,append = T)
  
  par(mfrow=c(2,1))
  max.val <- max(data.results.A[,-c(1:3)],data.results.Watt[,-c(1:3)],data.results.Taj[,-c(1:3)],data.results.fw[,-c(1:3)],na.rm=T)
  plot(density(data.results.A$ThetaAWatt),main=sprintf("Model %s \n Theta from SFS total",model[fil]),xlim=c(0,max.val),ylim=c(0,25000)); abline(v=mean(data.results.A$ThetaAWatt,na.rm=T))
  lines(density(data.results.A$ThetaAFuLi),col="blue"); abline(v=mean(data.results.A$ThetaAFuLi,na.rm=T),col="blue")
  lines(density(data.results.A$ThetaATaj),col="red"); abline(v=mean(data.results.A$ThetaATaj,na.rm=T),col="red")
  lines(density(data.results.A$ThetaAFayWu),col="green"); abline(v=mean(data.results.A$ThetaAFayWu,na.rm=T),col="green")
  legend("topright",legend=c("Watt","FuLi","Taj","FayWu"),lty = 1,col=c("black","blue","red","green"))
  
  plot(density(data.results.Watt[,-c(1:3)]),main=sprintf("Model %s \n Theta from relative SFS",model[fil]),xlim=c(0,max.val),ylim=c(0,25000)); abline(v=mean(data.results.Watt[,-c(1:4)],na.rm=T))
  lines(density(data.results.Fuli[,-c(1:3)]),col="blue"); abline(v=mean(data.results.Fuli[,-c(1:4)],na.rm=T),col="blue")
  lines(density(data.results.Taj[,-c(1:3)]),col="red"); abline(v=mean(data.results.Taj[,-c(1:4)],na.rm=T),col="red")
  lines(density(data.results.fw[,-c(1:3)]),col="green"); abline(v=mean(data.results.fw[,-c(1:4)],na.rm=T),col="green")
  legend("topright",legend=c("Watt","FuLi","Taj","FayWu"),lty = 1,col=c("black","blue","red","green"))
  
  #Plot the distribution of r_theta values per iteration
  par(mfrow=c(2,2))
  plot(density(data.results.Fuli[1,c(4:(4+nsam/popsize-1))]),
       xlim=c(0,1e-3),ylim=c(0,1e5),main=sprintf("%s\nrFu&Li ind. Distribution per iteration",model[fil]));
  for(i in 2:length(data.results.Fuli[,1])) {
    lines(density(data.results.Fuli[i,c(4:(4+nsam/popsize-1))]));
  }
  plot(density(data.results.Watt[1,c(4:(4+nsam/popsize-1))]),
       xlim=c(0,1e-3),ylim=c(0,1e5),main=sprintf("%s\nrWatterson ind. Distribution per iteration",model[fil]));
  for(i in 2:length(data.results.Watt[,1])) {
    lines(density(data.results.Watt[i,c(4:(4+nsam/popsize-1))]));
  }
  plot(density(data.results.Taj[1,c(4:(4+nsam/popsize-1))]),
       xlim=c(0,1e-3),ylim=c(0,1e5),main=sprintf("%s\nrTajima ind. Distribution per iteration",model[fil]));
  for(i in 2:length(data.results.Taj[,1])) {
    lines(density(data.results.Taj[i,c(4:(4+nsam/popsize-1))]));
  }
  plot(density(data.results.fw[1,c(4:(4+nsam/popsize-1))]),
       xlim=c(0,1e-3),ylim=c(0,1e5),main=sprintf("%s\nrFay&Wu ind. Distribution per iteration",model[fil]));
  for(i in 2:length(data.results.Taj[,1])) {
    lines(density(data.results.fw[i,c(4:(4+nsam/popsize-1))]));
  }
  
  ###############################################################################
  ## calculate the standard deviation for fil=1 of different tests of neutrality:
  
  #Comparing SFS (validation)
  # data.results.A$ThetaATaj -  data.results.A$ThetaAWatt for each row
  # data.results.A$ThetaATaj -  data.results.A$ThetaAFuLi for each row  
  # data.results.A$ThetaATaj -  data.results.A$ThetaAFayWu for each row  
  # data.results.A$ThetaAWatt -  data.results.A$ThetaAFayWu for each row  
  
   #Comparing SFS vs rSFS
  # data.results.A$ThetaAWatt -  data.results.Watt[,-c(1:3)] for each row and indiviual 
  # data.results.A$ThetaAFuLi -  data.results.Fuli[,-c(1:3)] for each row and indiviual 
  # data.results.A$ThetaATaj  -  data.results.Taj[,-c(1:3)] for each row and indiviual 
  # data.results.A$ThetaAFayWu - data.results.fw[,-c(1:3)] for each row and indiviual 
  
  #or in relation to a single SFS statistic (theta Watt)
  #SAME BEFORE: data.results.A$ThetaAWatt -  data.results.Watt[,-c(1:3)] for each row and indiviual 
  # data.results.A$ThetaAWatt -  data.results.Fuli[,-c(1:3)] for each row and indiviual 
  # data.results.A$ThetaAWatt -  data.results.Taj[,-c(1:3)] for each row and indiviual 
  # data.results.A$ThetaAWatt -  data.results.fw[,-c(1:3)] for each row and indiviual 
  
  #or direct comparison rSFS from different estimators
  # data.results.Taj[,-c(1:3)] - data.results.Watt[,-c(1:3)]
  # data.results.Taj[,-c(1:3)] - data.results.Fuli[,-c(1:3)]
  # data.results.Taj[,-c(1:3)] - data.results.fw[,-c(1:3)]
  # data.results.Watt[,-c(1:3)] - data.results.fw[,-c(1:3)]
  
  if(fil==1) {
    #Comparing SFS (validation)
    # data.results.A$ThetaATaj -  data.results.A$ThetaAWatt for each row
    sd.ATaj.AWatt <- sd(data.results.A$ThetaATaj - data.results.A$ThetaAWatt)
    # data.results.A$ThetaATaj -  data.results.A$ThetaAFuLi for each row  
    sd.ATaj.AFuLi <- sd(data.results.A$ThetaATaj - data.results.A$ThetaAFuLi)
    # data.results.A$ThetaATaj -  data.results.A$ThetaAFayWu for each row  
    sd.ATaj.AFayWu <- sd(data.results.A$ThetaATaj - data.results.A$ThetaAFayWu)
    # data.results.A$ThetaAWatt -  data.results.A$ThetaAFayWu for each row  
    sd.AWatt.AFayWu <- sd(data.results.A$ThetaAWatt - data.results.A$ThetaAFayWu)

    #Comparing SFS vs rSFS
    #data.results.A$ThetaAWatt -  data.results.Watt[,-c(1:3)])
    sd.AWatt.RWatt <- 0
    for(i in c(1:dim(data.freqs)[1])) {
      sd.AWatt.RWatt <- sd.AWatt.RWatt + sd(data.results.A$ThetaAWatt[i] - data.results.Watt[i,-c(1:3)])
    }
    sd.AWatt.RWatt <- sd.AWatt.RWatt / dim(data.freqs)[1]
    # data.results.A$ThetaAFuLi -  data.results.Fuli[,-c(1:3)] for each row and indiviual 
    sd.AFuLi.RFuLi <- 0
    for(i in c(1:dim(data.freqs)[1])) {
      sd.AFuLi.RFuLi <- sd.AFuLi.RFuLi + sd(data.results.A$ThetaAFuLi[i] - data.results.Fuli[i,-c(1:3)])
    }
    sd.AFuLi.RFuLi <- sd.AFuLi.RFuLi / dim(data.freqs)[1]
    # data.results.A$ThetaATaj  -  data.results.Taj[,-c(1:3)] for each row and indiviual 
    sd.ATaj.RTaj <- 0
    for(i in c(1:dim(data.freqs)[1])) {
      sd.ATaj.RTaj <- sd.ATaj.RTaj + sd(data.results.A$ThetaATaj[i] - data.results.Taj[i,-c(1:3)])
    }
    sd.ATaj.RTaj <- sd.ATaj.RTaj / dim(data.freqs)[1]
    # data.results.A$ThetaAFayWu - data.results.fw[,-c(1:3)] for each row and indiviual 
    sd.AFW.RFW <- 0
    for(i in c(1:dim(data.freqs)[1])) {
      sd.AFW.RFW <- sd.AFW.RFW + sd(data.results.A$ThetaAFayWu[i] - data.results.fw[i,-c(1:3)])
    }
    sd.AFW.RFW <- sd.AFW.RFW / dim(data.freqs)[1]
    
    #or in relation to a single SFS statistic (theta Watt)
    # data.results.A$ThetaAWatt -  data.results.Watt[,-c(1:3)] for each row and indiviual 
    # data.results.A$ThetaAWatt -  data.results.Fuli[,-c(1:3)] for each row and indiviual 
    sd.AWatt.RFuLi <- 0
    for(i in c(1:dim(data.freqs)[1])) {
      sd.AWatt.RFuLi <- sd.AWatt.RFuLi + sd(data.results.A$ThetaAWatt[i] - data.results.Fuli[i,-c(1:3)])
    }
    sd.AWatt.RFuLi <- sd.AWatt.RFuLi / dim(data.freqs)[1]
    # data.results.A$ThetaAWatt -  data.results.Taj[,-c(1:3)] for each row and indiviual 
    sd.AWatt.RTaj <- 0
    for(i in c(1:dim(data.freqs)[1])) {
      sd.AWatt.RTaj <- sd.AWatt.RTaj + sd(data.results.A$ThetaAWatt[i] - data.results.Taj[i,-c(1:3)])
    }
    sd.AWatt.RTaj <- sd.AWatt.RTaj / dim(data.freqs)[1]
    # data.results.A$ThetaAWatt -  data.results.fw[,-c(1:3)] for each row and indiviual 
    sd.AWatt.RFW <- 0
    for(i in c(1:dim(data.freqs)[1])) {
      sd.AWatt.RFW <- sd.AWatt.RFW + sd(data.results.A$ThetaAWatt[i] - data.results.fw[i,-c(1:3)])
    }
    sd.AWatt.RFW <- sd.AWatt.RFW / dim(data.freqs)[1]
    
    #or direct comparison rSFS from different estimators
    # data.results.Taj[,-c(1:3)] - data.results.Watt[,-c(1:3)]
    sd.RTaj.RWatt <- 0
    for(i in c(1:dim(data.freqs)[1])) {
      sd.RTaj.RWatt <- sd.RTaj.RWatt + sd(data.results.Taj[i,-c(1:3)] - data.results.Watt[i,-c(1:3)])
    }
    sd.RTaj.RWatt <- sd.RTaj.RWatt / dim(data.freqs)[1]
    # data.results.Taj[,-c(1:3)] - data.results.Fuli[,-c(1:3)]
    sd.RTaj.RFuLi <- 0
    for(i in c(1:dim(data.freqs)[1])) {
      sd.RTaj.RFuLi <- sd.RTaj.RFuLi + sd(data.results.Taj[i,-c(1:3)] - data.results.Fuli[i,-c(1:3)])
    }
    sd.RTaj.RFuLi <- sd.RTaj.RFuLi / dim(data.freqs)[1]
    # data.results.Taj[,-c(1:3)] - data.results.fw[,-c(1:3)]
    sd.RTaj.RFW <- 0
    for(i in c(1:dim(data.freqs)[1])) {
      sd.RTaj.RFW <- sd.RTaj.RFW + sd(data.results.Taj[i,-c(1:3)] - data.results.fw[i,-c(1:3)])
    }
    sd.RTaj.RFW <- sd.RTaj.RFW / dim(data.freqs)[1]
    
    # data.results.Watt[,-c(1:3)] - data.results.fw[,-c(1:3)]
    sd.RWatt.RFW <- 0
    for(i in c(1:dim(data.freqs)[1])) {
      sd.RWatt.RFW <- sd.RWatt.RFW + sd(data.results.Watt[i,-c(1:3)] - data.results.fw[i,-c(1:3)])
    }
    sd.RWatt.RFW <- sd.RWatt.RFW / dim(data.freqs)[1]
  }
  
  ######################################################################
  ## MATRIX WITH ALL RESULTS ##
  
  name.test <- c("ATaj.AWatt","ATaj.AFuLi","ATaj.AFaywu", "AWatt.AFayWu",
                 "AWatt.RWatt","AFuLi.RFuLi","ATaj.RTaj", "AFayWu.RFayWu",
                 "AWatt.RFuLi","AWatt.RTaj","AWatt.RFayWu",
                 "RTaj.RWatt","RTaj.RFuLi","RTaj.RFayWu","RWatt.RFayWu")
  test.matrix <- array(0,dim=c(dim(data.freqs)[1]*nind,15))
  for(i in c(1:dim(data.freqs)[1])) {                          
    test.matrix[((i-1)*nind+1):(i*nind),1] <-  (data.results.A$ThetaATaj[i] - data.results.A$ThetaAWatt[i])/sd.ATaj.AWatt
    test.matrix[((i-1)*nind+1):(i*nind),2] <-  (data.results.A$ThetaATaj[i] - data.results.A$ThetaAFuLi[i])/sd.ATaj.AFuLi
    test.matrix[((i-1)*nind+1):(i*nind),3] <-  (data.results.A$ThetaATaj[i] - data.results.A$ThetaAFayWu[i])/sd.ATaj.AFayWu
    test.matrix[((i-1)*nind+1):(i*nind),4] <-  (data.results.A$ThetaAWatt[i] - data.results.A$ThetaAFayWu[i])/sd.AWatt.AFayWu

    test.matrix[((i-1)*nind+1):(i*nind),5] <- ((data.results.A$ThetaAWatt[i] - data.results.Watt[i,-c(1:3)])/sd.AWatt.RWatt)
    test.matrix[((i-1)*nind+1):(i*nind),6] <- ((data.results.A$ThetaAFuLi[i] - data.results.Fuli[i,-c(1:3)])/sd.AFuLi.RFuLi)
    test.matrix[((i-1)*nind+1):(i*nind),7] <- ((data.results.A$ThetaATaj[i]  - data.results.Taj[i,-c(1:3)]) /sd.ATaj.RTaj)
    test.matrix[((i-1)*nind+1):(i*nind),8] <- ((data.results.A$ThetaAFayWu[i] - data.results.fw[i,-c(1:3)]) /sd.AFW.RFW)
    
    test.matrix[((i-1)*nind+1):(i*nind),9] <- ((data.results.A$ThetaAWatt[i] - data.results.Fuli[i,-c(1:3)])/sd.AWatt.RFuLi)
    test.matrix[((i-1)*nind+1):(i*nind),10] <- ((data.results.A$ThetaAWatt[i] - data.results.Taj[i,-c(1:3)]) /sd.AWatt.RTaj)
    test.matrix[((i-1)*nind+1):(i*nind),11] <- ((data.results.A$ThetaAWatt[i] - data.results.fw[i,-c(1:3)])  /sd.AWatt.RFW)
    
    test.matrix[((i-1)*nind+1):(i*nind),12]  <- ((data.results.Taj[i,-c(1:3)] - data.results.Watt[i,-c(1:3)]) /sd.RTaj.RWatt)
    test.matrix[((i-1)*nind+1):(i*nind),13]  <- ((data.results.Taj[i,-c(1:3)] - data.results.Fuli[i,-c(1:3)]) /sd.RTaj.RFuLi)
    test.matrix[((i-1)*nind+1):(i*nind),14] <- ((data.results.Taj[i,-c(1:3)] - data.results.fw[i,-c(1:3)])   /sd.RTaj.RFW)
    test.matrix[((i-1)*nind+1):(i*nind),15] <- ((data.results.Watt[i,-c(1:3)] - data.results.fw[i,-c(1:3)])  /sd.RWatt.RFW)
  }
  
  coord.test.matrix <- test.matrix
  colnames(coord.test.matrix) <- name.test
  if(fil==1) coord.test.matrix1 <- coord.test.matrix
  
  par(mfrow=c(2,2))
  for(i in c(1:8,5,9:15)) { # one (5) is repeated for coherence in grouping in plots types of tests
    plot(density(coord.test.matrix[,i]),main=sprintf("Model %s \n Test sd %s",model[fil],name.test[i]))
    abline(v=mean(coord.test.matrix[,i]))
    if(fil > 1) {
      lines(density(coord.test.matrix1[,i]),col="grey")
      abline(v=mean(coord.test.matrix1[,i]),col="grey")
      #legend("topright",legend=c("SNM"),lty=1,col="grey")
    }
  }
}
######################################################################
dev.off()












to.do <- function()
{
  ######################################################################
  ######################################################################
  #
  # FOLDED
  #
  #ESTIMATE RELATIVE NUCLEOTIDE DIVERSITY PER POP (ind) (that is, IN RELATION TO TOTAL, option -o 92)
  #define data frames for sending results
  columns.fuli <- c("NAME","nsam","scaffold","Start","End","Length",sprintf("ThetaRFuLi_%0.f",c(1:(npops))))
  data.results.Fuli <- data.frame(matrix(nrow=0,ncol=length(columns.fuli)))
  colnames(data.results.Fuli) <- columns.fuli
  
  columns.watt <- c("NAME","nsam","scaffold","Start","End","Length",sprintf("ThetaRWatt_%0.f",c(1:(npops))))
  data.results.Watt <- data.frame(matrix(nrow=0,ncol=length(columns.watt)))
  colnames(data.results.Watt) <- columns.watt
  
  columns.taj  <- c("NAME","nsam","scaffold","Start","End","Length",sprintf("ThetaRTaj_%0.f",c(1:(npops))))
  data.results.Taj  <- data.frame(matrix(nrow=0,ncol=length(columns.taj)))
  colnames(data.results.Taj) <- columns.taj
  
  #data rSFS
  data.freqs <- read.table(file=sprintf("%s/SNM.Results.mstatspop-o92_%.0f_fold.txt",path_file,npops))
  
  w.fuli <- weight.fuli.folded(nsam)
  w.watt <- weight.watt.folded(nsam)
  w.taj  <- weight.taj.folded(nsam)
  
  w.phi <- phi.i(nsam) #w.watt
  w.psi <- psi.fold.ij(nsam,popsize)
  
  #for each gene/iteration extract the sfs per lineage and calculate their theta statistics
  init.mat <- grep("rSFS",data.freqs[1,])
  for(num.gene in 1:dim(data.freqs)[1]) {
    p.len.seq <- grep("Length:",data.freqs[num.gene,]) + 1
    sfs.line <- matrix(as.numeric(data.freqs[num.gene,((init.mat+1):(init.mat+(npops*floor(nsam/2))))]),byrow=T,nrow=npops,ncol=(floor(nsam/2)))
    
    Theta.FuLi.line <- round(apply(sfs.line,1,Calc.Theta.folded,w=w.fuli,phi=w.phi,psi=w.psi)/data.freqs[num.gene,p.len.seq],5)
    Theta.Watt.line <- round(apply(sfs.line,1,Calc.Theta.folded,w=w.watt,phi=w.phi,psi=w.psi)/data.freqs[num.gene,p.len.seq],5)
    Theta.Taji.line <- round(apply(sfs.line,1,Calc.Theta.folded,w=w.taj ,phi=w.phi,psi=w.psi)/data.freqs[num.gene,p.len.seq],5)
    
    data.results.Fuli <- rbind(data.results.Fuli,data.frame(num.gene,nsam,(data.freqs[num.gene,c(4)]),(data.freqs[num.gene,c(6)]),(data.freqs[num.gene,c(8)]),data.freqs[num.gene,p.len.seq],t(Theta.FuLi.line)))
    data.results.Watt <- rbind(data.results.Watt,data.frame(num.gene,nsam,(data.freqs[num.gene,c(4)]),(data.freqs[num.gene,c(6)]),(data.freqs[num.gene,c(8)]),data.freqs[num.gene,p.len.seq],t(Theta.Watt.line)))
    data.results.Taj  <- rbind(data.results.Taj, data.frame(num.gene,nsam,(data.freqs[num.gene,c(4)]),(data.freqs[num.gene,c(6)]),(data.freqs[num.gene,c(8)]),data.freqs[num.gene,p.len.seq],t(Theta.Taji.line)))
  }
  colnames(data.results.Fuli) <- columns.fuli
  colnames(data.results.Watt) <- columns.watt
  colnames(data.results.Taj)  <- columns.taj
  
  data.results.Fuli <- data.results.Fuli[!duplicated(data.results.Fuli),]
  data.results.Watt <- data.results.Watt[!duplicated(data.results.Watt),]
  data.results.Taj  <- data.results.Taj [!duplicated(data.results.Taj) ,]
  
  write.table(file=sprintf("%s/Results_RDiversity_%.0f_RThetaFuli_fold.txt",path_file,npops),x=data.results.Fuli,quote=F,row.names=F)
  write.table(file=sprintf("%s/Results_RDiversity_%.0f_RThetaWatt_fold.txt",path_file,npops),x=data.results.Watt,quote=F,row.names=F)
  write.table(file=sprintf("%s/Results_RDiversity_%.0f_RThetaTaji_fold.txt",path_file,npops),x=data.results.Taj, quote=F,row.names=F)
  
  #ESTIMATE REAL NUCLEOTIDE DIVERSITY PER POP (ind) (that is, NO RELATIVE TO TOTAL, option -o 1)
  #define data frames for sending results
  columns.A  <- c("NAME","nsam","scaffold","Start","End","Length","ThetaAFuLi","ThetaAWatt","ThetaATaj")#sprintf("ThetaATaj_%0.f",c(1:npops)))
  data.results.A  <- data.frame(matrix(nrow=0,ncol=length(columns.A)))
  colnames(data.results.A) <- columns.A
  
  data.freqs <- read.table(file=sprintf("%s/SNM.Results.mstatspop-o1_fold.txt_selected_columns.txt",path_file),header=T)
  
  #for each gene extract the sfs per lineage and calculate their theta statistics
  init.pi.values <- grep("Theta.Taj.",colnames(data.freqs))
  init.wat.values <- grep("Theta.Wat.",colnames(data.freqs))
  init.fl.values <- grep("Theta.Fu.Li.",colnames(data.freqs))
  
  p.len.seq <- grep("Eff_length2_pop",colnames(data.freqs))[1]#take only one. All are equal under option -u 0
  for(num.gene in 1:dim(data.freqs)[1]) {
    Theta.A.line <- c(data.freqs[num.gene,init.fl.values],data.freqs[num.gene,init.wat.values],data.freqs[num.gene,init.pi.values])/data.freqs[num.gene,p.len.seq]
    data.results.A  <- rbind(data.results.A, data.frame(num.gene,nsam,(data.freqs[num.gene,c(3)]),(data.freqs[num.gene,c(4)]),(data.freqs[num.gene,c(5)]),data.freqs[num.gene,p.len.seq],Theta.A.line[1],Theta.A.line[2],Theta.A.line[3]))
  }
  colnames(data.results.A)  <- columns.A
  data.results.A  <- data.results.A[!duplicated(data.results.A),]
  write.table(file=sprintf("%s/Results_ADiversity_fold.txt",path_file),x=data.results.A, quote=F,row.names=F)
  
  # Difference between relative and absolute pi within ind
  write.table(file=sprintf("%s/Results_proportions_%.0f_fold.txt",path_file,npops),sprintf("npops=%.0f. Times ThetaFuLi is larger than ThetaRel.FuLi: %.3f",npops,mean(data.results.A[,7])/mean(apply(data.results.Fuli[,-c(1:6)],1,mean))),quote=F,row.names=F, append = F)
  write.table(file=sprintf("%s/Results_proportions_%.0f_fold.txt",path_file,npops),sprintf("npops=%.0f. Times ThetaWat is larger than ThetaRel.Watt: %.3f" ,npops,mean(data.results.A[,8])/mean(apply(data.results.Watt[,-c(1:6)],1,mean))),quote=F,row.names=F,append = T)
  write.table(file=sprintf("%s/Results_proportions_%.0f_fold.txt",path_file,npops),sprintf("npops=%.0f. Times ThetaTaj is larger than ThetaRel.Taj: %.3f"  ,npops,mean(data.results.A[,9])/mean(apply(data.results.Taj[,-c(1:6)],1,mean))),quote=F,row.names=F,append = T)
}
