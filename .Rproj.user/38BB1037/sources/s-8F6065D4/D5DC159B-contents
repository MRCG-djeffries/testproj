#!/usr/bin/env Rscript
run_HPV6_HPC=function(){
  require(data.table)
  require(ICCbin)
  # 1 dec 2020 - this is Eds simulation
  args = commandArgs(trailingOnly=TRUE)
  pathin="~/HPV/"
  # data in
  d=readRDS("~/HPV/runfile2.rds")
  Nd=nrow(d)

  vill_var=1
  comp_var=0.5
  hous_var=0.5
  ICC=0.1
  Nsim=2
  total_sampling = matrix(data = 0, nrow = Nsim*Nd,ncol = 4)
  Pout=matrix(data=0,nrow=Nsim,ncol=Nd)
  pathout=paste0("~/HPV/HPV3/run_",args[1])


  setDTthreads(1)
  for (j in 1 : Nd){
      lower_age =  d[j,"lower_age"]
      upper_age =  d[j,"upper_age"]
      prob_val= d[j,"prob_val"]
      village_percent = d[j,"village_percent"]
      compound_percent = d[j,"compound_percent"]
      house_percent = d[j,"house_percent"]
      subject_percent = d[j,"subject_percent"]
      q=readRDS(paste0(pathin,"females.rds"))
      setkey(q, Village,Compound,Household)      
      q=q[Age>=lower_age & Age<upper_age,]
      vcounts=q[, .N, by=.(Village)] 
      compund_in_village=q[, .N, by=.(Village,Compound)]
      household_level = q[, .N, by=.(Village,Compound,Household)]
      for ( i in 1 : Nsim){
        cat(paste0("Sample # = ",i,"\n"))
        L=simHPV(prob_val,village_percent,compound_percent,house_percent,subject_percent,vill_var,comp_var,hous_var,ICC,vcounts,compund_in_village,household_level,q)
        Pout[i,j]=L$prob_estimate
        total_sampling[i+Nsim*(j-1),]=L$nums
      }
  }
  saveRDS(list(Tot=total_sampling,P=Pout),paste0(pathout,".rds"))
}


simHPV=function(prob_val,village_percent,compound_percent,house_percent,subject_percent,vill_var,comp_var,hous_var,ICC,vcounts,compund_in_village,household_level,q){
  # village_percent is the sampling effort for villages
  subs_per_house=NULL
  vill_var=vill_var/100
  comp_var=comp_var/100
  hous_var=hous_var/100
  prob_val = prob_val/100
  Infected = vector(mode="integer",length = length(q$Village))
  # simulate dataset
  vsamp=vcounts$Village
  NC = 0
  NH = 0
  NS = 0  
  for ( i in 1 :length(vsamp)){
    village_p=max(0.00001,rnorm(1,mean=prob_val,sd=vill_var))
    csamp = compund_in_village$Compound[compund_in_village$Village==vsamp[i]]
    for ( j in 1 :  length(csamp)){
      compound_p=max(0.00001,rnorm(1,mean=village_p,sd=comp_var))
      hsamp=household_level$Household[household_level$Village==vsamp[i] & household_level$Compound==csamp[j]]
      for (k in 1 : length(hsamp)){
        m=q[Village==vsamp[i] & Compound==csamp[j] & Household == hsamp[k],which=TRUE]
        #m=which(q$Village==vsamp[i] & q$Compound==csamp[j] & q$Household == hsamp[k])
        house_p=max(0.00001,rnorm(1,mean=compound_p,sd=hous_var))
        if (length(m)==1){
          Infected[m]=rbinom(1,1,house_p)
        }else{
          Infected[m]=rcbin(prop = house_p, prvar = 0, noc = 1, csize = length(m), csvar = 0, rho = ICC)$y
        }
      }
    }
  }
  
  denom_estimate = 0
  numer_estimate = 0
  # sample poportionally from village
  pvillage = vcounts$N/(sum(vcounts$N))
  nvillage = floor(length(vcounts$Village)*village_percent/100)
  vsamp = sample(vcounts$Village,size=nvillage,replace=FALSE,prob=pvillage)
  for ( i in 1 : length(vsamp)){
    comp=compund_in_village$N[compund_in_village$Village==vsamp[i]]
    comp_name = compund_in_village$Compound[compund_in_village$Village==vsamp[i]]
    pcompound = comp/sum(comp)
    ncompound = ceiling(length(comp_name)*compound_percent/100) # always make sure there is at least 1 compound to sample
    NC=NC + ncompound
    csamp = sample(comp_name,size=ncompound,replace=FALSE,prob=pcompound)
    for ( j in 1 : length(csamp)){
      hous = household_level$N[household_level$Village==vsamp[i] & household_level$Compound==csamp[j]]
      hous_name = household_level$Household[household_level$Village==vsamp[i] & household_level$Compound==csamp[j]]
      phouse = hous/sum(hous)
      nhouse = ceiling(length(hous_name)*house_percent/100) # always make sure there is at least 1 subject to sample
      NH = NH + nhouse
      hsamp = sample(hous_name,size=nhouse,replace=FALSE,prob=phouse) # these are households sample
      for ( k in 1 : length(hsamp)){
        m=q[Village==vsamp[i] & Compound==csamp[j] & Household == hsamp[k],which=TRUE]
        #m=which(q$Village==vsamp[i] & q$Compound==csamp[j] & q$Household == hsamp[k])
        nsubject = ceiling(length(m)*subject_percent/100)
        subs_selected = sample(m,size=nsubject,replace=FALSE)
        NS = NS + length(subs_selected)
        numer_estimate = numer_estimate + sum(Infected[subs_selected])
        denom_estimate = denom_estimate + nsubject
        subs_per_house=c(subs_per_house,nsubject)
      }
    }
  }
  return(L=list(prob_estimate = 100*numer_estimate/denom_estimate,N=denom_estimate,S=subs_per_house,nums=c(nvillage,NC,NH,NS)))
}

sample = function(x, size, replace = F, prob = NULL) {
  if (length(x) == 1) return(x)
  base::sample(x, size = size, replace = replace, prob = prob)
}

run_HPV6_HPC()
