interpretHPV=function(){
  library(flextable)
  library(officer)
  library(magrittr)
  d=readRDS("~/HPV/runfile3.rds")
  Nd=nrow(d)
  Nsim = 20  
  Nrun = 500
  lower_age=rep(0,Nd)
  upper_age=lower_age
  prob_val=lower_age
  village_percent=lower_age
  compound_percent=lower_age
  house_percent=lower_age
  subject_percent=lower_age
  for (j in 1 : Nd){
    lower_age[j] =  d[j,"lower_age"]
    upper_age[j] =  d[j,"upper_age"]
    prob_val[j]= d[j,"prob_val"]
    village_percent[j] = d[j,"village_percent"]
    compound_percent[j] = d[j,"compound_percent"]
    house_percent[j] = d[j,"house_percent"]
    subject_percent[j] = d[j,"subject_percent"]
  }
  
  Pout=NULL
  NV=NULL
  NC=NULL
  NH=NULL
  NS=NULL

  for (j in 1 : Nrun){
    cat(paste0("samp=",j,"\n"))
    L=readRDS(file=paste0("~/HPV/HPV4/run_",j,".rds"))
    Pout=rbind(Pout,L$P)
    Nv=matrix(L$Tot[,1],nrow=20,ncol=168)
    Nc=matrix(L$Tot[,2],nrow=20,ncol=168)
    Nh=matrix(L$Tot[,3],nrow=20,ncol=168)
    Ns=matrix(L$Tot[,4],nrow=20,ncol=168)
    NV=rbind(NV,Nv)
    NC=rbind(NC,Nc)
    NH=rbind(NH,Nh)
    NS=rbind(NS,Ns)

  }
  meanP=colMeans(Pout)
  dum=colQ(Pout,0.025,0.975)
  lowP=dum$LQ
  uppP=dum$UQ
  meanNV=colMeans(NV)
  meanNC=colMeans(NC)
  meanNH=colMeans(NH)
  meanNS=colMeans(NS)
  nsim = rep(10000,Nd)
  samp=1:Nd
  lowage=lower_age
  uppage=upper_age
  villeffort = village_percent
  compeffort=compound_percent
  houseffort=house_percent
  subjeffort=subject_percent
  df=data.frame(scenario=samp,lower_age=lowage,upper_age=uppage,
                vill_eff=villeffort,comp_eff=compeffort,
                house_eff=house_percent,sub_eff=subjeffort,
                lower_95=lowP,upper_95=uppP,point_P=meanP,actual_P=prob_val,
                numV= meanNV,numC=meanNC,numH=meanNH,numS=meanNS)
  FF=flextable(df)
  FF=colformat_num(FF,j=10,digits=3)
  FF=color(FF,j=c(8,9),color="red")
  FF= set_caption(FF, caption = "Simulation")
  FF= autofit(FF)
  save_as_html(FF,path="~/HPV/sim2.html")
  #doc =read_docx()%>%body_add_flextable(value = FF,split=TRUE)%>%
  #body_end_section_landscape()  %>% print(target = "~/HPV/sim.docx" )  
  # fileout = tempfile(fileext = ".docx")
  # fileout = "~/HPV/sim.docx" 
  # print(doc, target = fileout)
}

colQ = function(X,lo,up){
  L=rep(0,ncol(X))
  U=rep(0,ncol(X))
  for ( i in 1 : ncol(X)){
    L[i]=quantile(X[,i],lo)
    U[i]=quantile(X[,i],up)
  }
  return(list(LQ=L,UQ=U))
}
