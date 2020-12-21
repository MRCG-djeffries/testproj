interpretHPV=function(){
  library(flextable)
  library(officer)
  library(magrittr)
  # to get repository git config --get remote.origin.url
  # test key
  # again
  #z
  # how to add krsa key
  # ssh-keygen -t rsa -C "djeffries@mrc.gm"
  # use defaults - just hit return
  # check exists ls ~/.ssh/*.pub
  # copy key (sudo apt-get install xclip -y if necessary)
  # xclip -selection clipboard < ~/.ssh/id_rsa.pub
  # go here https://github.com/settings/keys
  # log in and paste key
  # git remote set-url origin git@github.com:MRCG-djeffries/testproj.git
  # git remote show origin
  # changed from windows
  # runs the CV version of prev variabiltiy
  # change 2
  # change 3
  # git push -u origin main
  # now from mac
  # now from ubuntu
  d=readRDS("~/HPV/runfile4.rds")
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
    L=readRDS(file=paste0("~/HPV/HPV5/run_",j,".rds"))
    Pout=rbind(Pout,L$P)
    Nv=matrix(L$Tot[,1],nrow=20,ncol=84)
    Nc=matrix(L$Tot[,2],nrow=20,ncol=84)
    Nh=matrix(L$Tot[,3],nrow=20,ncol=84)
    Ns=matrix(L$Tot[,4],nrow=20,ncol=84)
    NV=rbind(NV,Nv)
    NC=rbind(NC,Nc)
    NH=rbind(NH,Nh)
    NS=rbind(NS,Ns)
    
  }
  # test the CI on RR
  z=1-Pout[,1]/Pout[,5]
  hist(z,500)
  dum=quantile(z,c(0.025,0.975))
  LS=100*dum[1]
  US=100*dum[2]
  rv=1/100
  ru=5/100
  rr=rv/ru
  Nv=3619
  Nu=3619
  L=exp(log(rr)-1.96*sqrt((1-rv)/(Nv*rv) + (1-ru)/(Nu*ru)))
  U=exp(log(rr)+1.96*sqrt((1-rv)/(Nv*rv) + (1-ru)/(Nu*ru)))
  cat(paste0("95% CI on VB=80% is (using simple taylor series formula):",100*(1-U)," to ",100*(1-L),"\n"))
  cat(paste0("95% CI on VB=80% is (using simulated data):",LS," to ",US,"\n"))
  
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
  save_as_html(FF,path="~/HPV/sim4.html")
  #doc =read_docx()%>%body_add_flextable(value = FF,split=TRUE)%>%
  #body_end_section_landscape()  %>% pri2nt(target = "~/HPV/sim.docx" )  
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
