

"M1.fun"<-function(mu,Alpha,Gamma,b,Rho)
{
mueta<-  ((1+Alpha)/((Gamma+b*mu^{-Rho}))^{Alpha+2})*
    rbind(c(0,0,0),c(0,0,-log(mu)*mu^{-Rho}),
          c(0,-mu^{-Rho}*log(mu),b*mu^{-Rho}*(log(mu))^2))
  return(mueta)
}

"M2.fun"<-function(mu,Alpha,Gamma,b,Rho)
{
  mueta<-  ((1+Alpha)/((Gamma+b*mu^{-Rho}))^{Alpha+2})*
    rbind(c(1,mu^{-Rho},b*mu^{-Rho}*log(mu)),
          c(mu^{-Rho},mu^{-2*Rho},-b*mu^{-2*Rho}*log(mu)),
          c(-b*mu^{-Rho}*log(mu),-b*mu^{-2*Rho}*log(mu),
            b^2*mu^{-2*Rho}*(log(mu))^2))
  return(mueta)
}


"Phy0"<-function(mu,Alpha,b,Rho,k,i,Theta)
{
  Vec<-((Alpha^4+2*Alpha^3+2*Alpha^2+1)/(1+Alpha)^3)*
    M2.fun(mu=i/(k+1),Alpha=Alpha,
           Gamma=Gamma,b=b,Rho=Rho)-
    ((Alpha^2*Theta)/(1+Alpha)^2)*
    M1.fun(mu=i/(k+1),Alpha=Alpha,
           Gamma=Gamma,b=b,Rho=Rho)
  return(Vec)
}

"JAlpha"<-function(Alpha,Gamma,Rho,b,mu,est=c("gamma","b","Rho"))
{
  if(est=="gamma")
  {Ja<-(1+Alpha)/((Gamma+b*mu^{-Rho}))^{Alpha+2}
  }else if(est=="b"){
  Ja<-((1+Alpha)*mu^{-Rho})/((Gamma+b*mu^{-Rho}))^{Alpha+2} 
  }else{
Ja<-((1+Alpha)*b*mu^{-Rho}*log(mu))/((Gamma+b*mu^{-Rho}))^{Alpha+2} 
  }
}

"Theta.fun"<-function(Gamma,Rho,b,i,k)Gamma+b*(i/(k+1))^{-Rho}


"IFfun"<-function(Alpha,b,Rho,k,i0,est=c("gamma","b","rho"),t_i0)
{
  Theta<-Theta.fun(Gamma,Rho,b,i0,k)
  Phy0Inv<-1/mean(
    sapply(1:(k-1),function(i)
      Phy0(mu=i/(k+1),Alpha,b,Rho,k,i=i,Theta)))
  
  Ja_123<-sapply(i0,function(i)
    JAlpha(Alpha=Alpha,Gamma=Gamma,
           Rho=Rho,b=b,mu=i/(k+1),est=est))
  
  if(!Alpha==0){
  Sq.brac<-((Alpha*Theta)/(1+Alpha)^2)+
            (t_i0-Theta)*exp(-Alpha*t_i0/Theta)
  }else{
    Sq.brac<-t_i0-Theta
  }
  
  IF<-(1/(k-1))*Ja_123*Phy0Inv*Sq.brac
  
   
  return(IF)
}


"k"=100;"Alpha"=c(0,0.3,0.5,1);"b"=2;"Rho"=-1;
i0=70; "est"="gamma";"Gamma"=1
IF1<-sapply(Alpha,function(alpha)
  sapply(1:30,function(t0)
  IFfun(Alpha=alpha,b,Rho,k,i0,est="gamma",
        t_i0=t0)))
par(mar=c(4,4,1,1))
plot(1:nrow(IF1),IF1[,1],type="l",
     xlab=expression(t[0]),ylab="IF",
     ylim=c(0,0.1))
lapply(2:4,function(i)
  lines(1:nrow(IF1),IF1[,i],type="l",
        lty=i))
legend("topright",col=rep(1,4),lty=1:4,
       legend=c(expression(paste(alpha, " = ", 0)),
                expression(paste(alpha, " = ", 0.3)),
                expression(paste(alpha, " = ", 0.5)),
                expression(paste(alpha, " = ", 1))))
