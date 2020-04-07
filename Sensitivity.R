

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




"Sensitivity"<-function(Alpha,b,Rho,k,i0,
                        est=c("gamma","b","rho"))
{
  Theta<-Theta.fun(Gamma,Rho,b,i0,k)
  Phy0Inv<-1/mean(
    sapply(1:(k-1),function(i)
      Phy0(mu=i/(k+1),Alpha,b,Rho,k,i=i,Theta)))
  
  Ja_123<-sapply(i0,function(i)
    JAlpha(Alpha=Alpha,Gamma=Gamma,
           Rho=Rho,b=b,mu=i/(k+1),est=est))
  Coeff<-(Theta*(Alpha^2+(1+Alpha)^2*exp(-(1+Alpha)))/
            ((Alpha*(1+Alpha)^2)*(k-1)))
  
  Sens<-Coeff*sqrt(t(Ja_123)%*%Phy0Inv%*%Ja_123)
  return(Sens)
}

"k"=c(10,30,50,70,150,200,300,500);
"Alpha"=seq(0.01,1,by=0.02);"b"=1.2;
"Rho"=-0.5; "est"="gamma";"Gamma"=0.5

SEN1<-sapply(Alpha,function(alpha)
 sapply(k,function(t0)
    Sensitivity(Alpha=alpha,b,Rho,k=t0,
                i0=t0/2,est="gamma"))
)

#SEN1<-SEN1[,-1]

plot(y=SEN1[1,],x=Alpha,
     xlab=expression(alpha),
     ylab=expression(s[j]),type="l",
     ylim=c(0,6))
#Ind<-c(2,4,8)
lapply(2:4,
       function(i)lines(x=Alpha,
                        y=SEN1[i,],
                        lty=i))
legend("topright",lty=1:4,
       legend=k[1:4],
       title="k")


###################################
plot(y=SEN1[5,],x=Alpha,
     xlab=expression(alpha),
     ylab=expression(s[j]),type="l",
     lty=5,
     ylim=c(0,0.35))
#Ind<-c(2,4,8)
lapply(6:8,
       function(i)lines(x=Alpha,
                        y=SEN1[i,],
                        lty=i))
legend("topright",lty=5:8,
       legend=k[-c(1:4)],
       title="k")




