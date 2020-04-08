tlibrary(MASS)
library(actuar);
library(evd)

############################################################
## Estimating rho  according to Fraga Alves
################################################
"rho.fraga.alves"<-function(x)
{
  n<-length(x)
  x<-sort(x,decreasing = TRUE)
  U<-numeric(n)
  for (i in 1:(n-1))   U[i]=i*log(x[i]/x[i+1])
  
  k1=min((n-1),round(n^0.995))
  Uk=U[1:k1]
  
  M<-numeric(3)
  T<-numeric(3)
  
  xk <- x[1:k1]/x[k1+1]
  M[1]<-sum(log(xk))/k1
  M[2]<-sum(log(xk)^2)/k1
  M[3]<-sum(log(xk)^3)/k1
  T[1]=(log(M[1])-1/2*log(M[2]/2))/(1/2*log(M[2]/2)-1/3*log(M[3]/6))
  T[2]=((M[1])^(1/2)-(M[2]/2)^(1/4))/((M[2]/2)^(1/4)-(M[3]/6)^(1/6))
  T[3]=((M[1])^(1)-(M[2]/2)^(1/2))/((M[2]/2)^(1/2)-(M[3]/6)^(1/3))
  r.tau.tilde=-abs((3*(T[1]-1))/(T[1]-3))
  
  list(rho=r.tau.tilde)
}


## Function to compute the ERM_M estimator of the
## tail index
"ERM.MD" <- function(data,alpha,Gam.stat)
{
  # MDPD estimators for the GPD.
  # data: the data vector.
  # alpha: the tuning constant alpha.
  # Gam.stat: starting value for gamma. 
  n<-length(data)
  X <- sort(data)
  rho<-rho.fraga.alves(X)$rho               
  gam1<-sapply( 1:(n-1),
                function(i)(1/i)*sum(log(X[n:(n-i+1)]))-log(X[n-i]))
  Z <-(1:(n-1))*(log(X[n:2])-log(X[(n-1):1]))
  startv<-c(median(gam1,na.rm=T),Gam.stat);
  Hill.theta <- matrix(nrow=n, ncol=2)
  Hill.theta[n,] <- startv
  
  ## Computes the MDPE value for the ERM
  ## k: No of top order statistics
  ## 
  "H.gen"<-function(start,Zj,alpha,k,rho)
  {
    #k<-k-1
    assign("start",start,pos=1)
    assign("rho",rho,pos=1)
    assign("k",k,pos=1)
    assign("alpha",alpha,pos=1) 
    "theta.fun"<-function(start,rho,k)
    {
      ## Computes the mean of the ERM
      Gamma<-start[1];bnk<-start[2]
      return(Gamma+bnk*((1:k)/(k+1))^(-rho))
    }
    
    First<-1/((1+alpha)*(theta.fun(start=start,rho=rho,k=k)^alpha))
    Second<-(1+1/alpha)*(1/(theta.fun(start,rho,k)^alpha))*
      exp(-((alpha*Zj)/(theta.fun(start=start,rho=rho,k=k))))
    Hk_iner<-First-Second
    Hk<-mean(Hk_iner)
    return(Hk)
  }  
  
  
  #All<-  t(sapply((n-1):10,function(k)
  
  for (k in (n-1):10)  
  {
   # J <- ((1:k)/(k+1))
    Zj <- Z[1:k]
    assign("Zj",Zj,pos=1)
    # assign("start",start,pos=1)
    assign("rho",rho,pos=1)
    assign("k",k,pos=1)
    assign("alpha",alpha,pos=1)
    
    opt <- nlminb(start=Hill.theta[k+1,],objective=H.gen,
                  k=k,Zj=Zj,rho=rho,alpha=alpha,
                  lower=c(0,0),upper=c(Inf,Inf))
    estimator<-as.vector(opt$par)
    Hill.theta[k,]<-estimator
    # Hill.theta[k,]
  }
  #))
  return(Hill.theta[-c(1:9),])      
}

### Sample usage of ERM_M estimator 
### in estimating the tail index of a sample
## generated from a Pareto with tail index=0.5
Shape1=2 ## tail index=1/Shape
Shape2=0.5
n<-200
eps<-0.05 ## Percentage of contamination

### Generate data from the Frechet distribution with
## tail index 1/Shape1
## and contaminate it with another sample 
## from Frechet distribution with large
## tail index, 1/Shape2.
U.Data<-rfrechet(n=round(n*(1-eps)),shape=Shape1,scale=1)
C.Data<-rfrechet(n=round(n*eps),shape=Shape2,scale=1)
Data<-c(U.Data,C.Data)


### Robustness parameter
Alpha<-c(0.2,1.0)


### Compute the tail indices
Gamma.ERM<-sapply(Alpha,function(alpha)
                  ERM.MD(data=Data,alpha=alpha,
                         Gam.stat=1)[,1])

### Plot the tail indices
par(mar=c(4,5,1,1))
plot(x=1:nrow(Gamma.ERM),y=Gamma.ERM[,1],type="l",
     xlab="k",ylab=expression(hat(gamma)),
     ylim=c(0,1))
### Draw the a horizontal line of the 
### tail index parameter value
abline(h=1/Shape1,col="red")
lines(x=1:nrow(Gamma.ERM),y=Gamma.ERM[,2],lty=2)
legend("bottomleft",legend=
         as.expression(lapply(Alpha, function(x)
  bquote(alpha==.(x)))), lty = 1:2)
