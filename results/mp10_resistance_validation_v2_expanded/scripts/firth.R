# Jeffreys-prior penalized logistic likelihood; profile inference for MP10 coefficient.
# No external package installation or environment changes.
firth_profile=function(X,y,j){
 stopifnot(qr(X)$rank==ncol(X),length(unique(y))==2)
 loglik=function(b){eta=drop(X%*%b);p=plogis(eta);w=p*(1-p);I=crossprod(X*sqrt(w));de=determinant(I,logarithm=TRUE);if(de$sign<=0||!is.finite(as.numeric(de$modulus)))return(-1e100);sum(y*eta-pmax(eta,0)-log1p(exp(-abs(eta))))+.5*as.numeric(de$modulus)}
 fit=optim(rep(0,ncol(X)),function(b)-loglik(b),method='BFGS',control=list(maxit=2000,reltol=1e-12),hessian=TRUE);stopifnot(fit$convergence==0,all(is.finite(fit$par)),min(eigen(fit$hessian,symmetric=TRUE,only.values=TRUE)$values)>0)
 prof=function(v){f=optim(fit$par[-j],function(a){b=fit$par;b[j]=v;b[-j]=a;-loglik(b)},method='BFGS',control=list(maxit=2000,reltol=1e-11));if(f$convergence!=0)stop('profile optimization failed');-f$value}
 peak=-fit$value;cut=qchisq(.95,1);root=function(v)2*(peak-prof(v))-cut;hat=fit$par[j]
 bracket=function(dir){z=hat+dir;step=1;while(root(z)<0&&step<128){step=step*2;z=hat+dir*step};if(root(z)<0)stop('profile CI unbounded');sort(c(hat,z))}
 ci=c(uniroot(root,bracket(-1),tol=1e-7)$root,uniroot(root,bracket(1),tol=1e-7)$root)
 list(beta=hat,low=ci[1],high=ci[2],P=pchisq(max(0,2*(peak-prof(0))),1,lower.tail=FALSE),converged=TRUE,fitted=plogis(drop(X%*%fit$par)))
}
# Numerical validation only: complete separation in a 2x2 design has known Firth OR=441.
test=firth_profile(cbind(1,c(rep(0,10),rep(1,10))),c(rep(0,10),rep(1,10)),2);stopifnot(abs(exp(test$beta)-441)<.02,test$low<test$beta,test$high>test$beta,is.finite(test$P))
