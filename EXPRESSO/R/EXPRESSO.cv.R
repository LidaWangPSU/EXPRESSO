#' Title
#'
#' @param r.t 
#' @param r.v 
#' @param R.t 
#' @param R.v 
#' @param n.t 
#' @param n.v 
#' @param lambda_list_1 
#' @param group 
#' @param maxIter 
#' @param phi 
#'
#' @return
#' @export
#'
#' @examples
cv.lassosum<-function(r.t,r.v,R.t,R.v,n.t,n.v,lambda_list_1,group,maxIter,phi){
  init.vec<- array(0, dim = length(r.t))
  res.sol <- array(0, dim = c(length(lambda_list_1), length(r.t)))
  
  s.vec<-rep(1,length(r.t))
  v.vec <- 1/s.vec
  r2_cv<-c()
  bic_cv<-c()
  
  alpha_0<-0.5
  
  
  for(ii in 1:length(lambda_list_1)){
    #print(ii)
    alpha.ori <- c(alpha_0,alpha_0*phi)
    lambda.ori <- c(lambda_list_1[ii],lambda_list_1[ii]*phi)
    alpha <- alpha.ori*lambda.ori
    lambda <- lambda.ori*(1-alpha.ori)/2;
    
    init.vec <- rep(0,length(r.t))
    
    maxIter=1000
    tmp <- fast.lasso.sum.ess(r.t, s.vec, R.t,n.t, group, lambda, alpha, init.vec, maxIter)
    res.cv <- as.vector(tmp$beta)
    res.sol[ii,]<-res.cv
    
    k<-length(which(res.cv!=0))
    
    r2_cv[ii]   <- 1-mse.eval(res.cv, r.v, v.vec, R.v, n.v)
    
    bic_cv[ii]   <- bic.eval(res.cv, r.v, v.vec, R.v, n.v)
    
    #print(c("# coef:",ii,length(which(res.cv!=0))))
  
    
  }
  
  
  return(list(r2_cv,bic_cv))
}



#' Title
#'
#' @param res.lasso 
#' @param u.vec 
#' @param v.vec 
#' @param cov.mat 
#' @param n 
#'
#' @return
#' @export
#'
#' @examples
bic.eval <- function(res.lasso,u.vec,v.vec,cov.mat,n) {
  ix.var <- which(res.lasso!=0);
  if(length(ix.var)==0) return(NA);
  beta.vec <- res.lasso[ix.var];
  u.vec <- u.vec[ix.var];
  cov.mat <- cov.mat[ix.var,ix.var];
  L <- (t(beta.vec)%*%cov.mat%*%beta.vec-2*sum(u.vec*beta.vec)+1)
  if(L<=0) return(NA);
  if(L>=1) return(NA);
  bic <- log(n)*length(ix.var)-n*(1-L);
  return(bic);
}



#' Title
#'
#' @param res.lasso 
#' @param u.vec 
#' @param v.vec 
#' @param cov.mat 
#' @param n 
#'
#' @return
#' @export
#'
#' @examples
mse.eval<- function(res.lasso,u.vec,v.vec,cov.mat,n) {
  ix.var <- which(res.lasso!=0);
  if(length(ix.var)==0) return(NA);
  beta.vec <- res.lasso[ix.var];
  u.vec <- u.vec[ix.var];
  cov.mat <- cov.mat[ix.var,ix.var];
  L <- (t(beta.vec)%*%cov.mat%*%beta.vec-2*sum(u.vec*beta.vec)+1)
  if(L<=0) return(NA);
  if(L>=1) return(1);
  mse <- (L);
  return(mse);
}

#' Title
#'
#' @param res.lasso 
#' @param u.vec 
#' @param v.vec 
#' @param cov.mat 
#' @param n 
#'
#' @return
#' @export
#'
#' @examples
cv.eval<- function(res.lasso,u.vec,v.vec,cov.mat,n) {
  ix.var <- which(res.lasso!=0);
  if(length(ix.var)==0) return(NA);
  beta.vec <- res.lasso[ix.var];
  u.vec <- u.vec[ix.var];
  v.vec<-v.vec[ix.var]
  cov.mat <- cov.mat[ix.var,ix.var];
  
  r2<-sum(beta.vec*u.vec)/(t(beta.vec)%*%cov.mat%*%beta.vec)
  
  return(r2);
}




#' Title
#'
#' @param b.vec 
#' @param cor.analysis 
#' @param ss 
#' @param n 
#' @param group 
#' @param maxIter 
#' @param fold 
#'
#' @return
#' @export
#'
#' @examples
EXPRESSO.cv<-function(b.vec,cor.analysis,ss,n,group,maxIter,fold){
  r<-b.vec
  s.vec<-rep(1,length(r))
  v.vec <- 1/s.vec
  r2.mat<-cor.analysis
  nlambda<-25
  
  no.group<-max(group)
  if (no.group > 1) {
    phi.vec <- c(seq(1/3, 1, by = 1/3))
    phi.vec <- phi.vec[length(phi.vec):1]
  }
  if (no.group == 1)
    phi.vec <- 1
  
  
  
  nt<-ceiling(n*(fold-1)/fold)
  nv<-n-nt
  
  set.seed(123)
  cv_fold_testing = caret::createFolds(c(1:n), k = fold, list = TRUE)
  
  set.seed(123)
  cv_fold_tuning = caret::createFolds(c(1:n), k = fold, returnTrain = TRUE)
  
  u.vec.v.new<-data.frame()
  u.vec.t.new<-data.frame()
  R_v<-list()
  R_t<-list()
  
  for(f in 1:fold){

    temp1<-ss[cv_fold_tuning[[f]],]
    tmp1<-colSums(temp1)/nt
    u.vec.t.new<-rbind(u.vec.t.new,tmp1)
    R_t[[f]]<-cor(temp1)
    temp2<-ss[cv_fold_testing[[f]],]
    R_v[[f]]<-cor(temp2)
    tmp2<-colSums(temp2)/nv
    u.vec.v.new<-rbind(u.vec.v.new,tmp2)
    
  }
  cor_ss<-cor((cbind(r,t(u.vec.t.new),t(u.vec.v.new))))

  r2_mse<-c()
  r2_bic<-c()
  r2_p<-c()
  
  out<-data.frame()
  for(k in 1:length(phi.vec)){
  phi<-phi.vec[k]
  for (f in 1:fold) {
    r2<-data.frame()
    bic<-data.frame()
    lambda_list<-data.frame()
    for(ff in 1:fold){
      if(ff==f){next}
      temp1<-ss[-union(cv_fold_testing[[f]],cv_fold_testing[[ff]]),]
      R.t<-cor(temp1)
      r.t<-colSums(temp1)/(n-2*nv)
      r.v<-as.numeric(u.vec.v.new[ff,])
      R.v<-R_v[[ff]]

      alpha_0<-0.5
      lambda.max<-max(abs(r.t))/((n-2*nv)*alpha_0)
      #lambda.min <- as.numeric(quantile(abs(r.t), probs = c(0.5)))/((n-2*nv)*alpha_0)
      lambda.min <- lambda.max*0.1
      lambda_list_1<- exp(seq(log(lambda.max), log(lambda.min), length.out = nlambda))
      lambda_list<-rbind(lambda_list,lambda_list_1)
      
      res<-cv.lassosum(r.t,r.v,R.t,R.v,(n-2*nv),nv,lambda_list_1,group,maxIter,phi)
      r2<-rbind(r2,res[[1]])
      bic<-rbind(bic,res[[2]])
    }
    r2[is.na(r2)]=0
    bic[is.na(bic)]=1e5
    
    
    
    lbd1<- mean(lambda_list[,which.max(colSums(r2))])*(n-2*nv)
    lbd2<- mean(lambda_list[,which.min(colSums(bic))])*(n-2*nv)
    
    
    #####pseudo variable
    mu1<-rep(0,length(r.t))
    
    lbd_p<-c()
  
    for(tt in 1:5){
      set.seed(tt)
      r_pseudo_ss<-mvrnorm(nt,mu1, R_t[[f]])
      r_pseudo<-colSums(r_pseudo_ss)/nt
      lbd_p[tt]<-max(abs(r_pseudo))
    }
    lbd3<-max(lbd_p)/alpha_0
  

    R.t<-R_t[[f]]
    r.t<-u.vec.t.new[f,]
    r.v<-u.vec.v.new[f,]
    R.v<-R_v[[f]]
    init.vec<- array(0, dim = length(r.t))

    alpha.ori <- c(0.5,0.5*phi)
    lbd1<-c(lbd1,lbd1*phi)
    lbd2<-c(lbd2,lbd2*phi)
    lbd3<-c(lbd3,lbd3*phi)
    
    alpha1 <- alpha.ori*lbd1/(n-nv)
    alpha2 <- alpha.ori*lbd2/(n-nv)
    alpha3 <- alpha.ori*lbd3/(n-nv)
    lambda1 <- lbd1*(1-alpha.ori)/2
    lambda2 <- lbd2*(1-alpha.ori)/2;
    lambda3 <- lbd3*(1-alpha.ori)/2;
    
    
    tmp1 <- fast.lasso.sum.ess(as.numeric(r.t), s.vec, R.t, (n-nv), group, lambda1, alpha1, init.vec, maxIter)
    
    res.cv1 <- as.vector(tmp1$beta)
    
    tmp2 <- fast.lasso.sum.ess(as.numeric(r.t), s.vec, R.t, (n-nv), group, lambda2, alpha2, init.vec, maxIter)
    
    res.cv2 <- as.vector(tmp2$beta)
    
    tmp3 <- fast.lasso.sum.ess(as.numeric(r.t), s.vec, R.t, (n-nv), group, lambda3, alpha3, init.vec, maxIter)
    
    res.cv3 <- as.vector(tmp3$beta)
    
    r2_mse[f] <- 1-mse.eval(res.cv1, r.v, v.vec, R.v, nv)
    r2_bic[f] <- 1-mse.eval(res.cv2, r.v, v.vec, R.v, nv)
    r2_p[f] <- 1-mse.eval(res.cv3, r.v, v.vec, R.v, nv)
    
    #print(c("# mse:","fold",f,length(which(res.cv1!=0))))
    #print(c("# bic:","fold",f,length(which(res.cv2!=0))))
    #print(c("# pseudo:","fold",f,length(which(res.cv3!=0))))

  }
  
  
  r_2_mse<-mean(r2_mse,na.rm=T)
  r_2_bic<-mean(r2_bic,na.rm=T)
  r_2_p<-mean(r2_p,na.rm=T)
  if(is.na(r_2_mse)){
    r_2_mse<-0
  }
  if(is.na(r_2_bic)){
    r_2_bic<-0
  }

  
  if(is.na(r_2_p)){
    r_2_p<-r_2_bic
  }
  out_temp1<-c(phi,r_2_mse,r_2_bic,r_2_p)
  out<-rbind(out,out_temp1)
  
  }
  colnames(out)<-c("phi","mse_r2","bic_r2","pseudo_r2")

return(out)
}
