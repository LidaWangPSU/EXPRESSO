
#' Title
#'
#' @param b.vec 
#' @param cor.analysis 
#' @param ss 
#' @param n 
#' @param group 
#' @param maxIter 
#' @param fold 
#' @param phi 
#'
#' @return
#' @export
#'
#' @examples
EXPRESSO.BIC<-function(b.vec,cor.analysis,ss,n,group,maxIter,fold,phi){
  nt<-ceiling(n*(fold-1)/fold)
  nv<-n-nt
  r<-b.vec
  s.vec<-rep(1,length(r))
  v.vec <- 1/s.vec
  r2.mat<-cor.analysis
  
  
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
  
  lambda_list<-data.frame()
  bic<-data.frame()

  for (f in 1:fold) {

  
  temp1<-ss[-c(cv_fold_testing[[f]]),]
  R.t<-cor(temp1)
  r.t<-colSums(temp1)/(n-nv)
  r.v<-u.vec.v.new[f,]
  R.v<-R_v[[f]]
  #beta.full<-beta_full[f,]
  init.vec<- array(0, dim = length(r.t))
  
  alpha_0<-0.5
  lambda.max<-max(abs(r))/((n-nv)*alpha_0)
  lambda.min <- lambda.max*0.1
  lambda_list_1<- exp(seq(log(lambda.max), log(lambda.min), length.out = 25))
  lambda_list<-rbind(lambda_list,lambda_list_1)
  
  bic_cv<-c() 
  for(ii in 1:length(lambda_list_1)){
    print(ii)
    alpha_0<-0.5
    alpha.ori <- c(alpha_0,alpha_0*phi)
    lambda.ori <- c(lambda_list_1[ii],lambda_list_1[ii]*phi)
    alpha <- alpha.ori*lambda.ori
    lambda <- lambda.ori*(1-alpha.ori)/2;
    
    init.vec <- rep(0,length(r.t))
    
    maxIter=1000
    tmp <- fast.lasso.sum.ess(r.t, s.vec, R.t,nt, group, lambda, alpha, init.vec, maxIter)
    res.cv <- as.vector(tmp$beta)
    
    print(c("# coef:",ii,length(which(res.cv!=0))))
    
    bic_cv[ii]   <- bic.eval(res.cv, r.v, v.vec, R.v, nv)
    
  }
  
  bic<-rbind(bic,bic_cv)
  }
  
  #print(which.min(colSums(bic)))
  

  lbd<- mean(lambda_list[,which.min(colSums(bic))])*(n-nv)/n
  if(is.na(lbd)){
  lbd<- mean(lambda_list[,10])*(n-nv)/n
  }
  
  lbd<-c(lbd,lbd*phi)
  alpha_0<-0.5
  alpha.ori <- c(alpha_0,alpha_0*phi)
  lambda.ori <- c(lambda_list_1[ii],lambda_list_1[ii]*phi)
  alpha <- alpha.ori*lambda.ori
  lambda <- lambda.ori*(1-alpha.ori)/2;
  
  
  R.all<-cor(ss)
  
  tmp1 <- fast.lasso.sum.ess(b.vec, s.vec, r2.mat, n, group, lambda, alpha, init.vec, maxIter)
  
  res.cv1 <- as.vector(tmp1$beta)
  
  listM<-list(beta=tmp1$beta,l1_penalty=alpha[1],l2_penalty=lambda[1])
  
  
  
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
#' @param phi 
#'
#' @return
#' @export
#'
#' @examples
EXPRESSO.MSE<-function(b.vec,cor.analysis,ss,n,group,maxIter,fold,phi){
  nt<-ceiling(n*(fold-1)/fold)
  nv<-n-nt
  r<-b.vec
  s.vec<-rep(1,length(r))
  v.vec <- 1/s.vec
  r2.mat<-cor.analysis
  
  
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
    #R_t[[f]]<-diag(1/tmp1)%*%cov(temp1)%*%diag(1/tmp1)/n
    temp2<-ss[cv_fold_testing[[f]],]
    R_v[[f]]<-cor(temp2)
    tmp2<-colSums(temp2)/nv
    u.vec.v.new<-rbind(u.vec.v.new,tmp2)
    
  }
  
  lambda_list<-data.frame()
  mse<-data.frame()
  
  for (f in 1:fold) {
    
    
    temp1<-ss[-c(cv_fold_testing[[f]]),]
    R.t<-cor(temp1)
    r.t<-colSums(temp1)/(n-nv)
    r.v<-u.vec.v.new[f,]
    R.v<-R_v[[f]]
    #beta.full<-beta_full[f,]
    init.vec<- array(0, dim = length(r.t))
    
    alpha_0<-0.5
    lambda.max<-max(abs(r))/((n-nv)*alpha_0)

    lambda.min <- lambda.max*0.1
    
    lambda_list_1<- exp(seq(log(lambda.max), log(lambda.min), length.out = 25))
    lambda_list<-rbind(lambda_list,lambda_list_1)
    
    mse_cv<-c() 
    for(ii in 1:length(lambda_list_1)){
      print(ii)
      alpha_0<-0.5
      alpha.ori <- c(alpha_0,alpha_0*phi)
      lambda.ori <- c(lambda_list_1[ii],lambda_list_1[ii]*phi)
      alpha <- alpha.ori*lambda.ori
      lambda <- lambda.ori*(1-alpha.ori)/2;
      
      init.vec <- rep(0,length(r.t))
      
      maxIter=1000
      tmp <- fast.lasso.sum.ess(r.t, s.vec, R.t,nt, group, lambda, alpha, init.vec, maxIter)
      res.cv <- as.vector(tmp$beta)
      print(c("# coef:",ii,length(which(res.cv!=0))))
      
      mse_cv[ii]   <- 1-mse.eval(res.cv, r.v, v.vec, R.v, nv)
      #mse_cv[ii]   <- cv.eval(res.cv, r.v, v.vec, R.v, nv)
      
    }
    
    mse<-rbind(mse,mse_cv)
  }
  
  print(which.max(colSums(mse)))
  

  
  lbd<- mean(lambda_list[,which.max(colSums(mse))])*(n-nv)/n
  if(is.na(lbd)){
    lbd<- mean(lambda_list[,20])*(n-nv)/n
  }
  
  
  lbd<-c(lbd,lbd*phi)
  alpha_0<-0.5
  alpha.ori <- c(alpha_0,alpha_0*phi)
  lambda.ori <- c(lambda_list_1[ii],lambda_list_1[ii]*phi)
  alpha <- alpha.ori*lambda.ori
  lambda <- lambda.ori*(1-alpha.ori)/2;
  
  
  R.all<-cor(ss)
  
  tmp1 <- fast.lasso.sum.ess(b.vec, s.vec, r2.mat, n, group, lambda, alpha, init.vec, maxIter)
  
  res.cv1 <- as.vector(tmp1$beta)
  
  listM<-list(beta=tmp1$beta,l1_penalty=alpha[1],l2_penalty=lambda[1])
  
  
  
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
#' @param phi 
#'
#' @return
#' @export
#'
#' @examples
EXPRESSO.PSEUDO<-function(b.vec,cor.analysis,ss,n,group,maxIter,fold,phi){
  nt<-ceiling(n*(fold-1)/fold)
  nv<-n-nt
  r<-b.vec
  s.vec<-rep(1,length(r))
  v.vec <- 1/s.vec
  r2.mat<-cor.analysis
  init.vec<- array(0, dim = length(r))
  maxIter<-1000
  R.all<-cor(ss)
  lbd_p<-c()
  mu1<-rep(0,length(r))
  alpha_0<-0.5
  
  for(tt in 1:5){
    set.seed(tt)
    r_pseudo_ss<-mvrnorm(n,mu1, R.all)
    r_pseudo<-colSums(r_pseudo_ss)/n
    lbd_p[tt]<-max(abs(r_pseudo))
  }
  lbd<-max(lbd_p)/(n*alpha_0)
  
  lbd<-c(lbd,lbd*phi)
  
  alpha.ori <- c(alpha_0,alpha_0*phi)
  lambda.ori <- c(lbd,lbd*phi)
  alpha <- alpha.ori*lambda.ori
  lambda <- lambda.ori*(1-alpha.ori)/2;
  
  
  
  
  tmp1 <- fast.lasso.sum.ess(b.vec, s.vec, r2.mat, n, group, lambda, alpha, init.vec, maxIter)
  
  res.cv1 <- as.vector(tmp1$beta)
  
  listM<-list(beta=tmp1$beta,l1_penalty=alpha[1],l2_penalty=lambda[1])
  
  
  
}
