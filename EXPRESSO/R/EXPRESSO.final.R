#' Title
#'
#' @param root 
#' @param impute 
#'
#' @return
#' @export
#'
#' @examples
read_plink_custom <- function(root, impute = c('none', 'avg', 'random')) {
  if(impute == 'random') {
    stop("The 'impute' random option has not been implemented.", call. = FALSE)
  }
  
  ## structure from https://github.com/gabraham/plink2R/blob/master/plink2R/R/plink2R.R
  proot <- path.expand(root)
  
  bedfile <- paste(proot, ".bed", sep="")
  famfile <- paste(proot, ".fam", sep="")
  bimfile <- paste(proot, ".bim", sep="")
  
  ## Could change this code to use data.table
  bim <- read.table(bimfile, header=FALSE, sep="", stringsAsFactors=FALSE)
  fam <- read.table(famfile, header=FALSE, sep="", stringsAsFactors=FALSE)
  ## Set the dimensions
  geno <- BEDMatrix::BEDMatrix(bedfile, n = nrow(fam), p = nrow(bim))
  
  ## Convert to a matrix
  geno <- as.matrix(geno)
  if(impute == 'avg') {
    ## Check if any are missing
    geno_na <- is.na(geno)
    if(any(geno_na)) {
      means <- colMeans(geno, na.rm = TRUE)
      geno[geno_na] <- rep(means, colSums(geno_na))
    }
  }
  colnames(geno) <- bim[,2]
  rownames(geno) <- paste(fam[,1], fam[, 2], sep=":")
  
  list(bed=geno, fam=fam, bim=bim)
}



#' Title
#'
#' @param sumstatFile 
#' @param annoFile 
#' @param windowFile 
#' @param refFile 
#' @param out_path 
#' @param minMaf 
#' @param refIndv 
#' @param maxIter 
#' @param gene.vec 
#' @param append 
#'
#' @return
#' @export
#'
#' @examples
EXPRESSO<- function(sumstatFile, annoFile,windowFile,refFile,out_path,minMaf=0.05,refIndv=NULL,maxIter=1000,gene.vec=NULL,append=TRUE) {
  fname.out1 <- paste0(out_path,'-weight.mse.txt');
  fname.out2 <- paste0(out_path,'-weight.bic.txt');
  fname.out3 <- paste0(out_path,'-weight.pseudo.txt');
  cat(file=fname.out1,c('gene','snp','weight'),'\n',sep='\t',append=append);
  cat(file=fname.out2,c('gene','snp','weight'),'\n',sep='\t',append=append);
  cat(file=fname.out3,c('gene','snp','weight'),'\n',sep='\t',append=append);
  
  cv.out1 <- paste0(out_path,'-cv.mse.txt');
  cv.out2 <- paste0(out_path,'-cv.bic.txt');
  cv.out3 <- paste0(out_path,'-cv.pseudo.txt');
  cat(file=cv.out1,c('gene','window','phi','r2'),'\n',sep='\t',append=append);
  cat(file=cv.out2,c('gene','window','phi','r2'),'\n',sep='\t',append=append);
  cat(file=cv.out3,c('gene','window','phi','r2'),'\n',sep='\t',append=append);
  
  
  
  window <- list();
  for(ii in 1:length(windowFile)) {
    
    if(file.exists(windowFile[[ii]])){
    window[[ii]] <- as.data.frame(fread(windowFile[[ii]],header=TRUE))
    }else{
      window[[ii]]<-NA}
  }
  names(window) <- basename(windowFile);
  
  if(is.null(annoFile)){
  anno<-as.data.frame(t(c(0,0,0,0,0,0)))
  colnames(anno)<-c("snp", "chr","pos","ref","alt","overlap")
  }else{
  anno <- as.data.frame(fread(annoFile,header=T));
  colnames(anno)<-c("snp", "chr","pos","ref","alt","overlap")
  
  }
  sumstat1 <- as.data.frame(fread(sumstatFile,header=T));
  colnames(sumstat1)<-c("gene","symbol","chr","pos","ref","alt","af","beta.est","n","snp_id")
  head(sumstat1)
  
  ix.common <- which(sumstat1$af>minMaf & sumstat1$af<(1-minMaf));
  sumstat1 <- sumstat1[ix.common,];
  
  if(is.null(gene.vec))
    gene.vec <- unique(sumstat1$gene);
  
  
  ix.min.mse<-c()
  ix.min.bic<-c()
  ix.min.pseudo<-c()
  
  
  ref = read_plink_custom(refFile, impute = 'avg')
  genoref = ref$bed
  genoref = as.data.frame( t(genoref) )
  
  genoref[1:5,1:10]
  colnames(genoref)<-gsub(":","_",colnames(genoref))
 
  pos.ref1<-paste0(ref$bim$V1,":",ref$bim$V4)
  refalt.ref1<-paste0(ref$bim$V5,"_",ref$bim$V6)
  
  for(ii in 1:length(gene.vec)) {
    
    sumstat<-sumstat1[which(sumstat1[,1]==gene.vec[ii]),]
    
    pos.sumstat <- paste(sumstat$chr,sumstat$pos,sep=":");
    pos.anno <- paste(anno$chr,anno$pos,sep=":");
    
    
    cat('Analysing',gene.vec[ii],'\n');
    ix.var <- which(sumstat[,1]==gene.vec[ii]);
    pos.var <- paste(sumstat$chr[ix.var],sumstat$pos[ix.var],sep=":");
    pos.var <-unique(pos.var)
    chr.var <- as.numeric(unique(gsub(":.*","",pos.var)));
    
    index<-na.omit(match(pos.var,pos.ref1))
    
    genoref_1<-genoref[index,]
    
    pos.ref<-pos.ref1[index]
    refalt.ref<-refalt.ref1[index]
    
    cor.ref <- cor(t(genoref_1))
    
    colnames(cor.ref) <- pos.ref;
    rownames(cor.ref) <- pos.ref
    
    
    pos.ss<-intersect(pos.var,pos.ref)
    
    ix.match <- match(pos.ss,pos.sumstat);
    
    if(length(ix.match)==0){next}
    
    ss.simu <- sumstat[ix.match,];
    ss.r<-ss.simu$beta.est
    
    n<-ceiling(median(ss.simu$n))
    
    set.seed(123)
    
    ss_all<-mvrnorm(n,ss.r,cor.ref)
    
    res.lasso <- list();
    cv.list<-data.frame()
    
    b.vec.list<-list()
    cor.analysis.list<-list()
    group.list<-list()
    pos.list<-list()
    ss.list<-list()
    
    for(jj in 1:length(window)) {
      if(is.na(window[[jj]])){next}
      cat('analyzing',names(window)[jj],'\n',sep=" ");
      a <- Sys.time();
      pos.ref<-gsub("X",23,pos.ref)
      pos.anno<-gsub("X",23,pos.anno)
      gene.window <- window[[jj]][which(window[[jj]]$gene==gene.vec[ii]),][1,]
      pos.var.in.window <- getVarInWindow(pos.var,gene.window);
      
      pos.analysis <- intersect(pos.var.in.window,pos.ss);
      
      ix.match <- match(pos.analysis,pos.ref);
      
      cor.analysis <- matrix(cor.ref[ix.match,ix.match],nrow=length(ix.match),ncol=length(ix.match));
      diag(cor.analysis) <- 1;
      
      
      ix.match <- match(pos.analysis,pos.ss);
      
      sumstat.analysis <- ss.simu[ix.match,];
      ss<-ss_all[,ix.match]
      
      
      ix.match <- match(pos.analysis, pos.anno);
      anno.analysis <- anno[ix.match,];
      ess <- rm.na(anno.analysis$overlap);
      group <- rep(1,length(pos.analysis));
      group[which(ess==1)] <- 2;
      
      b.vec <- as.numeric(sumstat.analysis$beta.est);
      snp<-sumstat.analysis$snp_id
      
      cat("b.vec:")
      print(summary(b.vec))
      if(length(b.vec)>1) {
        
        fold<-4
        res.lasso[[jj]] <- EXPRESSO.cv(b.vec,cor.analysis,ss,n,group,maxIter,fold);
        windowName<-rep(gene.window$type,nrow(res.lasso[[jj]]))
        temp<-cbind(windowName,res.lasso[[jj]])
        cv.list<-rbind(cv.list,temp)
        
        l1<-list(b.vec)
        names(l1)<-gene.window$type[1]
        b.vec.list<-append(b.vec.list, l1)
        
        l2<-list(cor.analysis)
        names(l2)<-gene.window$type[1]
        cor.analysis.list<-append(cor.analysis.list, l2)
        
        l3<-list(group)
        names(l3)<-gene.window$type[1]
        group.list<-append(group.list, l3)
        
        l4<-list(snp)
        names(l4)<-gene.window$type[1]
        pos.list<-append(pos.list, l4)
        
        l5<-list(ss)
        names(l5)<-gene.window$type[1]
        ss.list<-append(ss.list, l5)
      }
      b <- Sys.time();
      cat('window finishes in',format(b-a),'\n');
      
    }
    
    if(nrow(cv.list)==0){next}
    
    ix.min.mse[ii] <- which.max(cv.list$mse_r2);
    ix.min.bic[ii] <- which.max(cv.list$bic_r2);
    ix.min.pseudo[ii] <- which.max(cv.list$pseudo_r2);
    
    
    cv.out.mse<-cbind(gene.vec[ii],cv.list[ix.min.mse[ii],c(1,2,3)])
    cv.out.bic<-cbind(gene.vec[ii],cv.list[ix.min.bic[ii],c(1,2,4)])
    cv.out.pseudo<-cbind(gene.vec[ii],cv.list[ix.min.pseudo[ii],c(1,2,5)])
    
    
    write.table(file=cv.out1,cv.out.mse,sep='\t',append=TRUE,col.names=F,row.names=F,quote=F)
    write.table(file=cv.out2,cv.out.bic,sep='\t',append=TRUE,col.names=F,row.names=F,quote=F)
    write.table(file=cv.out3,cv.out.pseudo,sep='\t',append=TRUE,col.names=F,row.names=F,quote=F)
    
    #########computing weights
    
    
    b.vec.mse<-as.numeric(b.vec.list[[paste0(cv.out.mse$windowName)]])
    cor.analysis.mse<-as.matrix(cor.analysis.list[[paste0(cv.out.mse$windowName)]])
    group.mse<-as.numeric(group.list[[paste0(cv.out.mse$windowName)]])
    phi.mse<-as.numeric(cv.out.mse$phi)
    pos.mse<-pos.list[[paste0(cv.out.mse$windowName)]]
    ss.mse<-ss.list[[paste0(cv.out.mse$windowName)]]
    mse.out<-EXPRESSO.MSE(b.vec.mse,cor.analysis.mse,ss.mse,n,group.mse,maxIter,fold,phi.mse)
    
    
    if(length(mse.out$beta)>0) {
      ix.keep.b <- which(mse.out$beta!=0);
      beta.out.b <- mse.out$beta[ix.keep.b];
      #print(beta.out.b)
      gene.out.b <- rep(gene.vec[ii],length(beta.out.b));
      snp.out.b <- pos.mse[ix.keep.b];
      window.out.b <- rep(cv.out.mse$windowName,length(ix.keep.b));
      phi.out.b <- rep(phi.mse,length(ix.keep.b));
      l1.out.b <- rep(mse.out$l1_penalty,length(ix.keep.b));
      l2.out.b <- rep(mse.out$l2_penalty,length(ix.keep.b));
      res.out.b <- cbind(gene.out.b,snp.out.b,beta.out.b);
      write.table(file=fname.out1,res.out.b,sep='\t',append=TRUE,col.names=F,row.names=F,quote=F)

    }
    
    
    
    
    
    b.vec.bic<-as.numeric(b.vec.list[[paste0(cv.out.bic$windowName)]])
    cor.analysis.bic<-as.matrix(cor.analysis.list[[paste0(cv.out.bic$windowName)]])
    group.bic<-as.numeric(group.list[[paste0(cv.out.bic$windowName)]])
    phi.bic<-as.numeric(cv.out.bic$phi)
    pos.bic<-pos.list[[paste0(cv.out.bic$windowName)]]
    ss.bic<-ss.list[[paste0(cv.out.bic$windowName)]]
    bic.out<-EXPRESSO.BIC(b.vec.bic,cor.analysis.bic,ss.bic,n,group.bic,maxIter,fold,phi.bic)
    
    
    if(length(bic.out$beta)>0) {
      ix.keep.b <- which(bic.out$beta!=0);
      beta.out.b <- bic.out$beta[ix.keep.b];
      #print(beta.out.b)
      gene.out.b <- rep(gene.vec[ii],length(beta.out.b));
      snp.out.b <- pos.bic[ix.keep.b];
      window.out.b <- rep(cv.out.bic$windowName,length(ix.keep.b));
      phi.out.b <- rep(phi.bic,length(ix.keep.b));
      l1.out.b <- rep(bic.out$l1_penalty,length(ix.keep.b));
      l2.out.b <- rep(bic.out$l2_penalty,length(ix.keep.b));
      res.out.b <- cbind(gene.out.b,snp.out.b,beta.out.b);
      write.table(file=fname.out2,res.out.b,sep='\t',append=TRUE,col.names=F,row.names=F,quote=F)

    }else if(length(mse.out$beta)>0){
      write.table(file=fname.out2,res.out.b,sep='\t',append=TRUE,col.names=F,row.names=F,quote=F)

    }
    
    
    
    b.vec.pseudo<-as.numeric(b.vec.list[[paste0(cv.out.pseudo$windowName)]])
    cor.analysis.pseudo<-as.matrix(cor.analysis.list[[paste0(cv.out.pseudo$windowName)]])
    group.pseudo<-as.numeric(group.list[[paste0(cv.out.pseudo$windowName)]])
    phi.pseudo<-as.numeric(cv.out.pseudo$phi)
    pos.pseudo<-pos.list[[paste0(cv.out.pseudo$windowName)]]
    ss.pseudo<-ss.list[[paste0(cv.out.pseudo$windowName)]]
    pseudo.out<-EXPRESSO.PSEUDO(b.vec.pseudo,cor.analysis.pseudo,ss.pseudo,n,group.pseudo,maxIter,fold,phi.pseudo)
    
    if(length(pseudo.out$beta)==0&length(bic.out$beta)!=0){
      pseudo.out<-bic.out
    }
    
    if(length(pseudo.out$beta)==0&length(mse.out$beta)!=0){
      pseudo.out<-mse.out
    }
    
    if(length(pseudo.out$beta)>0) {
      ix.keep.b <- which(pseudo.out$beta!=0);
      beta.out.b <- pseudo.out$beta[ix.keep.b];
      #print(beta.out.b)
      gene.out.b <- rep(gene.vec[ii],length(beta.out.b));
      snp.out.b <- pos.pseudo[ix.keep.b];
      window.out.b <- rep(cv.out.pseudo$windowName,length(ix.keep.b));
      phi.out.b <- rep(phi.pseudo,length(ix.keep.b));
      l1.out.b <- rep(pseudo.out$l1_penalty,length(ix.keep.b));
      l2.out.b <- rep(pseudo.out$l2_penalty,length(ix.keep.b));
      res.out.b <- cbind(gene.out.b,snp.out.b,beta.out.b);
      write.table(file=fname.out3,res.out.b,sep='\t',append=TRUE,col.names=F,row.names=F,quote=F)

    }else if(length(bic.out$beta)>0){
      write.table(file=fname.out3,res.out.b,sep='\t',append=TRUE,col.names=F,row.names=F,quote=F)
    }else if(length(mse.out$beta)>0){
      write.table(file=fname.out3,res.out.b,sep='\t',append=TRUE,col.names=F,row.names=F,quote=F)
    }
    
  }
  
}
## tabix.range <- get.tabix.range(pos.var);
## ##     cfh <- readSingleChromosomeVCFToMatrixByRange(fileName, "1:196621007-196716634")
## a1 <- Sys.time();refMat <- readSingleChromosomeVCFToMatrixByRange(refFile[chr.var], tabix.range);b1 <- Sys.time();
## a2 <- Sys.time();refMat1 <- readVCFToMatrixByRange(refFile[chr.var], tabix.range,"");b2 <- Sys.time();

#' get windows;
#'
#' @param pos.var the variant positions surrouding the transcript;
#' @param gene.window the windows defined for the gene;
#' @export 
getVarInWindow <- function(pos.var,gene.window) {
  pos <- as.numeric(gsub(".*:","",pos.var));
  window.start <- gene.window$pos_start;
  window.end <- gene.window$pos_end;
  ix <- 0;
  for(ii in 1:length(pos.var)) {
    ix[ii] <- (length(which(pos[ii]>=window.start & pos[ii]<=window.end))>0);
  }
  ##print(ix);
  return(pos.var[which(ix==1)]);
}

#' cor2cov
#'
#' @param cor.mat correlation matrix;
#' @param sd.vec standard deviation;
#' @export
cor2cov <- function (cor.mat, sd.vec) {
  return(t(t(cor.mat * sd.vec) * sd.vec))
}
