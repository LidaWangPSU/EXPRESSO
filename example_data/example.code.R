#setwd()#yourlocal path

refFile <- "gtex_chr22.hm3"
annoFile <- "Whole_Blood_annotation_chr22.txt"
sumstatFile<- "Whole_Blood_GTEx_hm3_chr22_sumstat_beta.txt"

##Find gene list##
gene_subset = unique(read.table(sumstatFile, header = TRUE)$gene) 

windowFile <- c( "Whole_Blood_1mb_chr22.txt",
                 "Whole_Blood_250kb_chr22.txt",
                 "Whole_Blood_loop_chr22.txt",
                 "Whole_Blood_domain_chr22.txt",
                 "Whole_Blood_tad_chr22.txt",
                 "Whole_Blood_pchic_chr22.txt")

out_path= "output/EXPRESSO"


res.tmp <- EXPRESSO(sumstatFile,annoFile,windowFile,refFile,out_path,minMaf=0.05,maxIter=1000,gene.vec = gene_subset,append=F)

