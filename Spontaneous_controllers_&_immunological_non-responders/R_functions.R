## Project: Gut Microbiome in Spontaneous HIV-1 Controllers and Immunological Non Responders: a metagenomics study in 1,559 people with HIV
## Author: Ángela del Castillo Izquierdo
## Last time edited: 2026-10-07
##
## Association functions used by scripts 1-4.
##
## Packages:
##   glmnet   4.1.10
##   reshape2 1.4.5
##   R.utils  2.13.0
##
## qtrans()     inverse-normal (rank) transform. "-999" is treated as missing.
## lm_phenotypes()  many phenotypes x many taxa, linear model.
## lm_phenotype()   one phenotype x many taxa, linear model.
## lrm_phenotypes() one phenotype x many taxa, logistic model (presence/absence).
## Every model inverse-normal transforms the outcome, the phenotype and the
## covariates, then fits lm() or glm(family = binomial). FDR is Benjamini-Hochberg
## across the taxa in that call.

library(glmnet)

## The Normal Quantile Transformation (inverse rank)
qtrans<-function(x){
  k<-!is.na(x)
  k<-which(x!="-999")
  ran<-rank(as.numeric(x[k]))
  y<-qnorm((1:length(k)-0.5)/length(k))
  x[k]<-y[ran]
  x
}

## linear model
lm_phenotypes<-function(y_mat,x_mat,cov_mat,covar){
  require(reshape2)
  require(R.utils)

  my_lm<-function(y,x){
    beta    <- NA
    se      <- NA
    p.value <- NA
    N       <- NA
    
    y_uniq_N        <- NA
    x_uniq_N        <- NA
    y_non_zero_N    <-NA
    x_non_zero_N    <-NA
    y_non_zero_rate <-NA
    x_non_zero_rate <-NA
    
    lm_input<-data.frame(Y = y, X = x,cov_mat[,covar]) %>% sapply(as.numeric) %>% na.omit
    N        <- nrow(lm_input)
    
    y_uniq_N   <- length(unique(lm_input[,"Y"])) 
    x_uniq_N   <- length(unique(lm_input[,"X"]))
    
    y_non_zero_N <- sum(!isZero(lm_input[,"Y"])) #number of participants with this bacteria
    x_non_zero_N <- sum(!isZero(lm_input[,"X"]))
    
    y_rate<-y_non_zero_N/length(unique(colnames(y_mat))) #rate of participants with this bacteria
    x_rate<-x_non_zero_N/N
    
    lm_input <- apply(lm_input, 2, qtrans) %>% as.data.frame
    try(lm_res <- summary(lm(Y~.,data = lm_input)), silent = T)
    indv<-'X'
    
    try(beta    <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),1],silent = T)
    try(se      <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),2], silent = T)
    try(p.value <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),4],silent = T)
    
    
    try(return(list(beta = beta,
                    se = se,
                    p.value = p.value,
                    N = N,
                    y_uniq_N     = y_uniq_N,
                    x_uniq_N     = x_uniq_N,
                    y_non_zero_N = y_non_zero_N,
                    x_non_zero_N = x_non_zero_N,
                    y_rate = y_rate,
                    x_rate= x_rate)),
        silent = T)
    
  }
  
  
  y_x<-sapply( 
    as.data.frame(y_mat),
    function(x) Map(function(a,b) my_lm(a,b),
                    list(x),
                    as.data.frame(x_mat)
    )
  )
  y_x.unlist <- matrix(unlist(y_x), ncol = 10, byrow = T)
  
  # beta matrix
  y_x.beta<-matrix(y_x.unlist[,1],ncol = ncol(y_mat), byrow = F)
  colnames(y_x.beta)<-colnames(y_mat)
  rownames(y_x.beta)<-colnames(x_mat)
  y_x_edge_beta  <- reshape2::melt(y_x.beta)
  
  y_x_edge<-data.frame(y_x_edge_beta[,c(2,1)],
                       as.data.frame(y_x.unlist),
                       fdr.p = p.adjust(y_x.unlist[,3], method = "fdr"),
                       bonferroni.p = p.adjust(y_x.unlist[,3], method = "bonferroni"))
  
  
  colnames(y_x_edge)<-c("Y", "X", "Beta","SE", "p","N","y_uniq_N","x_uniq_N", "y_non_zero_N", "x_non_zero_N", "y_rate","x_rate","fdr.p","bonferroni.p")
  
  return(y_x_edge)
}

## linear model with one phenotype
lm_phenotype<-function(y_mat,x_mat,cov_mat,covar){
  require(reshape2)
  require(R.utils)

  my_lm<-function(y,x){
    beta    <- NA
    se      <- NA
    p.value <- NA
    N       <- NA
    y_uniq_N        <- NA
    x_uniq_N        <- NA
    Phenotype1_sample_number <- NA
    Phenotype2_sample_number <- NA
    taxa_0_sample_number <- NA
    taxa_1_sample_number <- NA
    Phenotype1_sample_rate <- NA
    Phenotype2_sample_rate <- NA
    taxa_0_sample_rate <- NA
    taxa_1_sample_rate <- NA
    y_n_N    <-NA
    x_n_N    <-NA
    y_sample_rate <-NA
    x_sample_rate <-NA
    
    lm_input<-data.frame(Y = y, X = x, cov_mat[,covar]) %>% sapply(as.numeric) %>% na.omit
    colnames(lm_input)[1:2]=c("taxa","phenotype")
    N        <- nrow(lm_input)
    
    y_uniq_N   <- length(unique(lm_input[,"taxa"]))
    x_uniq_N   <- length(unique(lm_input[,"phenotype"]))
    
    Phenotype1_sample_number <- sum(lm_input[,"phenotype"]==c("1"),na.rm = T)
    Phenotype2_sample_number <- sum(lm_input[,"phenotype"]==c("0"),na.rm = T)
    
    taxa_0_sample_number = sum(lm_input[,"taxa"]==c("0"),na.rm = T)
    taxa_1_sample_number = sum(lm_input[,"taxa"]==c("1"),na.rm = T)
    
    taxa_0_sample_rate = taxa_0_sample_number/N
    taxa_1_sample_rate = taxa_1_sample_number/N
    
    Phenotype1_sample_rate<-Phenotype1_sample_number/N
    Phenotype2_sample_rate<-Phenotype2_sample_number/N
    
    y_n_N <- sum(!is.na(lm_input[,"taxa"]))
    x_n_N <- sum(!is.na(lm_input[,"phenotype"]))
    
    y_sample_rate<-y_n_N/N
    x_sample_rate<-x_n_N/N
    
    lm_input <- apply(lm_input, 2, qtrans) %>% as.data.frame
    try(lm_res <- summary(lm(taxa~.,data = lm_input)), silent = T)
    indv<-c("phenotype")
    
    try(beta    <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),1],silent = T)
    try(se      <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),2], silent = T)
    try(p.value <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),4],silent = T)
    
    
    try(return(list(beta = beta,
                    se = se,
                    p.value = p.value,
                    N = N,
                    y_uniq_N     = y_uniq_N,
                    x_uniq_N     = x_uniq_N,
                    y_n_N = y_n_N,
                    x_n_N = x_n_N,
                    y_sample_rate = y_sample_rate,
                    x_sample_rate= x_sample_rate,
                    Phenotype1_sample_number = Phenotype1_sample_number,
                    Phenotype2_sample_number = Phenotype2_sample_number,
                    taxa_0_sample_number = taxa_0_sample_number,
                    taxa_1_sample_number = taxa_1_sample_number,
                    Phenotype1_sample_rate = Phenotype1_sample_rate,
                    Phenotype2_sample_rate = Phenotype2_sample_rate,
                    taxa_0_sample_rate = taxa_0_sample_rate,
                    taxa_1_sample_rate = taxa_1_sample_rate)),
        silent = T)
    
  }
  
  
  y_x<-sapply( 
    as.data.frame(y_mat),
    function(x) Map(function(a,b) my_lm(a,b),
                    list(x),
                    as.data.frame(x_mat)
    )
  )
  y_x.unlist <- matrix(unlist(y_x), ncol = 18, byrow = T)
  
  # beta matrix
  y_x.beta<-matrix(y_x.unlist[,1],ncol = ncol(y_mat), byrow = F)
  colnames(y_x.beta)<-colnames(y_mat)
  rownames(y_x.beta)<-colnames(x_mat)
  y_x_edge_beta  <- reshape2::melt(y_x.beta)
  
  y_x_edge<-data.frame(y_x_edge_beta[,c(2,1)],
                       as.data.frame(y_x.unlist),
                       fdr.p = p.adjust(y_x.unlist[,3], method = "fdr"),
                       bonferroni.p = p.adjust(y_x.unlist[,3], method = "bonferroni"))
  
  
  colnames(y_x_edge)<-c("taxa", "phenotype", "Beta","SE", "p","N","y_uniq_N","x_uniq_N", 
                        "y_n_N", "x_n_N", "y_sample_rate","x_sample_rate",
                        "Phenotype1_sample_number","Phenotype2_sample_number", "taxa_0_sample_number", "taxa_1_sample_number",
                        "Phenotype1_sample_rate","Phenotype2_sample_rate", "taxa_0_sample_rate", "taxa_1_sample_rate",
                        "fdr.p","bonferroni.p")
  
  return(y_x_edge)
}

## logistic regression model
lrm_phenotypes<-function(y_mat,x_mat,cov_mat,covar, direction = c(1,1)){
  # mat0: phenotypes, mat1: microbiome, mat2: covar
  # direction: 1 means samples in row and variables in column; 2 means samples in column and variables in row
  
  require(reshape2)
  require(glmnet)
  require(R.utils)

  my_lm<-function(y,x){
    beta    <- NA
    se      <- NA
    p.value <- NA
    N       <- NA
    
    y_uniq_N        <- NA
    x_uniq_N        <- NA
    y_non_zero_N    <-NA
    x_non_zero_N    <-NA
    y_non_zero_rate <-NA
    x_non_zero_rate <-NA
    y_0_N<-NA
    y_1_N<-NA
    
    lm_input<-data.frame(Y = y, X = x,cov_mat[,covar]) %>% sapply(as.numeric) %>% na.omit

    N        <- nrow(lm_input) # number of participants
    
    y_uniq_N   <- length(unique(lm_input[,"Y"])) # number of variables in Y
    x_uniq_N   <- length(unique(lm_input[,"X"])) # number of variables in X
    
    y_non_zero_N <- sum(!isZero(lm_input[,"Y"])) # number of 1 in Y (no taxa - presence in total)
    x_non_zero_N <- sum(!isZero(lm_input[,"X"])) # number of 1 in X (no phenotype of interest)
    
    y_non_zero_rate<-y_non_zero_N/N
    x_non_zero_rate<-x_non_zero_N/N
    
    y_0_N<-sum(lm_input[,"Y"]==0) # number of taxa that is absent
    y_1_N<-sum(lm_input[,"Y"]==1) # number of taxa that is present
    
    lm_input_tmp <- apply(lm_input[,c(2:ncol(lm_input))], 2, qtrans) %>% as.data.frame
    lm_input<-data.frame(Y = lm_input[,1],lm_input_tmp)
    
    try(lm_res <- summary(glm(Y~., data = lm_input, family = 'binomial')))
    
    indv<-'X'
    
    try(beta    <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),1],silent = T)
    try(se      <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),2], silent = T)
    try(p.value <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),4],silent = T)
    # Wald p-value can underflow to 0 when |z| is very large (strong association / separation)
    if (!is.na(p.value) && (p.value <= 0 || p.value < .Machine$double.xmin))
      p.value <- .Machine$double.xmin

    try(return(list(beta = beta,
                    se = se,
                    p.value = p.value,
                    N = N,
                    y_uniq_N     = y_uniq_N,
                    x_uniq_N     = x_uniq_N,
                    y_non_zero_N = y_non_zero_N,
                    x_non_zero_N = x_non_zero_N,
                    y_non_zero_rate = y_non_zero_rate,
                    x_non_zero_rate= x_non_zero_rate,
                    y_0_N = y_0_N,
                    y_1_N = y_1_N)),
        silent = T)
    
  }
  
  
  y_x<-sapply( 
    as.data.frame(y_mat),
    function(x) Map(function(a,b) my_lm(a,b),
                    list(x),
                    as.data.frame(x_mat)
    )
  )
  y_x.unlist <- matrix(unlist(y_x), ncol = 12, byrow = T)
  
  # beta matrix
  y_x.beta<-matrix(y_x.unlist[,1],ncol = ncol(y_mat), byrow = F)
  colnames(y_x.beta)<-colnames(y_mat)
  rownames(y_x.beta)<-colnames(x_mat)
  y_x_edge_beta  <- reshape2::melt(y_x.beta)
  
  y_x_edge<-data.frame(y_x_edge_beta[,c(2:1)],
                       as.data.frame(y_x.unlist),
                       fdr.p = p.adjust(y_x.unlist[,3], method = "fdr"),
                       bonferroni.p = p.adjust(y_x.unlist[,3], method = "bonferroni"))
  
  
  colnames(y_x_edge)<-c("Y", "X", "Beta","SE", "p","N","y_uniq_N","x_uniq_N", "y_non_zero_N", "x_non_zero_N", "y_non_zero_rate","x_non_zero_rate","y_0_N","y_1_N","fdr.p","bonferroni.p")
  return(y_x_edge)
}


