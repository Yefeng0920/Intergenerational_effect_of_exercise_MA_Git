#-----------------------------------------------------------------------------------------#
# Functions used throughout analyses and processing of script for
# Intergenerational effects of physical exercise on cognition and brain: 
# a multilevel meta-analysis of mean and variance
# Author: Yefeng Yang
# Date update: August 2020
#-----------------------------------------------------------------------------------------#


#============================================================================ 
# 'Title: I2 based on Wolfgang and Nakagawa & Santos, 2012'
#============================================================================ 

## Following: (i) Wolfgang Viechtbauer: http://www.metafor-project.org/doku.php/tips:i2_multilevel_multivariate
##            (ii) Nakagawa S, Santos E S A. Methodological issues and advances in biological meta-analysis[J]. Evolutionary Ecology, 2012, 26(5): 1253-1274.


I2 <- function(model, method = c("Wolfgang", "Shinichi")) {
  
  ## evaluate choices
  method <- match.arg(method)
  
  # Wolfgang's method
  if (method == "Wolfgang") {
    W <- solve(model$V)
    X <- model.matrix(model)
    P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
    I2_total <- sum(model$sigma2)/(sum(model$sigma2) + (model$k - model$p)/sum(diag(P)))
    I2_each <- model$sigma2/(sum(model$sigma2) + (model$k - model$p)/sum(diag(P)))
    names(I2_each) = paste0("I2_", model$s.names)
    
    # putting all together
    I2s <- round(c(I2_total = I2_total, I2_each)*100, 3)
    
    # or our way (Shinichi)
  } else {
    # sigma2_v = typical sampling error variance
    sigma2_v <- sum(1/model$vi) * (model$k - 1)/(sum(1/model$vi)^2 - sum((1/model$vi)^2))
    I2_total <- sum(model$sigma2)/(sum(model$sigma2) + sigma2_v)  #s^2_t = total variance
    I2_each <- model$sigma2/(sum(model$sigma2) + sigma2_v)
    names(I2_each) = paste0("I2_", model$s.names)
    
    # putting all together
    I2s <- round(c(I2_total = I2_total, I2_each)*100, 3)
  }
  return(I2s)
}


#============================================================================ 
# 'Title: R2 based on Nakagawa & Schielzeth, 2013'
#============================================================================ 

## Following: Nakagawa S, Schielzeth H. A general and simple method for obtaining R2 from generalized linear mixed‐effects models[J]. Methods in ecology and evolution, 2013, 4(2): 133-142.


R2 <- function(model) {
  # fixed effect variance
  fix <- var(as.numeric(as.vector(model$b) %*% t(as.matrix(model$X))))
  # marginal
  R2m <- fix / (fix + sum(model$sigma2))
  # Rm <- round(R2m, 3)
  # conditional
  R2c <- (fix + sum(model$sigma2) - model$sigma2[length(model$sigma2)]) /
    (fix + sum(model$sigma2))
  R2s <- c(R2_marginal = R2m, R2_coditional = R2c)
  return(R2s)
}




#============================================================================ 
# 'Title: Converting percentage differences'
#============================================================================ 

per.transformation <- function(model, digits = 1){
  ifelse(model$b>0,
         per <- round((exp(model$b)-1)*100, digits), 
         per <- round((1-exp(model$b))*(-100), digits))
}



#============================================================================ 
# 'Title: Conceptual distribution function (tailed for function)'
#============================================================================ 

rnorm.exercise <- function(n, mean, sd, ll = -Inf, ul = Inf){
    qnorm(runif(n, pnorm(ll, mean, sd), pnorm(ul, mean, sd)), mean, sd)
}



#============================================================================ 
# 'Title: model selection with rma.mv in MuMIn (made by *Kamil Barto)'
#============================================================================ 

## Following: # *Kamil Bartoń: https://CRAN.R-project.org/package=MuMIn /


formula.rma.mv <- function(x, ...) return(eval(getCall(x)$mods))

makeArgs.rma.mv <- function(obj, termNames, comb, opt, ...) {
  ret <- MuMIn:::makeArgs.default(obj, termNames, comb, opt)
  names(ret)[1L] <- "mods"
  ret
}

nobs.rma.mv <- function(object, ...) attr(logLik(object), "nall")

coefTable.rma.mv <- function(model, ...) MuMIn:::.makeCoefTable(model$b, model$se, 
                                                                coefNames = rownames(model$b))



#============================================================================ 
# 'Title: variance-covariance matrix (made by Lagisz Malgorzata & Shinichi Nakagawa)'
#============================================================================ 

## Following: Noble D W A, Lagisz M, O'dea R E, et al. Nonindependence and sensitivity analyses in ecological and evolutionary meta‐analyses[J]. Molecular Ecology, 2017, 26(9): 2410-2425.


make_VCV_matrix <- function(data, V, cluster, obs, type=c("vcv", "cor"), rho=0.5){
  
  if (missing(data)) 
    stop("Must specify dataframe via 'data' argument.")
  if (missing(V)) 
    stop("Must specify name of the variance variable via 'V' argument.")
  if (missing(cluster)) 
    stop("Must specify name of the clustering variable via 'cluster' argument.")
  if (missing(obs)) 
    obs <- 1:length(V)   
  if (missing(type)) 
    type <- "vcv" 
  
  new_matrix <- matrix(0,nrow = dim(data)[1],ncol = dim(data)[1]) 
  rownames(new_matrix) <- data[ ,obs]
  colnames(new_matrix) <- data[ ,obs]
  shared_coord <- which(data[ ,cluster] %in% data[duplicated(data[ ,cluster]), cluster]==TRUE)
  combinations <- do.call("rbind", tapply(shared_coord, data[shared_coord,cluster], function(x) t(combn(x,2))))
  
  if(type == "vcv"){
    # calculate covariance values between  values at the positions in shared_list and place them on the matrix
    for (i in 1:dim(combinations)[1]){
      p1 <- combinations[i,1]
      p2 <- combinations[i,2]
      p1_p2_cov <- rho * sqrt(data[p1,V]) * sqrt(data[p2,V])
      new_matrix[p1,p2] <- p1_p2_cov
      new_matrix[p2,p1] <- p1_p2_cov
    }
    diag(new_matrix) <- data[ ,V]
  }
  
  if(type == "cor"){
    for (i in 1:dim(combinations)[1]){
      p1 <- combinations[i,1]
      p2 <- combinations[i,2]
      p1_p2_cov <- rho
      new_matrix[p1,p2] <- p1_p2_cov
      new_matrix[p2,p1] <- p1_p2_cov
    }
    diag(new_matrix) <- 1 
  }
  
  return(new_matrix)
}