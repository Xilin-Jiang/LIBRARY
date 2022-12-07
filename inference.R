# an example inference of mixture of dirichlet distribution

###########################################
# implementing baseline LDA: mean-field VB
###########################################
comp_lda_lb <- function(para){
  # compute the lower bound for the whole dataset

  # terms for cross-entropy
  term1 <- para$M*(lgamma(sum(para$alpha)) - sum(lgamma(para$alpha)) ) +
    sum( para$E_lntheta %*% (para$alpha - 1) )
  term2 <- sapply(1:para$M,
                  function(s) sum( para$E_zn[[s]] %*% para$E_lntheta[s,] )) %>% sum
  # term3 <- sum(sapply(1:para$M, function(s) sum(para$E_zn[[s]] * log(para$beta_w[[s]]) ) ) )
  # term3 could be fully vectorized using para$unlist_zn and para$beta_w_full
  # term3 <- sum(para$unlist_zn * log(para$beta_w_full) )
  term3 <- sum( log(para$beta_w_full^para$unlist_zn) ) # use this method to avoid numeric issue

  # terms for entropy
  term4 <- sum(sapply(1:para$M,
                      function(s) lgamma(sum(para$alpha_z[s,])) - sum(lgamma(para$alpha_z[s,])))) +
    sum((para$alpha_z - 1)*para$E_lntheta)
  # term5 <- sum(sapply(1:para$M, function(s) sum( para$E_zn[[s]]*log(para$E_zn[[s]]) ) ) )
  # term5 could be fully vectorized using para$unlist_zn
  # term5 <- sum( para$unlist_zn*log(para$unlist_zn) )
  term5 <- sum( log(para$unlist_zn^para$unlist_zn) ) # similarly, avoid numeric issue

  return((term1 + term2 + term3 - term4 - term5))
}

comp_E_zn <- function(para){
  # compute E-step for zn
  # two layer of apply: outer layer is for documents s=1,..,M; saplly layer is for words in documents n=1,2,...Ns
  para$E_zn <-sapply(1:para$M,
                     function(s)
                       (para$beta_w[[s]] %*% diag(exp(para$E_lntheta[s,]))  )/
                       (para$beta_w[[s]]) %*% t(exp(para$E_lntheta[rep(s,para$K),])),
                     simplify = FALSE)
  # update the variables to facilitate computation
  para$alpha_z <- sapply(1:para$M, function(s) para$alpha + colSums(para$E_zn[[s]])) %>% t
  para$unlist_zn <- do.call(rbind, para$E_zn)
  return(para)
}

comp_E_lntheta <- function(para){
  # compute E-step for ln(theta)
  para$E_lntheta <- sapply(1:para$M,
                           function(s) digamma(para$alpha_z[s,]) -
                             digamma(sum(para$alpha_z[s,]))) %>% t
  return(para)
}

update_beta_basic_lda <- function(para){
  # compute M-step for beta: basic case, direct maximize the upper bound
  para$beta <- sapply(1:para$D, function(j) colSums(para$unlist_zn[para$ds_list[[j]]$id,]) ) %>% t
  # normalize beta
  para$beta <- sapply(1:para$K, function(i) para$beta[,i]/sum(para$beta[,i]))
  # this beta_w parameter save the beta for each word: it is a list of M elements and each contain a K*Ns matrix
  # para$beta_w <- lapply(para$w, function(w) para$beta[w$Ds_id,,drop=FALSE] )
  para$beta_w_full <- para$beta[para$unlist_Ds_id$Ds_id,,drop=FALSE]
  para$beta_w <- lapply(para$patient_lst, function(x) para$beta_w_full[x,,drop=F])
  return(para)
}

update_alpha <- function(para){
  # compute M-step for alpha; optimize the dirichlet with Newton-Raphson method
  para$lb_alpha <- function(alpha){
    para$M*(lgamma(sum(alpha)) - sum(lgamma(para$alpha)) ) +
      sum(colSums(para$E_lntheta) * (para$alpha - 1) )
  }

  para$grad_alpha <- function(alpha){
    para$M*(digamma(sum(alpha)) - digamma(para$alpha))  +
      colSums( para$E_lntheta )
  }

  para$hess_alpha <- function(alpha){
    para$M*( trigamma(sum(alpha)) - diag(trigamma(alpha)) )
  }
  # para$optim_alpha <- maxBFGS(para$lb_alpha, grad = para$grad_alpha, start = para$alpha,
  #                             constraints=list(ineqA = diag(nrow = para$K), ineqB = matrix(0,nrow = para$K)) )
  para$optim_alpha <- maxNR(para$lb_alpha, grad = para$grad_alpha, hess = para$hess_alpha, start = para$alpha) # rep(10^(-5), para$K))

  para$alpha <- para$optim_alpha$estimate

  return(para)
}

# example
rec_data <- read.csv("rec2subjectAbove500occur_include_death_ICDA2N.csv")
ds_list <- read.csv("listAbove500include_deaths_ICDA2N.csv")
topic_num <- 10
para <- topic_init_baseline(rec_data, ds_list, topic_num)
# set the number of update
para$max_itr <- 500
para$lb <- data.frame("Iteration" = 0,"Lower_bound" = comp_lda_lb(para))
para$alpha <- rep(1, para$K)
para$tol <- 10^(-7) # tolerance of lower bound step
para$itr_check_lb <- 1
for(itr in 1:para$max_itr){
  print(paste0("Interation: ",itr))
  para <- comp_E_zn(para)
  para <- comp_E_lntheta(para)
  para <- update_beta_basic_lda(para)
  # print(paste0("Current Lower bound at beta ", comp_lda_lb(para)))
  # para <- update_alpha(para) # we do not need to update alpha
  # print(paste0("Current Lower bound at alpha ", comp_lda_lb(para)))
  para$lb[nrow(para$lb) + 1,] <- c(itr, comp_lda_lb(para))
  if(itr %% para$itr_check_lb  ==0){
    print(paste0("Current Lower bound ", pull(filter(para$lb, Iteration == itr), Lower_bound), " at iteration: ",itr))
    curr_lb <- pull(filter(para$lb, Iteration == itr), Lower_bound)
    prev_lb <- pull(filter(para$lb, Iteration == (itr -para$itr_check_lb )), Lower_bound)
    if(is.finite((curr_lb - prev_lb)) & (curr_lb - prev_lb)/abs(prev_lb) < para$tol ){
      print(paste0("Optimization converged at step ", itr))
      break
    }
  }
}


# save the parameter
save(para, file = paste0("~/Desktop/comorbidity/Results/","Run_",Sys.Date(),".RData"))


##############################################
# testing for subtypes functions: Mixture of dirichlet distribution
##############################################
comp_MixDir_pi <- function(parMixDir) {
  # update pi
  parMixDir$pi <- colSums(parMixDir$z_n)/sum(parMixDir$z_n)
  return(parMixDir)
}

comp_MixDir_zn <- function(parMixDir) {
  # update z_n
  z_n <- sapply(1:parMixDir$K, function(i) parMixDir$pi[i]*ddirichlet(parMixDir$X, parMixDir$alpha[i,]))
  parMixDir$z_n <- z_n/rowSums(z_n)
  return(parMixDir)
}

logL_full_MixDir <- function(alphal, parMixDir, l) {
  # term1 is the multivariate beta distribution multiplied by sum of gammas
  term1 <- (sum(lgamma(alphal)) - lgamma(sum(alphal))) * parMixDir$GammaSumN[l]
  term2 <- (alphal - 1) %*% parMixDir$GammaLogZn[l, ]
  # print(c(term2 , term1, lgamma(sum(alphal)), alphal) )
  return(term2 - term1)
}
grad_logL_mixDir <- function(alphal, parMixDir, l) {
  # term1 is the multivariate beta distribution multiplied by sum of gammas
  term1 <- (digamma(sum(alphal)) - digamma(alphal)) * parMixDir$GammaSumN[l]
  term2 <- parMixDir$GammaLogZn[l, ]
  # print(c(term2 , term1, lgamma(sum(alphal)), alphal) )
  return(term1 + term2)
}
########################
# maybe we will need gradient information
########################

comp_MixDir_alpha <- function(parMixDir){
  # update each alpha in turn, using full likelihood
  parMixDir$GammaSumN <- colSums(parMixDir$z_n)
  parMixDir$GammaLogZn <- crossprod(parMixDir$z_n, log(parMixDir$X))
  for(l in 1:parMixDir$K){
    parMixDir$alpha[l,] <- optim(par = parMixDir$alpha[l,],
                                 fn = function(x) logL_full_MixDir(x, parMixDir, l),
                                 gr = function(x) grad_logL_mixDir(x, parMixDir, l),
                                 control = list(fnscale = -1),
                                 lower = rep(1, parMixDir$M),  upper = rep(10, parMixDir$M),
                                 method="L-BFGS-B")$par
  }
  return(parMixDir)
}

# compute likelihood for convergence
logL_MixDir <-  function(parMixDir) {
  terms <- 0
  for(l in 1:parMixDir$K){
    terms <- terms + parMixDir$pi[l] * ddirichlet(parMixDir$X,parMixDir$alpha[l,])
  }
  return(sum(log(terms)))
}

fit_MixDir <- function(loadings, K){
  parMixDir <- list()
  parMixDir$K <- K
  parMixDir$N <- dim(loadings)[1]
  parMixDir$M <- dim(loadings)[2]
  parMixDir$pi <- runif(parMixDir$K)
  parMixDir$pi <- parMixDir$pi/sum(parMixDir$pi)
  parMixDir$X <- loadings
  parMixDir$alpha <- matrix(runif(parMixDir$K * parMixDir$M), nrow = parMixDir$K, ncol = parMixDir$M)
  parMixDir$max_itr <- 200
  parMixDir$tol <- 10^(-6)
  parMixDir$lb <- data.frame("Iteration" = 0,"Lower_bound" = logL_MixDir(parMixDir))
  for(itr in 1:parMixDir$max_itr){
    # print(paste0("Interation: ",itr))
    parMixDir <- comp_MixDir_zn(parMixDir)
    # print(paste0("after zn", logL_MixDir(parMixDir)))
    parMixDir <- comp_MixDir_pi(parMixDir)
    # print(paste0("after pi", logL_MixDir(parMixDir)))
    parMixDir <- comp_MixDir_alpha(parMixDir)
    # print(paste0("after alpha", logL_MixDir(parMixDir)))
    # para <- update_alpha(para) # we use non-informative alpha
    parMixDir$lb[nrow(parMixDir$lb) + 1,] <- c(itr, logL_MixDir(parMixDir))
    if(itr %% 5 ==0){
      curr_lb <- pull(filter(parMixDir$lb, Iteration == itr), Lower_bound)
      prev_lb <- pull(filter(parMixDir$lb, Iteration == (itr -para$itr_check_lb )), Lower_bound)
      # print(paste0("Current Lower bound ", curr_lb, " at iteration: ",itr))
      try({
        if(is.finite((curr_lb - prev_lb)) & (curr_lb - prev_lb)/abs(prev_lb) < para$tol ){
          print(paste0("Optimization converged at step ", itr))
          break
        }
      })
    }
  }
  return(parMixDir)
}

########################
# treeLDA
########################
topic_init_tree <- function(rec_data, ds_list, topic_num, tree_leveles){
  # arrange the data stack by individual is very important as we need to make sure the matrix could be rejoined into a single matrix
  first_incidence_age <- rec_data %>%
    arrange(eid)

  # plot the number distribution of indiviudal diseases
  df_number_records <- first_incidence_age %>%
    group_by(eid) %>%
    summarise(n())

  para <- list()
  para$eid <- df_number_records$eid

  # add death to it
  para$list_above500occu <- ds_list
  # para$list_above500occu <- read.csv(paste0("listAbove500.csv"))
  para$D <- dim(para$list_above500occu)[1] # disease number
  para$M <- length(para$eid) # subject number
  para$K <- topic_num # start with 10 component
  para$L <- tree_leveles # for ICD10 it is 4 layers
  # also need to compute record number per individual for cvb
  para$Ns <- df_number_records$`n()`

  code2id <- function(x){
    return( match(x, para$list_above500occu$diag_icd10))
  }

  # here I am rounding the disease time to year for computation efficiency
  para$unlist_Ds_id <- first_incidence_age %>%
    mutate(Ds_id = code2id(diag_icd10)) %>%
    select(-diag_icd10) %>%
    mutate(age_diag = round(age_diag))

  # the patient_list provide the column index for efficiently breaking down matrix into list of matrices
  para$patient_lst <- para$unlist_Ds_id %>%
    mutate(id = row_number()) %>%
    select(eid, id) %>%
    group_by(eid) %>%
    group_split(keep = F) %>%
    lapply(pull)

  para$w <- para$unlist_Ds_id %>%
    group_by(eid) %>%
    group_split(keep = F)

  # this list is splited by disease
  para$ds_list <- para$unlist_Ds_id %>%
    select(-eid) %>%
    mutate(id = row_number()) %>%
    group_by(Ds_id) %>%
    group_split()

  print(object.size(para$w), unit = "MB" , standard = "SI")

  # initiate beta
  para$eta <- rgamma(para$D,shape = 100, rate = 100)

  # each column is a topic; D*K matrix
  para$beta <- t(rdirichlet(para$K, para$eta))

  # initiate alpha
  para$alpha <- rgamma(para$K, shape = 50, rate = 10)

  # this beta_w parameter save the beta for each word: it is a list of M elements and each contain a K*Ns matrix
  para$beta_w_full <- para$beta[para$unlist_Ds_id$Ds_id,,drop=FALSE]
  para$beta_w <- lapply(para$patient_lst, function(x) para$beta_w_full[x,,drop=F])

  # Matrix of M*K
  para$E_lntheta <- t(sapply(rgamma(para$M,shape = 100, rate = 100), function(x) x*(digamma(para$alpha) - digamma(sum(para$alpha))) ))

  # update E_zn: list of M; each element is matrix of Ns*K
  para <- comp_E_zn(para)

  # eventually reset alpha to be non-informative
  para$alpha <- rep(1, para$K)

  # prepare the treeLDA
  para <- create_structure(para)
  return(para)
}

# topics with Tree structures
# for point estimate of beta, we don't need to change any other updates including the lower bound

# The key is created a similar structure of the beta_w; each word is mapped to multiple layer.
create_structure <- function(para){
  # there are three metrics to compute: 1. NodeDS index matrix of (L)-by-D: elements refereing parent for disease j at layer l
  # 2. a disease list: treeDS[[l]][[cl]] the records associated with each node, tree equivalent of para$ds_list (using the list of diseases code from 1:para$D that are under that node)
  # 3. the beta list: treeBeta[[l]] is the vector (matrix to include i) of betas at this layer
  # 4. para$Cl[[l]] save the number of nodes at each layer
  # using Z for filling the empty
  filled_list <- as.character(para$list_above500occu$diag_icd10)
  # make sure death is separate from the tree
  if("Death" %in% filled_list){
    filled_list[which(filled_list == "Death")] <- "ZZZZ"}
  # healthy is also seperate from other branches
  if("Healthy" %in% filled_list){
    filled_list[which(filled_list == "Healthy")] <- "YYYY"}

  para$NodeDS <- matrix(nrow = para$L, ncol = para$D)
  para$treeDS <- list()
  para$treeBeta <- list()
  para$Cl <- list()
  for(l in 1:para$L){
    if(l == para$L){
      nodes_l <- unlist(lapply(filled_list, function(x) (substring(x, 1,l+4))))
    }else{
      nodes_l <- unlist(lapply(filled_list, function(x) (substring(x, 1,l)))) # longest ICD10 is 7 character
    }
    node_set <- unique(nodes_l)
    para$Cl[[l]] <- length(node_set) # save the number of children at each layer
    para$NodeDS[l,] <- match(nodes_l, node_set)
    para$treeDS[[l]] <- list()
    for(nd in 1:length(node_set)){
      para$treeDS[[l]][[nd]] <- bind_rows(para$ds_list[which(node_set[nd] == nodes_l)])
    }
  }
  # initialize treeBeta
  for(l in 1:para$L){
    para$treeBeta[[l]] <- matrix(nrow = para$Cl[[l]], ncol = para$K)
  }
  return(para)
}

update_beta_treeLDA <- function(para){
  # compute the beta
  for(l in 1:para$L){
    # compute M-step for beta: basic case, direct maximize the upper bound
    para$treeBeta[[l]] <- sapply(1:para$Cl[[l]], function(c) colSums(para$unlist_zn[para$treeDS[[l]][[c]]$id,]) ) %>% t
    # the normaliztion is the key: we need to normalize with respect to common parent
    if(l == 1){
      para$treeBeta[[l]] <- sapply(1:para$K, function(i) para$treeBeta[[l]][,i]/sum(para$treeBeta[[l]][,i]))
    }else{
      for(c in 1:para$Cl[[l-1]]){
        parent_idx <- which(para$NodeDS[l-1,] == c) # the indices at the parent layer
        childeren_idx <-  unique(para$NodeDS[l, parent_idx])# the child under this c
        para$treeBeta[[l]][childeren_idx, ] <- sapply(1:para$K, function(i) para$treeBeta[[l]][childeren_idx,i]/
                                                        sum(para$treeBeta[[l]][childeren_idx,i]))
      }
    }
  }

  # compute the beta at leaves: using para$NodeDS going to go down the tree
  for(j in 1:para$D){
    beta_along_tree <- sapply(1:para$L, function(l) para$treeBeta[[l]][para$NodeDS[l,j],,drop = F])
    para$beta[j,] <- sapply(1:para$K, function(i) prod(beta_along_tree[i,]))
  }

  # this beta_w parameter save the beta for each word: it is a list of M elements and each contain a K*Ns matrix
  # para$beta_w <- lapply(para$w, function(w) para$beta[w$Ds_id,,drop=FALSE] )
  para$beta_w_full <- para$beta[para$unlist_Ds_id$Ds_id,,drop=FALSE]
  para$beta_w <- lapply(para$patient_lst, function(x) para$beta_w_full[x,,drop=F])
  return(para)
}

# example
rec_data <- read.csv("DiseaseAbove500occur_include_death_ICDA2N.csv")
ds_list <- read.csv("listAbove500include_deaths_ICDA2N.csv")
topic_num <- 10
tree_leveles <- 4
para <- topic_init_tree(rec_data, ds_list, topic_num, tree_leveles)


# set the number of update
para$max_itr <- 500
para$lb <- data.frame("Iteration" = 0,"Lower_bound" = comp_lda_lb(para))
para$alpha <- rep(.1, para$K)
para$tol <- 10^(-6) # tolerance of lower bound step
para$itr_check_lb <- 5
for(itr in 1:para$max_itr){
  print(paste0("Interation: ",itr))
  para <- CVB0_E_zn(para)
  para <- update_beta_treeLDA(para)
  # print(paste0("Current Lower bound at beta ", comp_lda_lb(para)))
  # para <- update_alpha(para) # we do not need to update alpha
  # print(paste0("Current Lower bound at alpha ", comp_lda_lb(para)))
  if(mean(para$alpha) < 1){ # in this case, abandone the approximated lower bound, using the Mean-field lower bound.
    para <- comp_E_lntheta(para)
    para$lb[nrow(para$lb) + 1,] <- c(itr, comp_lda_lb(para))
  }else{
    para$lb[nrow(para$lb) + 1,] <- c(itr, CVB_lb(para))
  }
  if(itr %% para$itr_check_lb  ==0){
    print(paste0("Current Lower bound ", pull(filter(para$lb, Iteration == itr), Lower_bound), " at iteration: ",itr))
    curr_lb <- pull(filter(para$lb, Iteration == itr), Lower_bound)
    prev_lb <- pull(filter(para$lb, Iteration == (itr -para$itr_check_lb )), Lower_bound)
    if(is.finite((curr_lb - prev_lb)) & abs(curr_lb - prev_lb)/abs(prev_lb) < para$tol ){
      print(paste0("Optimization converged at step ", itr))
      break
    }
  }
}

# save the parameter
save(para, file = paste0("~/Desktop/comorbidity/Results/","Run_treeLDA_",Sys.Date(),".RData"))





