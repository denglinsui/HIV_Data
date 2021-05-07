#' Run the procedures to detect the gene assosiated HIV-1 drug resistence
#'
 
Run_Procedures <- function (X, y, q) {
  # Log-transform the drug resistance measurements.
  y = log(y)
  
  # Remove patients with missing measurements.
  missing = is.na(y)
  y = y[!missing]
  X = X[!missing,]
  
  # Remove predictors that appear less than 3 times.
  X = X[,colSums(X) >= 5]
  
  # Remove duplicate predictors.
  X = X[,colSums(abs(cor(X)-1) < 1e-4) == 1]
  
  # Run the knockoff filter.
  knock.gen = function(x) create.fixed(x, method='equi')
  result = knockoff.filter(X, y, fdr=q, 
                           knockoffs=knock.gen, 
                           statistic=stat.glmnet_lambdasmax)
  knockoff_selected = names(result$selected)
  
  # Run the multilayer filter
  # MULTILAYER KNOCKOFF FILTER OPTIONS
  fdr.multi = c(q,q)
  knockoff_type = "fixed_equi"
  statistic_type = "group_LSM" 
  FDP_hat_type = "kn+"
  groups.tmp <- get_position(colnames(X))
  all.gen <- unique(groups.tmp)
  groups <- matrix(sapply(groups.tmp, 
                          function(x){which(all.gen==x)}),
                   ncol=1)
  groups <- cbind(1:nrow(groups),groups)
  
  # RUN MULTILAYER KNOCKOFF FILTER
  output = multilayer_knockoff_filter(X, y, groups, fdr.multi,
                                      knockoff_type,
                                      statistic_type, FDP_hat_type)
  knockoff_selected_multi <- colnames(X)[output$S_hat]
  
  # Run the knockoff filter(Second Order).
  knock.gen = function(x) create.second_order(x, method='equi')
  result = knockoff.filter(X, y, fdr=q, knockoffs=knock.gen, statistic=stat.glmnet_lambdasmax)
  #knockoff_selected_2ndO = names(result$selected)
  
  # Run BHq.
  p = ncol(X)
  lm.fit = lm(y ~ X - 1) # no intercept
  p.values = coef(summary(lm.fit))[,4]
  names(p.values) <- substring(names(p.values),2)
  cutoff = max(c(0, which(sort(p.values) <= q * (1:p) / p)))
  bhq_selected = names(which(p.values <= q * cutoff / p))
  
  # Run p-filter
  p_filter = pfilter(p.values, fdr.multi, groups)
  PF =  colnames(X)[which(p_filter==1)]
  # PF Result
  #PF = PF(knockoff_selected, PF(knockoff_selected_2ndO, bhq_selected))
  list(Knockoff = knockoff_selected,
       MKF = knockoff_selected_multi,
       BHq = bhq_selected, 
       PF = PF)
}

#' Flatten a matrix to a vector with names from concatenating row/column names.
flatten_matrix <- function(M, sep='.') {
  x <- c(M)
  names(x) <- c(outer(rownames(M), colnames(M),
                      function(...) paste(..., sep=sep)))
  x
}

#' Flatten the design matrix X
Flatten_X <- function(Drug){
  gene_df <- Drug$gene_df
  muts = c(LETTERS, 'i', 'd')
  pos_start = which(names(gene_df) == 'P1')
  
  X = outer(muts, as.matrix(gene_df[,pos_start:ncol(gene_df)]),
            Vectorize(grepl))
  X = aperm(X, c(2,3,1))
  dimnames(X)[[3]] <- muts
  X = t(apply(X, 1, flatten_matrix))
  mode(X) <- 'numeric'
  return(X)
}

RunProcedures.combine <- function (X, Y, q) {
  # Combine Y together
  n <- nrow(X)
  p <- ncol(X)
  num.task <- ncol(Y)
  drug.type <- rep(names(Y), each = nrow(X))
  X.new <- matrix(0, ncol = ncol(X)*num.task, nrow = num.task*nrow(X))
  for(i in 1:num.task){
    ind.x <- (1+(i-1)*n) : (i*n )
    ind.y <- (1+(i-1)*p) : (i*p )
    X.new[ind.x, ind.y] <- X
  }
  colnames(X.new) <- rep(colnames(X),times = num.task)
  
  y <- as.numeric(unlist(Y))
  
  # Log-transform the drug resistance measurements.
  y = log(y)
  
  # Remove patients with missing measurements.
  missing = is.na(y)
  y = y[!missing]
  X.new = X.new[!missing,]
  drug.type = drug.type[!missing]
  
  # Remove predictors that appear less than 3 times.
  X.new.ind <- colSums(X.new) >= 3
  X.new = X.new[,X.new.ind]
  
  # Remove duplicate predictors.
  X.new.ind <- colSums(abs(cor(X.new)-1) < 1e-4) == 1
  X.new = X.new[,X.new.ind]
  
  gene.groups.tmp <- colnames(X.new)
  gene.all.gen <- unique(gene.groups.tmp)
  gene.group <- matrix(sapply(gene.groups.tmp,
                              function(x){which(gene.all.gen==x)}),
                       ncol=1)
  # Run the group knockoff.
  
  # Run the multilayer filter
  fdr.multi = c(q,q)
  knockoff_type = "fixed_equi"
  statistic_type = "group_LSM" 
  FDP_hat_type = "kn+"
  
  # group_knockoff_single pos+direction
  
  output.group = multilayer_knockoff_filter(X.new, y, gene.group, c(q),
                                            knockoff_type,
                                            statistic_type, FDP_hat_type)
  knockoff_selected.group <- colnames(X.new)[output.group$S_hat]
  
  # MULTILAYER KNOCKOFF FILTER OPTIONS
  groups.tmp <- get_position(colnames(X.new))
  all.gen <- unique(groups.tmp)
  groups <- matrix(sapply(groups.tmp, 
                          function(x){which(all.gen==x)}),
                   ncol=1)
  groups <- cbind(gene.group,groups)
  
  # RUN MULTILAYER KNOCKOFF FILTER
  output = multilayer_knockoff_filter(X.new, y, groups, fdr.multi,
                                      knockoff_type,
                                      statistic_type, FDP_hat_type)
  knockoff_selected_MKF.group <- colnames(X.new)[output$S_hat]
  
  # Run BHq.
  p = ncol(X.new)
  lm.fit = lm(y ~ X.new - 1) # no intercept
  p.values = coef(summary(lm.fit))[,4]
  names(p.values) <- substring(names(p.values),6)
  cutoff = max(c(0, which(sort(p.values) <= q * (1:p) / p)))
  bhq_selected = names(which(p.values <= q * cutoff / p))
  
  # Run p-filter
  BH.group <- pfilter(p.values, q, gene.group)
  BH.selected <- colnames(X.new)[which(BH.group==1)]
  
  # Run p-filter
  p_filter = pfilter(p.values, fdr.multi, groups)
  PF =  colnames(X.new)[which(p_filter==1)]
  # PF Result
  #PF = PF(knockoff_selected, PF(knockoff_selected_2ndO, bhq_selected))
  list(Knockoff = knockoff_selected.group,
       MKF = knockoff_selected_MKF.group,
       BHq = BH.selected, 
       bhq = bhq_selected,
       PF = PF)
}


RunProcedures_3layer.combine <- function (X, Y, q) {
  # Combine Y together
  n <- nrow(X)
  p <- ncol(X)
  num.task <- ncol(Y)
  drug.type <- rep(names(Y), each = nrow(X))
  X.new <- matrix(0, ncol = ncol(X)*num.task, nrow = num.task*nrow(X))
  for(i in 1:num.task){
    ind.x <- (1+(i-1)*n) : (i*n )
    ind.y <- (1+(i-1)*p) : (i*p )
    X.new[ind.x, ind.y] <- X
  }
  colnames(X.new) <- rep(colnames(X),times = num.task)
  
  y <- as.numeric(unlist(Y))
  
  # Log-transform the drug resistance measurements.
  y = log(y)
  
  # Remove patients with missing measurements.
  missing = is.na(y)
  y = y[!missing]
  X.new = X.new[!missing,]
  drug.type = drug.type[!missing]
  
  # Remove predictors that appear less than 3 times.
  X.new.ind <- colSums(X.new) >= 3
  X.new = X.new[,X.new.ind]
  
  # Remove duplicate predictors.
  X.new.ind <- colSums(abs(cor(X.new)-1) < 1e-4) == 1
  X.new = X.new[,X.new.ind]
  
  gene.groups.tmp <- colnames(X.new)
  gene.all.gen <- unique(gene.groups.tmp)
  gene.group <- matrix(sapply(gene.groups.tmp,
                              function(x){which(gene.all.gen==x)}),
                       ncol=1)
  # Run the group knockoff.
  
  # Run the multilayer filter
  fdr.multi = c(q,q,q)
  knockoff_type = "fixed_equi"
  statistic_type = "group_LSM" 
  FDP_hat_type = "kn+"
  
  # group_knockoff_single pos+direction
  
  output.group = multilayer_knockoff_filter(X.new, y, gene.group, c(q),
                                            knockoff_type,
                                            statistic_type, FDP_hat_type)
  knockoff_selected.group <- colnames(X.new)[output.group$S_hat]
  
  # MULTILAYER KNOCKOFF FILTER OPTIONS
  groups.tmp <- get_position(colnames(X.new))
  all.gen <- unique(groups.tmp)
  groups <- matrix(sapply(groups.tmp, 
                          function(x){which(all.gen==x)}),
                   ncol=1)
  groups <- cbind(c(1:nrow(groups)),gene.group,groups)
  
  # RUN MULTILAYER KNOCKOFF FILTER
  output = multilayer_knockoff_filter(X.new, y, groups, fdr.multi,
                                      knockoff_type,
                                      statistic_type, FDP_hat_type)
  knockoff_selected_MKF.group <- colnames(X.new)[output$S_hat]
  
  # Run BHq.
  p = ncol(X.new)
  lm.fit = lm(y ~ X.new - 1) # no intercept
  p.values = coef(summary(lm.fit))[,4]
  names(p.values) <- substring(names(p.values),6)
  cutoff = max(c(0, which(sort(p.values) <= q * (1:p) / p)))
  bhq_selected = names(which(p.values <= q * cutoff / p))
  
  # Run p-filter
  BH.group <- pfilter(p.values, q, gene.group)
  BH.selected <- colnames(X.new)[which(BH.group==1)]
  
  # Run p-filter
  p_filter = pfilter(p.values, fdr.multi, groups)
  PF =  colnames(X.new)[which(p_filter==1)]
  # PF Result
  #PF = PF(knockoff_selected, PF(knockoff_selected_2ndO, bhq_selected))
  list(Knockoff = knockoff_selected.group,
       MKF = knockoff_selected_MKF.group,
       BHq = BH.selected, 
       bhq = bhq_selected,
       PF = PF)
}


