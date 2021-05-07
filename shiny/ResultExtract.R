#' This file includes the codes to extract information in result.

#' Get the position of gene
get_position <- function(x)
  sapply(regmatches(x, regexpr('[[:digit:]]+', x)), as.numeric)

#' Extract the position of result
comp_pos <- function(results){
  comparisons_pos <- lapply(results, function(drug) {
    lapply(drug, function(selected) {
      positions = unique(get_position(selected)) # remove possible duplicates
      discoveries = length(positions)
      false_discoveries = length(setdiff(positions, tsm_df$Position))
      list(true_discoveries = discoveries - false_discoveries,
           false_discoveries = false_discoveries,
           fdp = false_discoveries / max(1, discoveries),
           dp = (discoveries - false_discoveries) / length(tsm_df$Position))
    })
  })
  
  return(comparisons_pos)
}

#' Extract the information to plot result of position
plot_pos <- function(comparisons_pos, tsm_df){
  plot_data_pos <- NULL
  for (drug in names(comparisons_pos)) {
    plot_data <- do.call(cbind, comparisons_pos[[drug]])
    
    plot_data_pos <- rbind(plot_data_pos,
                           data.frame(KF = 
                                        as.numeric(plot_data[,'Knockoff']),
                                      MKF =
                                        as.numeric(plot_data[,'MKF']),
                                      BH = as.numeric(plot_data[,'BHq']),
                                      PF = as.numeric(plot_data[,'PF']),
                                      Criteria = rownames(plot_data),
                                      drug = rep(drug,nrow(plot_data))))
  }
  return(plot_data_pos)
}


#' Extract the position of posmut
comp_posmut <- function(results, posmut){
  comparisons_posmut <- lapply(results, function(drug) {
    lapply(drug, function(selected) {
      selected <- unique(selected)
      discoveries = length(selected)
      false_discoveries = length(setdiff(selected, posmut))
      list(true_discoveries = discoveries - false_discoveries,
           false_discoveries = false_discoveries,
           fdp = false_discoveries / max(1, discoveries),
           dp = (discoveries - false_discoveries) / length(posmut))
    })
  })
  
  return(comparisons_posmut)
}

#' Extract the information to plot result of posmut
plot_posmut <- function(comparisons_posmut){
  plot_data_posmut <- NULL
  for (drug in names(comparisons_posmut)) {
    plot_data <- do.call(cbind, comparisons_posmut[[drug]])
    
    plot_data_posmut <- rbind(plot_data_posmut,
                              data.frame(KF = 
                                           as.numeric(plot_data[,'Knockoff']),
                                         MKF =
                                           as.numeric(plot_data[,'MKF']),
                                         BH = as.numeric(plot_data[,'BHq']),
                                         PF = as.numeric(plot_data[,'PF']),
                                         Criteria = rownames(plot_data),
                                         drug = rep(drug,nrow(plot_data))))
  }
  
  return(plot_data_posmut)
}

#' Extract the position of result in joint detection
comp_pos_cmb <- function(result.combine, tsm_df){
  comparisons_pos.combine <- 
    lapply(result.combine, function(selected) {
      positions = unique(get_position(selected)) # remove possible duplicates
      discoveries = length(positions)
      false_discoveries = length(setdiff(positions, tsm_df$Position))
      list(true_discoveries = discoveries - false_discoveries,
           false_discoveries = false_discoveries,
           fdp = false_discoveries / max(1, discoveries),
           dp = (discoveries - false_discoveries) / length(tsm_df$Position))
    })
  
  return(comparisons_pos.combine)
}

#' Extract the posmut of result in joint detection
comp_posmut_cmb <- function(result.combine, posmut){
  comparisons_posmut.combine <- 
    lapply(result.combine, function(selected) {
      selected <- unique(selected)
      discoveries = length(selected)
      false_discoveries = length(setdiff(selected, posmut))
      list(true_discoveries = discoveries - false_discoveries,
           false_discoveries = false_discoveries,
           fdp = false_discoveries / max(1, discoveries),
           dp = (discoveries - false_discoveries) / length(posmut))
    })
  
  return(comparisons_posmut.combine)
}

#' Extract the information to plot result in joint detection
plot.combine <- function(plot_data.tmp, drug_class){
  plot_data.combine <- data.frame(
    GKF =  
      as.numeric(plot_data.tmp[,'Knockoff']),
    MKF = as.numeric(plot_data.tmp[,'MKF']),
    GBH = as.numeric(plot_data.tmp[,'BHq']),
    BH = as.numeric(plot_data.tmp[,'bhq']),
    PF = as.numeric(plot_data.tmp[,'PF']),
    Criteria = rownames(plot_data.tmp),
    drug = rep(drug_class,nrow(plot_data.tmp)))
  
  return(plot_data.combine)
}