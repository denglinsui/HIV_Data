#' This file includes the codes for data visualization
#'

plot_data <- function(data, criteria, fdr=NULL){
  data <- 
    data %>% 
    filter( Criteria %in% criteria) 
  data <- melt(data,  id=c("Criteria","drug","Layer"))
  data <- reshape2::dcast(data, drug+variable+Layer~Criteria,
                          value.var = "value")
  p <- 
    ggplot(data, aes(x=get(criteria[1]), y =get(criteria[2]), 
                     color = variable))+
    geom_point()+
    xlab(criteria[1])+
    ylab(criteria[2])+
    facet_wrap(.~Layer)
  
  if(!is.null(fdr)){
    p <- p + geom_hline(aes(yintercept=fdr), linetype = 2)
  }
  
  return(p)
}

plot_data_class <- function(data, criteria, fdr=NULL){
  data <- 
    data %>% 
    filter( Criteria %in% criteria) 
  data <- melt(data,  id=c("Criteria","drug","Layer","DrugClass"))
  data <- reshape2::dcast(data, drug+variable+Layer+DrugClass~Criteria,
                          value.var = "value")
  p <- 
    ggplot(data, aes(x=get(criteria[1]), y =get(criteria[2]), 
                     color = variable))+
    geom_point()+
    xlab(criteria[1])+
    ylab(criteria[2])+
    facet_grid(Layer~DrugClass)
  
  if(!is.null(fdr)){
    p <- p + geom_hline(aes(yintercept=fdr), linetype = 2)
  }
  
  return(p)
}
