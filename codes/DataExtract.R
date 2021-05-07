#' Extract data from http://hivdb.stanford.edu
#' @param res_PI,res_NRTI,res_NNRTI Essential Information for each drug class;
#' @param tsm_PI,tsm_NRTI,tsm_NNRTI The ground truth 

drug_class_all <- c("PI","NRTI","NNRTI")
Drug_extract <- function(drug_class){
  base_url = 'http://hivdb.stanford.edu/pages/published_analysis/genophenoPNAS2006'
  gene_url = paste(base_url, 'DATA', paste0(drug_class, '_DATA.txt'), sep='/')
  tsm_url = paste(base_url, 'MUTATIONLISTS', 'NP_TSM', drug_class, sep='/')
  
  gene_df = read.delim(gene_url, na.string = c('NA', ''), stringsAsFactors = FALSE)
  
  tsm_df = read.delim(tsm_url, header = FALSE, stringsAsFactors = FALSE)
  names(tsm_df) = c('Position', 'Mutations')
  
  pos_start = which(names(gene_df) == 'P1')
  Drug_name = paste(names(gene_df)[4:(pos_start-1)],collapse = " ")
  Y = gene_df[,4:(pos_start-1)]

  return(list(gene_df = gene_df,
              tsm_df = tsm_df,
              Drug_name = Drug_name,
              Drug_class = drug_class,
              Y = Y))
}

res_PI <- Drug_extract("PI")
res_NRTI <- Drug_extract("NRTI")
res_NNRTI <- Drug_extract("NNRTI")

Drug_summary <- data.frame(DrugClass = drug_class_all,
                           Drug = c(res_PI$Drug_name,
                                    res_NRTI$Drug_name,
                                    res_NNRTI$Drug_name))
row.names(Drug_summary) <- NULL

Extract_tsm <- function(tsm_df){
  PosMut <- NULL
  for(i in 1:nrow(tsm_df)){
    Pos <- tsm_df$Position[i]
    Mut <- strsplit(tsm_df$Mutations[i], split = ' ')[[1]]
    PosMut <- c(PosMut, paste0("P",Pos,".",Mut))
  }  
  return(PosMut)
}

tsm_PI <- Extract_tsm(res_PI$tsm_df)
tsm_NRTI <- Extract_tsm(res_NRTI$tsm_df)
tsm_NNRTI <- Extract_tsm(res_NNRTI$tsm_df)

