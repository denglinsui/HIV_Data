source('~/project/HIV_Data/codes/LoadingFiles.R')

fdr <- 0.20
drug_class <- "NRTI" 
res <- get(paste0("res_",drug_class))
posmut <- get(paste0("tsm_",drug_class))

Criteria <- c("dp","fdp")

# Define Datasets for PIs
X <- Flatten_X(res)
gene_df <- res$gene_df
tsm_df <- res$tsm_df
Y = res$Y

# Run all algorithms (Run drug types in given drug class separately)
res.sep <- lapply(Y, function(y) Run_Procedures(X, y, fdr))

cmp_pos <- comp_pos(res.sep)
data_pos <- plot_pos(cmp_pos)

cmp_posmut <- comp_posmut(res.sep, posmut)
data_posmut <- plot_posmut(cmp_posmut)

data.sep <- rbind(data_pos %>% 
                    mutate(Layer="Layer: Position"),
                  data_posmut %>%
                    mutate(Layer="Layer: Position+Direction"))
p.sep <- plot_data(data.sep, Criteria, fdr)
p.sep
ggsave(paste0("Figure/Sep",drug_class,".eps"), width = 12, height = 6)

# Run all algorithms (Run drug types in given drug class jointly)
res.cmb <- RunProcedures.combine(X, Y, fdr) 

cmp_pos.cmb <- comp_pos_cmb(res.cmb, tsm_df)
data_pos.tmp <-  do.call(cbind, cmp_pos.cmb)
data_pos.cmb <- plot.combine(data_pos.tmp, drug_class)

cmp_posmut.cmb <- comp_posmut_cmb(res.cmb, posmut)
data_posmut.tmp <-  do.call(cbind, cmp_posmut.cmb)
data_posmut.cmb <- plot.combine(data_posmut.tmp, drug_class)

data.cmb <- rbind(data_pos.cmb %>% 
                    mutate(Layer="Layer: Position"),
                  data_posmut.cmb %>%
                    mutate(Layer="Layer: Position+Direction"))
p.cmb <- plot_data(data.cmb, Criteria, fdr)
p.cmb
ggsave(paste0("Figure/Cmb",drug_class,".eps"), width = 12, height = 6)

assign(paste0("res.sep_",drug_class),res.sep)
assign(paste0("res.com_",drug_class),res.cmb)
save(res.com_PI, res.com_NRTI, res.com_NNRTI,
     res.sep_PI, res.sep_NRTI, res.sep_NNRTI,
     tsm_PI, tsm_NRTI, tsm_NNRTI,
     tsm_df_PI, 
     tsm_df_NRTI, 
     tsm_df_NNRTI, 
     file="TmpResult/ResProcedeure.RData")


#=== Show the intersect of selected set of all algorithm
Intersect_Gene <- 
  res.cmb$Knockoff %>%
  intersect(res.cmb$BHq) %>%
  intersect(res.cmb$bhq)%>%
  intersect(res.cmb$PF) %>%
  sort()

Intersect_Positive <- Intersect_Gene[Intersect_Gene %in% posmut]
Intersect_Negative <- Intersect_Gene[!(Intersect_Gene %in% posmut)]

Intersect_Pos <- unique(get_position(Intersect_Gene))
Intersect_Pos_Positive <- Intersect_Pos[Intersect_Pos %in% tsm_df$Position]
Intersect_Pos_Negative <- Intersect_Pos[!(Intersect_Pos %in% tsm_df$Position)]

cat("drug class:", drug_class, "\n")
cat("PositiveGene:",Intersect_Positive,"\n",
    "NegativeGene:",Intersect_Negative,"\n")
cat("PositiveGenePos:",Intersect_Pos_Positive,"\n",
    "NegativeGenePos:",Intersect_Pos_Negative,"\n")
