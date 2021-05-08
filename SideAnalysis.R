source('~/project/HIV_Data/codes/LoadingFiles.R')

fdr = 0.20
drug_class <- "PI" 
res <- get(paste0("res_",drug_class))
posmut <- get(paste0("tsm_",drug_class))

Criteria <- c("dp","fdp")

# Define Datasets for PIs
X <- Flatten_X(res)
gene_df <- res$gene_df
tsm_df <- res$tsm_df
Y = res$Y

# Run all algorithms (Run drug types in given drug class jointly: 3Layer)
res_3l.cmb = RunProcedures_3layer.combine(X, Y, fdr) 

cmp_pos_3l.cmb <- comp_pos_cmb(res_3l.cmb, tsm_df)
data_pos_3l.tmp <-  do.call(cbind, cmp_pos_3l.cmb)
data_pos_3l.cmb <- plot.combine(data_pos_3l.tmp, drug_class)

cmp_posmut_3l.cmb <- comp_posmut_cmb(res_3l.cmb, posmut)
data_posmut_3l.tmp <-  do.call(cbind, cmp_posmut_3l.cmb)
data_posmut_3l.cmb <- plot.combine(data_posmut_3l.tmp, drug_class)

data_3l.cmb <- rbind(data_pos_3l.cmb %>% 
                       mutate(Layer="Layer: Position"),
                     data_posmut_3l.cmb %>%
                       mutate(Layer="Layer: Position+Direction"))
p_3l.cmb <- plot_data(data_3l.cmb, Criteria, fdr)
p_3l.cmb
ggsave(paste0("Figure/3_Layer",drug_class,".eps"), width = 12, height = 6)

# Investigate the relationship between quality of discoveries and the number of mutation 
X_PI <- Flatten_X(res_PI)
res_PI.cmb <- RunProcedures.combine(X_PI, 
                                   res_PI$Y, 
                                   fdr) 
X_NRTI <- Flatten_X(res_NRTI)
res_NRTI.cmb <- RunProcedures.combine(X_NRTI, 
                                     res_NRTI$Y,
                                     fdr) 

X_NNRTI <- Flatten_X(res_NNRTI)
res_NNRTI.cmb <- RunProcedures.combine(X_NNRTI, 
                                      res_NNRTI$Y, 
                                      fdr) 
PF_PI.cmb <- unique(res_PI.cmb$PF)
PF_NRTI.cmb <- unique(res_NRTI.cmb$PF)
PF_NNRTI.cmb <- unique(res_NNRTI.cmb$PF)

X.PI.count <- X_PI %>% 
  as.data.frame() %>%
  select(PF_PI.cmb) %>%
  colSums() %>%
  sort()
X.PI.judge <- names(X.PI.count) %in% tsm_PI
  
X.NRTI.count <- X_NRTI %>% 
  as.data.frame() %>%
  select(PF_NRTI.cmb) %>%
  colSums()%>%
  sort()
X.NRTI.judge <- names(X.NRTI.count) %in% tsm_NRTI

X.NNRTI.count <- X_NNRTI %>% 
  as.data.frame() %>%
  select(PF_NNRTI.cmb) %>%
  colSums()%>%
  sort()
X.NNRTI.judge <- names(X.NNRTI.count) %in% tsm_NNRTI

X.PI.judge <- names(X.PI.count) %in% tsm_PI
freq_data <- data.frame(name = c(names(X.PI.count),
                                 names(X.NRTI.count),
                                 names(X.NNRTI.count)),
                        count = c(X.PI.count,
                                  X.NRTI.count,
                                  X.NNRTI.count),
                        Positive = c(X.PI.judge,
                                     X.NRTI.judge,
                                     X.NNRTI.judge),
                        Class = rep(c("PI","NRTI","NNRTI"),
                                    c(length(PF_PI.cmb),
                                      length(PF_NRTI.cmb),
                                      length(PF_NNRTI.cmb))))

freq_p <- 
  ggplot(freq_data, aes(x=count, fill = Positive))+
  geom_density(alpha=0.5)+
  facet_grid(.~Class)
freq_p


freq_gt_data <- data.frame(name = c(tsm_PI,
                                 tsm_NRTI,
                                 tsm_NNRTI),
                        count = c(colSums(X_PI[,tsm_PI]),
                                  colSums(X_NRTI[,tsm_NRTI]),
                                  colSums(X_NNRTI[,tsm_NNRTI])),
                        Positive = "Ground_Truth",
                        Class = rep(c("PI","NRTI","NNRTI"),
                                    c(length(tsm_PI),
                                      length(tsm_NRTI),
                                      length(tsm_NNRTI))))

freq_com_data <- rbind(freq_data, freq_gt_data)
freq_p <- 
  ggplot(freq_com_data, aes(x=count, fill = Positive))+
  geom_density(alpha=0.5)+
  facet_grid(.~Class)
freq_p

freq_pos_data <- 
  freq_data %>%
  mutate(name = as.character(get_position(name))) %>%
  group_by(Class,name) %>%
  summarise(count = sum(count), 
            Positive = as.factor(as.logical(sum(Positive))))
freq_gt_pos_data <- 
  freq_gt_data %>%
  mutate(name = as.character(get_position(name))) %>%
  group_by(Class,name) %>%
  summarise(count = sum(count), 
            Positive = Positive)

freq_com_pos_data <- rbind(freq_pos_data,freq_gt_pos_data)
freq_com_pos_p <- 
  ggplot(freq_com_pos_data, aes(x=count, fill = Positive))+
  geom_density(alpha=0.5)+
  facet_grid(.~Class)
freq_com_pos_p

freq_final_data <- rbind(freq_com_data %>%
                           select(names(freq_com_pos_data)) %>%
                           as.data.frame(),
                         freq_com_pos_data)
freq_final_data <- data.frame(Class = levels(freq_com_data$Class)[c(freq_com_data$Class, freq_com_pos_data$Class)],
                              name = c(freq_com_data$name, freq_com_pos_data$name),
                              Positive = as.factor(c(freq_com_data$Positive, 
                                                     levels(freq_com_pos_data$Positive)[freq_com_pos_data$Positive])),
                              count = c(freq_com_data$count, freq_com_pos_data$count),
                              Position_Type = rep(c("PosMut","Pos"),
                                                  c(nrow(freq_com_data),nrow(freq_com_pos_data))))
levels(freq_final_data$Positive) <- c("False Detect","Ground Truth","True Detect") 

freq_final_p <- 
  ggplot(freq_final_data, aes(x=count, color = Positive))+
  geom_density(alpha=0.3)+
  facet_grid(Position_Type~Class)
freq_final_p

ggsave("Figure/Freq_Figure.eps", width = 12, height = 6)

# Combine the result of three drug classes
load("TmpResult/ResProcedeure.RData")
data.cmb <- NULL
for(drug_class in c("PI","NRTI","NNRTI")){
  res.cmb <- get(paste0("res.com_",drug_class))
  res <- get(paste0("res_",drug_class))
  tsm_df <- res$tsm_df
  posmut <- get(paste0("tsm_",drug_class))
  cmp_pos.cmb <- comp_pos_cmb(res.cmb, tsm_df)
  
  data_pos.tmp <-  do.call(cbind, cmp_pos.cmb)
  data_pos.cmb <- plot.combine(data_pos.tmp, drug_class)
  
  cmp_posmut.cmb <- comp_posmut_cmb(res.cmb, posmut)
  data_posmut.tmp <-  do.call(cbind, cmp_posmut.cmb)
  data_posmut.cmb <- plot.combine(data_posmut.tmp, drug_class)
  
  data.cmb.tmp <- rbind(data_pos.cmb %>% 
                      mutate(Layer="Layer: Position"),
                    data_posmut.cmb %>%
                      mutate(Layer="Layer: Position+Direction"))
  data.cmb.tmp$DrugClass <- rep(drug_class, nrow(data.cmb.tmp))
  
  data.cmb <- rbind(data.cmb, data.cmb.tmp)
}

p.cmb_class <- plot_data_class(data.cmb, Criteria, fdr)
p.cmb_class

ggsave(paste0("Figure/Combine.eps"), width = 12, height = 6)

data.sep <- NULL
for(drug_class in c("PI","NRTI","NNRTI")){
  res.sep <- get(paste0("res.sep_",drug_class))
  res <- get(paste0("res_",drug_class))
  tsm_df <- res$tsm_df
  posmut <- get(paste0("tsm_",drug_class))
  cmp_pos.cmb <- comp_pos_cmb(res.cmb, tsm_df)
  
  cmp_pos <- comp_pos(res.sep)
  data_pos <- plot_pos(cmp_pos)
  
  cmp_posmut <- comp_posmut(res.sep, posmut)
  data_posmut <- plot_posmut(cmp_posmut)
  
  data.sep.tmp <- rbind(data_pos %>% 
                      mutate(Layer="Layer: Position"),
                    data_posmut %>%
                      mutate(Layer="Layer: Position+Direction"))
  
  data.sep.tmp$DrugClass <- rep(drug_class, nrow(data.cmb.tmp))
  
  data.sep <- rbind(data.sep, data.sep.tmp)
}

p.sep_class <- plot_data_class(data.sep, Criteria, fdr)
p.sep_class
ggsave(paste0("Figure/Separate.eps"), width = 12, height = 6)

# Show the selected gene for each case

for(drug_class in c("PI","NRTI","NNRTI")){
  res.cmb <- get(paste0("res.com_",drug_class))
  res <- get(paste0("res_",drug_class))
  tsm_df <- res$tsm_df
  posmut <- get(paste0("tsm_",drug_class))
  
  Intersect_Gene <- 
    res.cmb$Knockoff %>%
    intersect(res.cmb$BHq) %>%
    intersect(res.cmb$bhq)%>%
    intersect(res.cmb$PF) %>%
    sort()
  cat("drug class:", drug_class, "\n")
  Intersect_Positive <- Intersect_Gene[Intersect_Gene %in% posmut]
  Intersect_Negative <- Intersect_Gene[!(Intersect_Gene %in% posmut)]
  cat("PositiveGene:",Intersect_Positive,"\n",
      "NegativeGene:",Intersect_Negative,"\n")
  Intersect_Pos <- unique(get_position(Intersect_Gene))
  Intersect_Pos_Positive <- Intersect_Pos[Intersect_Pos %in% tsm_df$Position]
  Intersect_Pos_Negative <- Intersect_Pos[!(Intersect_Pos %in% tsm_df$Position)]
  
  cat("PositiveGenePos:",Intersect_Pos_Positive,"\n",
      "NegativeGenePos:",Intersect_Pos_Negative,"\n")
}
