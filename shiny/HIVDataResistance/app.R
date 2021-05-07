#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)


# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Application title
    titlePanel("HIV Data Drug Resistance"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            selectInput("DrugClass", "Choose a drug class", 
                        choices = c("PI", "NRTI", "NNRTI")),
            selectInput("Case", "Will you treat different drug types within this drug class separately or integrantly:", 
                        choices = c("Separate","Combine"))
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            textOutput("plot_exlain"),
            plotOutput("ResultPlot"),
            textOutput("gene_select"),
            verbatimTextOutput("Select")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    output$ResultPlot <- renderPlot({
        library(ggplot2)
        library(dplyr)
        library(reshape2)
        fdr <- 0.2
        Criteria <- c("dp","fdp")
        drug_class <- input$DrugClass
        load("ResProcedeure.RData")
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
        
        #' Get the position of gene
        get_position <- function(x){
            sapply(regmatches(x, regexpr('[[:digit:]]+', x)), as.numeric)}
        
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
        
        posmut <- get(paste0("tsm_",drug_class))
        tsm_df <- get(paste0("tsm_df_",drug_class))
        
        if(input$Case == "Separate"){
            res.sep <- get(paste0("res.sep_",drug_class))
            
            cmp_pos <- comp_pos(res.sep)
            data_pos <- plot_pos(cmp_pos)
            
            cmp_posmut <- comp_posmut(res.sep, posmut)
            data_posmut <- plot_posmut(cmp_posmut)
            
            data.sep <- rbind(data_pos %>% 
                                  mutate(Layer="Layer: Position"),
                              data_posmut %>%
                                  mutate(Layer="Layer: Position+Direction"))
            p <- plot_data(data.sep, Criteria, fdr)
            p
        }else{
            res.cmb <- get(paste0("res.com_",drug_class))
            
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
            p <- plot_data(data.cmb, Criteria, fdr)
            p
        }
        
        p
    })
    output$Select <- renderPrint({ 
        library(dplyr)
        load("ResProcedeure.RData")
        #' Get the position of gene
        get_position <- function(x){
            sapply(regmatches(x, regexpr('[[:digit:]]+', x)), as.numeric)}
        
        drug_class <- input$DrugClass
        posmut <- get(paste0("tsm_",drug_class))
        tsm_df <- get(paste0("tsm_df_",drug_class))
        res.cmb <- get(paste0("res.com_",drug_class))
        
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
        
        Print_res <- 
            rbind(paste0(Intersect_Positive,collapse = " "),
              paste0(Intersect_Negative,collapse = " "),
              paste0(Intersect_Pos_Positive,collapse = " "),
              paste0(Intersect_Pos_Negative,collapse = " "))
        rownames(Print_res) <- c("PositiveGene","NegativeGene",
                                 "PositiveGenePos","NegativeGenePos")
        Print_res
    })
    output$plot_exlain <- renderText({ 
        "We run these procedures with significance level 0.20 and this figure shows the false discovery rate and discovery rate." 
    })
    output$gene_select <- renderText({ 
        "The gene position and gene mutation selected by Knockoff, BH procedure and p-filter together." 
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
