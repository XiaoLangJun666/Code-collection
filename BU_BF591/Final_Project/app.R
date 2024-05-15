library(shiny)
library(DT)
library(ggplot2)
library(GEOquery)
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(ggbeeswarm)
library(colourpicker)
library(fgsea)
options(shiny.maxRequestSize = 100*1024^2)
# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("BF591 Final Project Rshiny App"),
  br(),
  tabsetPanel(
    id = "main_tabs",
    tabPanel("Samples",
      sidebarLayout(
        sidebarPanel(
          width=2,
          fileInput('meta_data_matrix','Please Load GSE Meta Data Matrix!')
        ),
        mainPanel(
          width=10,
          tabsetPanel(
            id='second_tabs_samples',
            tabPanel("Summary",
                     tableOutput("meta_data_summary_table")),
            tabPanel('Table',
                     dataTableOutput("meta_data_table")),
            tabPanel('Plots',
                     sidebarLayout(
                       sidebarPanel(
                         width=2,
                         radioButtons('y_value_meta_data_plot','Choose the column for the y-axis',c('Age of onset','Age of death','h-v cortical score','h-v strialtal score','mRNA read counts','Duration'),selected='Age of death'),
                         radioButtons('meta_data_plot_type','Choose the plot type',c('Histogram plot','Density plot','Violin plot'),selected='Violin plot'),
                         actionButton("meta_data_plotButton","Plot",width='100%')
                       ),
                       mainPanel(
                         width=10,
                         div(style = "margin-left: 200px;", 
                             plotOutput("meta_plot", width='80%', height='60vh')
                         )
                       )
                     )))
          )
    )),
    tabPanel("Counts",
             sidebarLayout(
               sidebarPanel(
                 width=2,
                 br(),
                 
                 fileInput('data_count_matrix','Please Load Data Count Matrix!'),
                 br(),
                 br(),
                 sliderInput("slider_count_percentile_variance", 
                             "Choose the threshold for genes with at least X percentile of variance", 
                             min = 0, 
                             max = 100, 
                             value = 95, 
                             step = 1),
                 br(),
                 uiOutput("slider_count_nonzero_sample"),
                 br(),
                 actionButton("Count_Analysis","Analysis",width='100%')
               ),
               mainPanel(
                 width=10,
                 tabsetPanel(
                   id='second_tabs_counts',
                   tabPanel('Summary',
                            tableOutput("count_data_summary_table")),
                   tabPanel('Plots',
                            plotOutput('count_data_scatter_plot_1'),
                            plotOutput('count_data_scatter_plot_2')),
                   tabPanel('Heatmap',
                            plotOutput('count_data_heatmap_plot',width='100%',height='80vh')),
                   tabPanel("PCA",
                            sidebarLayout(
                              sidebarPanel(
                                width=2,
                                uiOutput('pca_pc_select_1'),
                                br(),
                                br(),
                                uiOutput('pca_pc_select_2'),
                                br(),
                                br(),
                                radioButtons('pca_plot_type','Choose the plot type:',c('Scatter Plot','Beeswarm Plot')),
                                br(),
                                uiOutput('pca_beeswarm_number'),
                                actionButton('PCA_Plot','Plot',width='100%')
                              ),
                              mainPanel(
                                width=8,
                                plotOutput('count_data_pca_plot',width='100%',height='60vh')
                              )
                            )
                            )
                 )
                 
               )
             )),
    tabPanel('DE',
             sidebarLayout(
               sidebarPanel(
                 width=3,
                 height='120vh',
                 fileInput('DE_data','Load differential expression results'),
                 HTML("A volcano plot can be generated with \"<b>log<sub>2</sub> fold-change</b>\" on the x-axis and \"<b>p-adjusted</b>\" on the y-axis."),
                 br(),
                 br(),
                 radioButtons('DE_x_value','Choose the column for the x-axis',c('baseMean','log2FoldChange','lfcSE','stat','pvalue','padj'),selected='log2FoldChange'),
                 uiOutput('DE_slider_1'),
                 radioButtons('DE_y_value','Choose the column for the y-axis',c('baseMean','log2FoldChange','lfcSE','stat','pvalue','padj'),selected='padj'),
                 uiOutput('DE_slider_2'),
                 colourInput('DE_color1','Base point color',value='#22577A'),
                 colourInput('DE_color2','Highlight point color',value='#FFCF56'),
                 actionButton("DE_plotButton", "Plot",width='100%')),
               mainPanel(
                 tabsetPanel(
                   id = "tabs",
                   tabPanel("Table", dataTableOutput("DE_table")),
                   tabPanel("Plot", plotOutput("DE_volcano",width='100%',height='80vh')),
                   tabPanel("Filtered Table",dataTableOutput('DE_table_filtered'))
                 )
                 
               )
             )),
    tabPanel('GSEA',
             sidebarLayout(
               sidebarPanel(
                 width=2,
                 radioButtons('fgsea_file_value','Choose the type of data type:',c('Differential Analysis Result','fgsea Result'),selected='Differential Analysis Result'),
                 uiOutput('fgsea_data_input'),
                 uiOutput('fgsea_database'),
                 uiOutput('fgsea_threshold_min'),
                 uiOutput('fgsea_threshold_max'),
                 actionButton("fgsea_analysis_start",'Analysis',width='100%')
                 # br(),
                 # br(),
                 # uiOutput('fgsea_data_button')
               ),
               mainPanel(
                 width=10,
                 tabsetPanel(
                   id="GSEA_tabset",
                   tabPanel('Top Results',
                            sidebarLayout(
                              sidebarPanel(
                                width=3,
                                uiOutput('fgsea_tab1_slider')
                              ),
                              mainPanel(
                                width=9,
                                plotOutput("fgsea_bar_plot",width='100%',height='80vh')
                              )
                            )),
                   tabPanel('Table',
                            sidebarLayout(
                              sidebarPanel(
                                width=3,
                                uiOutput('fgsea_table_slider'),
                                radioButtons('fgsea_NES_value','Choose the NES value that you want to display',c('Positive','Negative','All'),selected='All')
                              ),
                              mainPanel(
                                width=9,
                                dataTableOutput("fgsea_tab2_table")
                              )
                            )),
                   tabPanel('Plots',
                            sidebarLayout(
                              sidebarPanel(
                                width=3,
                                uiOutput('fgsea_scatter_plot_slider'),
                              ),
                              mainPanel(
                                plotOutput('fgsea_scatter_plot',width='100%',height='80vh')
                              )
                            ))
                 )
               )
             ))
  )

)

# Define server logic required to draw a histogram
server <- function(input, output) {
  #First part for meta data information
  load_meta_data_matrix <- reactive({
    file<-input$meta_data_matrix$datapath
    
    if(is.null(file))
      return(NULL)
    
    return(getGEO(filename=input$meta_data_matrix$datapath,GSEMatrix=TRUE))
  })
  
  draw_meta_data_summary_table<-function(meta_data){
    df <- data.frame(column_name=character(), type=character(), 'Mean (sd) or Distinct Values'=numeric(), stringsAsFactors=FALSE,check.names=FALSE)
    df <- rbind(df, data.frame(column_name='Sample number', type=class(length(meta_data$title)), 'Mean (sd) or Distinct Values'=length(meta_data$title), stringsAsFactors=FALSE,check.names=FALSE))
    df <- rbind(df, data.frame(column_name='GEO accession', type=class(meta_data$geo_accession[0]), 'Mean (sd) or Distinct Values'=paste0(meta_data$geo_accession[1],'-',meta_data$geo_accession[length(meta_data$geo_accession)]), stringsAsFactors=FALSE,check.names=FALSE))
    df <- rbind(df, data.frame(column_name='Pulic Date', type=class(unique(meta_data$status)), 'Mean (sd) or Distinct Values'=unique(meta_data$status), stringsAsFactors=FALSE,check.names=FALSE))
    df <- rbind(df, data.frame(column_name='Source', type=class(unique(meta_data$source_name_ch1)), 'Mean (sd) or Distinct Values'=unique(meta_data$source_name_ch1), stringsAsFactors=FALSE,check.names=FALSE))
    df <- rbind(df, data.frame(column_name='Organism', type=class(unique(meta_data$organism_ch1)), 'Mean (sd) or Distinct Values'=unique(meta_data$organism_ch1), stringsAsFactors=FALSE,check.names=FALSE))
    df <- rbind(df, data.frame(column_name='Diagnosis', type='Character', 'Mean (sd) or Distinct Values'=paste(unique(meta_data$`diagnosis:ch1`), collapse = ", "), stringsAsFactors=FALSE,check.names=FALSE))
    df <- rbind(df, data.frame(column_name='Age of onset', type=class(as.numeric(meta_data$`age of onset:ch1`)), 'Mean (sd) or Distinct Values'=paste0(round(mean(as.numeric(meta_data$`age of onset:ch1`),na.rm=TRUE),2),'(+/-',round(sd(as.numeric(meta_data$`age of onset:ch1`),na.rm=TRUE),2),')'), stringsAsFactors=FALSE,check.names=FALSE))
    df <- rbind(df, data.frame(column_name='Age of death', type=class(as.numeric(meta_data$`age of death:ch1`)), 'Mean (sd) or Distinct Values'=paste0(round(mean(as.numeric(meta_data$`age of death:ch1`),na.rm=TRUE),2),'(+/-',round(sd(as.numeric(meta_data$`age of death:ch1`),na.rm=TRUE),2),')'), stringsAsFactors=FALSE,check.names=FALSE))
    df <- rbind(df, data.frame(column_name='Duration', type=class(as.numeric(meta_data$`duration:ch1`)), 'Mean (sd) or Distinct Values'=paste0(round(mean(as.numeric(meta_data$`duration:ch1`),na.rm=TRUE),2),'(+/-',round(sd(as.numeric(meta_data$`duration:ch1`),na.rm=TRUE),2),')'), stringsAsFactors=FALSE,check.names=FALSE))
    df <- rbind(df, data.frame(column_name='h-v cortical score', type=class(as.numeric(meta_data$`h-v cortical score:ch1`)), 'Mean (sd) or Distinct Values'=paste0(round(mean(as.numeric(meta_data$`h-v cortical score:ch1`),na.rm=TRUE),2),'(+/-',round(sd(as.numeric(meta_data$`h-v cortical score:ch1`),na.rm=TRUE),2),')'), stringsAsFactors=FALSE,check.names=FALSE))
    df <- rbind(df, data.frame(column_name='h-v striatal score', type=class(as.numeric(meta_data$`h-v striatal score:ch1`)), 'Mean (sd) or Distinct Values'=paste0(round(mean(as.numeric(meta_data$`h-v striatal score:ch1`),na.rm=TRUE),2),'(+/-',round(sd(as.numeric(meta_data$`h-v striatal score:ch1`),na.rm=TRUE),2),')'), stringsAsFactors=FALSE,check.names=FALSE))
    df <- rbind(df, data.frame(column_name='mRNA read counts', type=class(as.numeric(meta_data$`mrna-seq reads:ch1`)), 'Mean (sd) or Distinct Values'=paste0(round(mean(as.numeric(meta_data$`mrna-seq reads:ch1`),na.rm=TRUE),2),' (+/-',round(sd(as.numeric(meta_data$`mrna-seq reads:ch1`),na.rm=TRUE),2),')'), stringsAsFactors=FALSE,check.names=FALSE))
    df <- rbind(df, data.frame(column_name='Vonsattel grade', type='Character', 'Mean (sd) or Distinct Values'=paste(unique(meta_data$`vonsattel grade:ch1`), collapse = ", "), stringsAsFactors=FALSE,check.names=FALSE))
    
    return(df)
  }
  
  draw_meta_data_table<-function(meta_data){
    return(pData(phenoData(meta_data)))
  }
  
  draw_meta_data_plot<-function(meta_data,y_value,plot_type){
    if (y_value=='Age of onset'){
      y_axis = as.numeric(as.character(meta_data$`age of onset:ch1`))
    }else if(y_value=='Age of death'){
      y_axis = as.numeric(as.character(meta_data$`age of death:ch1`))
    }else if(y_value=='mRNA read counts'){
      y_axis = as.numeric(as.character(meta_data$`mrna-seq reads:ch1`))
    }else if(y_value=='Duration'){
      y_axis = as.numeric(as.character(meta_data$`duration:ch1`))
    }else if(y_value=='h-v cortical score'){
      y_axis = as.numeric(as.character(meta_data$`h-v cortical score:ch1`))
    }else if(y_value=='h-v strialtal score'){
      y_axis = as.numeric(as.character(meta_data$`h-v striatal score:ch1`))
    }
    
    plot_data <- data.frame(
      diagnosis = as.factor(meta_data$`diagnosis:ch1`),
      y_axis = y_axis
    )
    plot_data <- na.omit(plot_data)
    
    if (plot_type=='Histogram plot'){
      meta_p<-ggplot(plot_data, aes(x = y_axis, fill = diagnosis)) +
        geom_histogram(bins = 10, color = "black") + 
        labs(title = "Histogram Plot",
             x = y_value,
             y = "Frequency") +
        theme_minimal()
    }else if(plot_type=='Density plot'){
      meta_p<-ggplot(plot_data, aes(x = y_axis, fill = diagnosis)) +
        geom_density(alpha = 0.5) + 
        scale_fill_brewer(palette = "Set1") + 
        labs(title = "Density Plot",
             x = y_value,
             y = "Density") +
        theme_minimal()
    }else{
      meta_p<-ggplot(plot_data, aes(x = diagnosis, y = y_axis,fill=diagnosis)) +
        geom_violin(trim = FALSE) +
        theme_minimal() +
        labs(title = "Violin Plot",
             x = "Diagnosis", 
             y = y_value)
    }
    return (meta_p)
  }
  
  output$meta_data_summary_table<-renderTable({draw_meta_data_summary_table(load_meta_data_matrix())},striped=TRUE,width='100%')
  
  output$meta_data_table<-renderDT({
    datatable(draw_meta_data_table(load_meta_data_matrix()),options = list(scrollX = TRUE))
  })
  
  output$meta_plot <- renderPlot({
    if(input$meta_data_plotButton){
      draw_meta_data_plot(load_meta_data_matrix(),input$y_value_meta_data_plot,input$meta_data_plot_type)
      
    }
  })
  
  #Second part for count data
  load_count_data_matrix<-reactive({
    file<-input$data_count_matrix$datapath
    
    if(is.null(file))
      return(NULL)
    count_matrix<-read.csv(input$data_count_matrix$datapath,sep='\t')
    
    return(list(data=count_matrix,ncol=ncol(count_matrix)-1))
  })
  
  
  output$slider_count_nonzero_sample<-renderUI({
    max_samples<-load_count_data_matrix()$ncol
    
    if(is.null(max_samples))
      return()
    
    sliderInput("slider_count_nonzero_sample", 
                "Choose the threshold for genes with at least X samples that are non-zero", 
                min = 0, 
                max = max_samples, 
                value =60, 
                step = 1)
  })
  
  
  count_data_filter<-function(count_data,threshold_variance,threshold_sample){
    gene_variance <- apply(count_data[-1],1,var)
    selected_var<-quantile(gene_variance,probs=threshold_variance/100)
    count_data_new<-count_data[apply(count_data[-1],1,var)>=selected_var,]
    count_data_new<-count_data_new[apply(count_data_new[-1],1,function(x) sum(x!=0))>threshold_sample,]
    
    df <- data.frame('Summary Item'=character(), Calculation=numeric() , stringsAsFactors=FALSE,check.names=FALSE)
    df <- rbind(df, data.frame('Summary Item'='Number of samples', Calculation=ncol(count_data)-1 , stringsAsFactors=FALSE,check.names=FALSE))
    df <- rbind(df, data.frame('Summary Item'='Total number of genes', Calculation=nrow(count_data) , stringsAsFactors=FALSE,check.names=FALSE))
    df <- rbind(df, data.frame('Summary Item'='Genes passing current filter', Calculation=paste0(nrow(count_data_new),' / ',nrow(count_data_new)/nrow(count_data),'%') , stringsAsFactors=FALSE,check.names=FALSE))
    df <- rbind(df, data.frame('Summary Item'='Genes not passing current filter', Calculation=paste0(nrow(count_data)-nrow(count_data_new),' / ',(nrow(count_data)-nrow(count_data_new))/nrow(count_data),'%'), stringsAsFactors=FALSE,check.names=FALSE))
    
    count_data_new[-1]<-lapply(names(count_data_new[-1]),function(x) count_data_new[-1][[x]]*10^6/vapply(count_data_new[-1],sum,1)[[x]])
    
    return(list(summary_df=df,filtered_data=count_data_new))
  }
  
  
  filtered_data_summary<-reactive({
    return(count_data_filter(load_count_data_matrix()$data,input$slider_count_percentile_variance,input$slider_count_nonzero_sample))
  })
  
  draw_count_data_scatter_plot<-function(count_data,filtered_data){

    scatter_data <- data.frame(gene_name = count_data)
    scatter_data$median_count<-apply(count_data[-1],1,median)
    scatter_data$variance <- apply(count_data[-1], 1, var)
    scatter_data$num_zeros <- apply(count_data[-1], 1, function(x) sum(x == 0))
    scatter_data$filtered <- count_data[[1]] %in% filtered_data[[1]]

    p1 <- ggplot(scatter_data, aes(x = median_count, y = variance, color = filtered)) +
      geom_point(alpha=0.5,size=0.5) +
      scale_color_manual(values = c("lightgray","darkred")) +
      scale_y_log10() +
      scale_x_log10()+
      ggtitle("Median Count vs Variance")

    p2 <- ggplot(scatter_data, aes(x = num_zeros, y =median_count , color = filtered)) +
      geom_point(alpha=0.5,size=1) +
      scale_color_manual(values = c("lightgray","darkgreen")) +
      scale_y_log10() +
      ggtitle("Median Count vs Number of Zeros")

    return(list(plot1=p1,plot2=p2))
  }
  
  filtered_data_scatter_plot<-reactive({
    return(draw_count_data_scatter_plot(load_count_data_matrix()$data,filtered_data_summary()$filtered_data))
  })
  
  heatmap_plot_count_data<-function(filtered_data){
    log_transformed_counts <- log1p(filtered_data[-1])
    p<-pheatmap(log_transformed_counts, 
             color = colorRampPalette(c("blue", "white", "red"))(50),
             scale = "row",
             clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean",
             clustering_method = "complete",
             show_rownames = TRUE,
             show_colnames = TRUE,
             legend = TRUE)
    return(p)
  }
  
  
  pca_Analysis<-function(filtered_data){
    
    pca_results <- prcomp(scale(t(filtered_data[-1])), center=FALSE, scale=TRUE)
    
    PCA_tibble<-as_tibble(pca_results$x)%>%
      mutate(Sample=ifelse(substr(rownames(pca_results$x),1,1)=='H','Huntington',ifelse(substr(rownames(pca_results$x),1,1)=='C','Control',NA)))
    
    prop_var <- pca_results$sdev^2 / sum(pca_results$sdev^2) * 100
    
    return(list(var=prop_var,pca_tibble=PCA_tibble,pca_number=ncol(PCA_tibble)-1))
  }
  
  output$pca_pc_select_1<-renderUI({
    pca_n<-pca_Analysis(filtered_data_summary()$filtered_data)$pca_number
    
    if(is.null(pca_n))
      return()
    
    selectInput("pcNumber_1", "Choose The First Principal Component You Want to Plot", 
                choices = 1:pca_n)
  })
  
  output$pca_pc_select_2<-renderUI({
    pca_n<-pca_Analysis(filtered_data_summary()$filtered_data)$pca_number
    
    if(is.null(pca_n))
      return()
    
    selectInput("pcNumber_2", "Choose The Second Principal Component You Want to Plot", 
                 choices = 1:pca_n)
  })
  
  
  output$pca_beeswarm_number<-renderUI({
    pca_n<-pca_Analysis(filtered_data_summary()$filtered_data)$pca_number
    if(is.null(pca_n))
      return()
    
    if(input$pca_plot_type=='Beeswarm Plot'){
      selectInput('Top_N_Beeswarm',"Choose The Top N PCs You Want to Plot In Beeswarm Plot",
                  choices=1:pca_n)
    }
  })
  
  pca_scatter_plot<-function(pp_var,pca_data,pc1,pc2){
    pc_1<-paste0('PC',pc1)
    pc_2<-paste0('PC',pc2)

    p<-ggplot(pca_data, aes(x = .data[[pc_1]], y = .data[[pc_2]], color = Sample)) +
      geom_point() +
      labs(x = paste0(pc_1,': ', round(pp_var[as.numeric(pc1)]), "% variance"),
           y = paste0(pc_2,': ', round(pp_var[as.numeric(pc2)]), "% variance"))+
      theme_classic()
    
    return(p)

  }
  
  pca_beeswarm_plot<-function(pca_data,n){
    pc_n<-paste0("PC",n)
    p<-pca_data%>%
      pivot_longer(PC1:pc_n,names_to="PC",values_to="projection")%>%
      mutate(PC=fct_relevel(PC,str_c("PC",1:n))) %>%
      ggplot(aes(x=PC,y=projection,color=Sample)) +
      geom_beeswarm() + labs(title="PCA Projection Plot") +
      theme(axis.text.x=element_text(angle=90,hjust=0,vjust=0.5))
    return(p)
  }
  
  output$count_data_summary_table<-renderTable({
    if(input$Count_Analysis){
      filtered_data_summary()$summary_df
    }
    },striped=TRUE,width='100%')
  
  output$count_data_scatter_plot_1<-renderPlot({
    if(input$Count_Analysis){
      filtered_data_scatter_plot()$plot1
    }})
  
  output$count_data_scatter_plot_2<-renderPlot({
    if(input$Count_Analysis){
      filtered_data_scatter_plot()$plot2
    } })
  output$count_data_heatmap_plot<-renderPlot({
      heatmap_plot_count_data(filtered_data_summary()$filtered_data)
  })
  
  output$count_data_pca_plot<-renderPlot({
    if(input$PCA_Plot){
      if (input$pca_plot_type=='Scatter Plot'){
        pca_scatter_plot(pca_Analysis(filtered_data_summary()$filtered_data)$var,pca_Analysis(filtered_data_summary()$filtered_data)$pca_tibble,input$pcNumber_1,input$pcNumber_2)
      }else{
        pca_beeswarm_plot(pca_Analysis(filtered_data_summary()$filtered_data)$pca_tibble,input$Top_N_Beeswarm)
      }
    }
  })
  
  
  ##DE part
  load_DE_data <- reactive({
    return(read.csv(input$DE_data$datapath,sep='\t', header = TRUE))
  })
  
  output$DE_slider_1<-renderUI({
    input_1<-input$DE_x_value
    
    sliderInput("DE_slider_x_value", 
                paste0("Choose the threshold for ",input_1,'.'), 
                min = round(min(load_DE_data()[[input_1]]))-1, 
                max = round(max(load_DE_data()[[input_1]]))+1, 
                value =round(max(load_DE_data()[[input_1]])/2), 
                step = 0.1)
  })
  
  output$DE_slider_2<-renderUI({
    input_2<-input$DE_y_value
    
    sliderInput("DE_slider_y_value", 
                paste0("Choose the threshold for ",input_2,'.'), 
                min = round(min(-log10(load_DE_data()[[input_2]])))-1, 
                max = round(max(-log10(load_DE_data()[[input_2]])))+1, 
                value =round(max(-log10(load_DE_data()[[input_2]]))/2), 
                step = 1)
  })
  
  DE_volcano_plot <-
    function(dataf, x_name, y_name, slider1,slider2, color1, color2) {
      p<-ggplot(dataf,aes(.data[[x_name]],-log10(.data[[y_name]])))+
        geom_point(aes(color=ifelse(
          is.na(.data[[y_name]])==TRUE,"NA",ifelse(
            .data[[y_name]]<10^(-slider1) & abs(.data[[x_name]])>=abs(slider2),"TRUE","FALSE")
        )
        )
        )+
        geom_vline(xintercept = c(-slider2, slider2), linetype = "dashed", color = "blue") + 
        geom_hline(yintercept = slider1, linetype = "dashed", color = "red") + 
        scale_color_manual(values=c("TRUE"=color2,"FALSE"=color1,'NA'='grey'),name=paste0(y_name,'<1*10^',slider1))+
        theme_bw()+
        theme(legend.position='bottom')+
        labs(x=x_name,y=paste0('-log10(',y_name,')'))
      return(p)
    }
  
  draw_DE_filtered_table <- function(dataf, slider1,slider2,x_name,y_name) {
    dataf_sorted<-dataf[dataf[[y_name]]<10^-(slider2),]
    dataf_sorted<-dataf_sorted[abs(dataf_sorted[[x_name]])>abs(slider1),]
    dataf_sorted$x_name<-formatC(dataf_sorted[[x_name]],digits=6)
    dataf_sorted$y_name<-formatC(dataf_sorted[[y_name]],digits=6)
    
    return(dataf_sorted)
  }
  
  
  output$DE_volcano <- renderPlot({
    if(input$DE_plotButton){
      DE_volcano_plot(load_DE_data(),input$DE_x_value,input$DE_y_value,input$DE_slider_y_value,input$DE_slider_x_value,input$DE_color1,input$DE_color2)
      
    }
  })
  
  
  output$DE_table<-DT::renderDataTable(DT::datatable(load_DE_data(),extensions='Buttons',class='display',
                                                     options=list(paging=TRUE,searching=TRUE,
                                                                  fixedColumns=TRUE,autoWidth=TRUE,
                                                                  ordering=TRUE,dom='Bfrtip',
                                                                  butttons=c('copy','csv'))))
  
  output$DE_table_filtered<-DT::renderDataTable(DT::datatable(draw_DE_filtered_table(load_DE_data(),input$DE_slider_x_value,input$DE_slider_y_value,input$DE_x_value,input$DE_y_value),extensions='Buttons',class='display',
                                                     options=list(paging=TRUE,searching=TRUE,
                                                                  fixedColumns=TRUE,autoWidth=TRUE,
                                                                  ordering=TRUE,dom='Bfrtip',
                                                                  butttons=c('copy','csv'))))
  ## fgsea
  
  output$fgsea_data_input<-renderUI({
    if (input$fgsea_file_value=='Differential Analysis Result'){
      fileInput('fgsea_data_fgsea','Load differential expression results')
    }else{
      fileInput('fgsea_data_fgsea','Load fgsea analysis results')
    }
  })
  
  # output$fgsea_data_button<-renderUI({
  #   if(input$fgsea_file_value=='Differential Analysis Result'){
  #       downloadButton('download_fgsea_data','Download fgsea data as CSV')
  #   } 
  # })
  
  load_fgsea_data <- reactive({
    return(read.csv(input$fgsea_data_fgsea$datapath,sep='\t', header = TRUE))
  })
  

  output$fgsea_database<-renderUI({
    if(input$fgsea_file_value=='Differential Analysis Result'){
      data_files<-list.files(path='msigdb_v2023.2.Hs_GMTs',full.names=FALSE)
      selectInput('fgsea_database_file','Choose a gmt file that you want to analysis:', choices=data_files)
    }
    
  })
  
  output$fgsea_threshold_min<-renderUI({
    if(input$fgsea_file_value=='Differential Analysis Result'){
      numericInput('min_num','Enter the min number for fgsea analysis:',value=15)
    }
  })
  
  output$fgsea_threshold_max<-renderUI({
    if(input$fgsea_file_value=='Differential Analysis Result'){
      numericInput('max_num','Enter the max number for fgsea analysis:',value=500)
    }
  })
  
  
  fgsea_analysis<-function(rank_data,dataset,min,max){
    rank_list<-as_tibble(rank_data)%>%
      arrange(log2FoldChange)%>%
      select(symbol,log2FoldChange)
  
    rnk_list<-setNames(rank_list$log2FoldChange,rank_list$symbol)
    
    pathways_fgsea <- fgsea::gmtPathways(paste0('msigdb_v2023.2.Hs_GMTs/',dataset))
    
    fgsea_results <- fgsea(pathways_fgsea, rnk_list, minSize =min, maxSize=max)%>%
                             as_tibble()%>%
                             arrange(padj)
    
    return(fgsea_results)
  }
  
  # output$download_fgsea_results<-downloadHandler(
  #   filename=function(){
  #     paste0('fgsea_results_',Sys.Date(),'.csv')
  #   },
  #   content=function(file){
  #     write.csv(fgsea_analysis(load_fgsea_data(),input$fgsea_database_file,input$min_num,input$max_num),file,row.names=FALSE)
  #   }
  # )
  fgsea_analysis_results<-reactive({
    if (input$fgsea_analysis_start){
      return(fgsea_analysis(load_fgsea_data(),input$fgsea_database_file,input$min_num,input$max_num))
    }
  })
  
  
  output$fgsea_tab1_slider<-renderUI({
    sliderInput('fgsea_tab1_padj',
                'Choose the top pathways to plot by padj',
                min=0,
                max=40,
                value=20,
                step=1)
  })
  
  
  draw_fgsea_bar_plot<-function(data,top_n){
    data<-data[order(data$padj,decreasing=FALSE),]
    data<-data[1:top_n,]
    data$color<-ifelse(data$NES>0,"Positive","Negative")
    data$short_pathway <- ifelse(nchar(data$pathway) > 80, 
                                    substr(data$pathway, 1, 80), 
                                    data$pathway)
    # top_NES<-top_NES%>%
    #   arrange(desc(NES))
    p<-ggplot(data,aes(x=NES,y = reorder(short_pathway, NES),fill=color))+
      geom_bar(stat='identity')+
      scale_fill_manual(values = c("Positive" = "red", "Negative" = "blue"))+
      theme_minimal() +
      labs(x='Normalized Enrichment Score (NES)')
    return(p)
  }
  
  output$fgsea_bar_plot<-renderPlot({
    draw_fgsea_bar_plot(fgsea_analysis_results(),input$fgsea_tab1_padj)
  })
  
  
  output$fgsea_table_slider<-renderUI({
    sliderInput('fgsea_tab2_padj',
                'Choose the threshold for adjusted p-value:',
                min=round(min(-log10(fgsea_analysis_results()$padj))-1),
                max=round(max(-log10(fgsea_analysis_results()$padj))+1),
                value=2,
                step=1)
  })
  
  draw_fgsea_padj_filter_table<-function(data,p_th, scope){
    data<-data%>%
      filter(padj<10^-(p_th))
    if (scope=='Positive'){
      data<-data%>%
        filter(NES>0)
    }else if(scope=='Negative'){
      data<-data%>%
        filter(NES<0)
    }
    
    return(data)
  }
  
  output$fgsea_tab2_table <- DT::renderDataTable({
    DT::datatable(
      draw_fgsea_padj_filter_table(fgsea_analysis_results(), input$fgsea_tab2_padj, input$fgsea_NES_value),
      extensions = 'Buttons',
      options = list(
        paging = FALSE, searching = TRUE,
        fixedColumns = TRUE, autoWidth = FALSE,
        ordering = TRUE, dom = 'Bfrtip',
        buttons = list(
          list(
            extend = 'copy',
            exportOptions = list(
              modifier = list(page = 'all')  # 
            )
          ),
          list(
            extend = 'csv',
            exportOptions = list(
              modifier = list(page = 'all',
                              search='none')  
            )
          )
        )
      )
    )
  })
  
  output$fgsea_scatter_plot_slider<-renderUI({
    sliderInput('fgsea_tab3_padj',
                'Choose the threshold for adjusted p-value:',
                min=round(min(-log10(fgsea_analysis_results()$padj))-1),
                max=round(max(-log10(fgsea_analysis_results()$padj))+1),
                value=2,
                step=1)
  })
  
  
  draw_fgsea_scatter_plot<-function(data,p_th){
    data$log10_padj<--log10(data$padj)
    
    p<-ggplot(data,aes(x=NES,y=log10_padj))+
      geom_point(aes(color=log10_padj>p_th),size=1)+
      scale_color_manual(values=c("TRUE"='red',"FALSE"='grey'))+
      labs(x = 'NES (Normalized Enrichment Score)', y = '-log10 Adjusted p-value') +
      theme_minimal()
    return(p)
  }
  
  output$fgsea_scatter_plot<-renderPlot({
    draw_fgsea_scatter_plot(fgsea_analysis_results(),input$fgsea_tab3_padj)
  })
}
  
# Run the application 
shinyApp(ui = ui, server = server)
