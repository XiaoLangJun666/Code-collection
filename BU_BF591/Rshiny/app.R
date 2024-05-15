## Author: Taylor Falk
## tfalk@bu.edu
## BU BF591
## Assignment 7
install.packages('colour1picker')
# Welcome to R Shiny. All that glitters is not gold.
library(shiny)
library(bslib)
library(ggplot2)
library(colourpicker) # you might need to install this
library(DT)

# Define UI for application that draws a histogram
ui <- fluidPage(
  titlePanel("BF591 Assignment 7"),
  p('To use this application, download the CSV deseq_res.csv from the data directory of this app\'s repository'),
  sidebarLayout(
    sidebarPanel(
      width=3,
      height='120vh',
      fileInput('data','Load differential expression results'),
      HTML("A volcano plot can be generated with \"<b>log<sub>2</sub> fold-change</b>\" on the x-axis and \"<b>p-adjusted</b>\" on the y-axis."),
      br(),
      br(),
      radioButtons('x_value','Choose the column for the x-axis',c('baseMean','log2FoldChange','lfcSE','stat','pvalue','padj'),selected='log2FoldChange'),
      radioButtons('y_value','Choose the column for the y-axis',c('baseMean','log2FoldChange','lfcSE','stat','pvalue','padj'),selected='padj'),
      colourInput('color1','Base point color',value='#22577A'),
      colourInput('color2','Highlight point color',value='#FFCF56'),
      sliderInput(inputId = "slider", min = -300, max = 0,
                  label = "How many pyramids are there?", value = -150, step = 1),
      actionButton("plotButton", icon = icon(name='fa car-burst', lib = "font-awesome"),"Plot",width='100%')),
    mainPanel(
      tabsetPanel(
        id = "tabs",
        tabPanel("Plot", plotOutput("volcano",width='100%',height='80vh')),
        tabPanel("Table", tableOutput("table"))
      )

      )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
    #' load_Data
    #'
    #' @details Okay this one is a little weird but bear with me here. This is 
    #' still a "function", but it will take no arguments. The `reactive({})` bit 
    #' says "if any of my inputs (as in, input$...) are changed, run me again". 
    #' This is useful when a user clicks a new button or loads a new file. In 
    #' our case, look for the uploaded file's datapath argument and load it with 
    #' read.csv. Return this data frame in the normal return() style.
    load_data <- reactive({
        # if (is.null(input$data) || is.null(input$data$filepath)) {
        #   # Return NULL if no file is uploaded
        #   return(NULL)
        # }
        return(read.csv(input$data$datapath, header = TRUE))
    })

    #' Volcano plot
    #'
    #' @param dataf The loaded data frame.
    #' @param x_name The column name to plot on the x-axis
    #' @param y_name The column name to plot on the y-axis
    #' @param slider A negative integer value representing the magnitude of
    #' p-adjusted values to color. Most of our data will be between -1 and -300.
    #' @param color1 One of the colors for the points.
    #' @param color2 The other colors for the points. Hexadecimal strings: "#CDC4B5"
    #'
    #' @return A ggplot object of a volcano plot
    #' @details I bet you're tired of these plots by now. Me too, don't worry.
    #' This is _just_ a normal function. No reactivity, no bells, no whistles.
    #' Write a normal volcano plot using geom_point, and integrate all the above
    #' values into it as shown in the example app. The testing script will treat
    #' this as a normal function.
    #'
    #' !!sym() may be required to access column names in ggplot aes().
    #'
    #' @examples volcano_plot(df, "log2fc", "padj", -100, "blue", "taupe")
    volcano_plot <-
        function(dataf, x_name, y_name, slider, color1, color2) {
            p<-ggplot(dataf,aes(.data[[x_name]],-log10(.data[[y_name]])))+
              geom_point(aes(color=ifelse(
                is.na(.data[[y_name]])==TRUE,"NA",ifelse(
                  .data[[y_name]]<10^(slider),"TRUE","FALSE")
              )
              )
              )+
              scale_color_manual(values=c("TRUE"=color2,"FALSE"=color1,'NA'='grey'),name=paste0(y_name,'<1*10^',slider))+
              theme_bw()+
              theme(legend.position='bottom')+
              labs(x=x_name,y=paste0('-log10(',y_name,')'))
            return(p)
        }

    #' Draw and filter table
    #'
    #' @param dataf Data frame loaded by load_data()
    #' @param slider Negative number, typically from the slider input.
    #'
    #' @return Data frame filtered to p-adjusted values that are less than
    #' 1 * 10^slider, columns for p-value and p-adjusted value have more digits
    #' displayed.
    #' @details Same as above, this function is a standard R function. Tests will
    #' evaluate it normally. Not only does this function filter the data frame to
    #' rows that are above the slider magnitude, it should also change the format
    #' of the p-value columns to display more digits. This is so that it looks
    #' better when displayed on the web page. I would suggest the function
    #' `formatC()`
    #'formatC()
    #' @examples draw_table(deseq_df, -210)
    #'    X  baseMean     log2FC     lfcSE      stat       pvalue         padj
    #'gene1 11690.780   9.852926 0.2644650  37.25607 8.45125e-304 1.54472e-299
    #'gene2  3550.435  -6.183714 0.1792708 -34.49369 9.97262e-261 9.11398e-257
    draw_table <- function(dataf, slider) {
        dataf_sorted<-dataf[dataf$pvalue<10^(slider),]
        dataf_sorted$pvalue<-formatC(dataf_sorted$pvalue,digits=6)
        dataf_sorted$padj<-formatC(dataf_sorted$padj,digits=6)
        names(dataf_sorted)[names(dataf_sorted)=='X']<-'gene'
        
        return(dataf_sorted)
    }

    #' These outputs aren't really functions, so they don't get a full skeleton,
    #' but use the renderPlot() and renderTabel() functions to return() a plot
    #' or table object, and those will be displayed in your application.
      
    output$volcano <- renderPlot({
      if(input$plotButton){
        volcano_plot(load_data(),input$x_value,input$y_value,input$slider,input$color1,input$color2)

        }
    })

    # Same here, just return the table as you want to see it in the web page
    # output$table <- renderDT({datatable(
    #   head(draw_table(load_data(),input$slider),10),
    #   option=list(
    #   autoWidth=TRUE,
    #   pageLength=10,
    #   lengthChange=FALSE,
    #   searching=FALSE,
    #   ordering=FALSE,
    #   info=FALSE
    # ),
    #   class='cell-border stripe')},) # replace this NULL
    output$table<-renderTable({draw_table(load_data(),input$slider)},striped=TRUE,width='100%',spacing='m')
    
}

# Run the application
shinyApp(ui = ui, server = server)
