#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(seqinr)
library(tidyverse)
library(cowplot)
library(DT)
library(Biostrings)
library(gtools)
library(fst)

#`%then%` <- shiny:::`%OR%`

source("functions.R", local = T)
options(warn=-1)


# Define UI for application that draws a histogram
ui <- fluidPage(
  tagList(
    singleton(
      tags$head(
        tags$style(type="text/css", ".dataTables_filter {display: none;    }")
      )
    )
  ),
   # Application title
   titlePanel(windowTitle="iPASS",title=div(img(src="Logo_iPASS.png", height='150'))),
   
   # Sidebar 
   sidebarLayout(
      sidebarPanel(
         fileInput("uploaded_file",
                  "File Upload",
                
                  placeholder = "Upload a fasta file"),
         withTags({hr()}),
         textAreaInput("seq_input", "Sequence Input", placeholder = "paste a DNA sequence", height = 180),
         actionButton("submit", "Submit"),
         checkboxInput("checkbox", label = "Optimize context", value = FALSE),
         withTags({
           div(class="header", checked=NA, style="text-align: justify;text-justify: inter-word;",
               hr(),
               h4("User guide:"),
              ul(
                li(b('"File upload"'),': Please upload a fasta file containing only the open reading frame to be analyzed. Only DNA base characters are allowed (G, A, T, C). Note that the first amber stop codon context to be analyzed by iPASS is at the third codon of the open reading frame (nucleotide context: -6 to -1 and +4 to +9, with TAG at +1, +2, +3).'),
                li(b('"Sequence input"'),': Please paste the open reading frame to be analyzed and use the "Submit" button to upload your sequence. Only DNA base characters are allowed (G, A, T, C). Note that the first amber stop codon context to be analyzed by iPASS is at the third codon of the open reading frame (nucleotide context: -6 to -1 and +4 to +9, with TAG at +1, +2, +3).'),
                li(b('"Optimize context"'),': Check this box to optimize the iPASS score and hence ncAA incorporation efficiency by synonymous exchange of codons flanking the amber stop codon.')
                ),
                h4('Output:'),
                p('An iPASS score of ≥ 1 should indicate above average relative ncAA incorporation efficiency. A minimal iPASS score difference of ca. 2.5 after amber stop codon context optimization usually improves ncAA incorporation efficiency. Please note that the iPASS tool has been developed to predict and optimize the suppression of amber stop codons in mammalian cells using the orthogonal',i('Methanosarcina mazei'), 'pyrrolysyl-tRNA synthetase/tRNA',sup('Pyl'),sub('CUA'),'pair. Hence, iPASS might not reliably predict ncAA incorporation efficiencies for other genetic code expansion strategies (e.g. tyrosyl amber suppressor tRNAs) or eukaryotic organisms (e.g. yeast).'),
                h4('Error report and feedback:'),
                p('Please report errors and bugs to s.bultmann(at)posteo.de. We would also appreciate to receive feedback on how well the tool predicts relative ncAA incorporation efficiencies in your applications and cell lines.'),
                h4("How to cite:"),
                p('The iPASS tool was developed by your colleagues. Please cite the publication in which iPASS has been described in your Material and Methods section: Bartoschek MD, Ugur E, Nguyen TA, Rodschinka G, Wierer M, Lang K, Bultmann S. Identification of permissive amber suppression sites for efficient non-canonical amino acid incorporation in mammalian cells. Nucleic Acids Res. 2021 Jun 21;49(11):e62. doi: 10.1093/nar/gkab132. PMID: 33684219; PMCID: PMC8216290.')
                
           )
         })
      ),
      
      # Show a plot 
      mainPanel(
        
         plotOutput("seqPlot",
                    dblclick = "seqPlot_dblclick",
                    click = clickOpts(
                      id = "seqPlot_plot_click",
                      clip = TRUE),
                    brush = brushOpts(
                      id = "seqPlot_brush",
                      resetOnNew = TRUE, 
                      direction = "x",
                      clip=TRUE,
                      delayType = "debounce"
                    )),
         withTags({div(class="header", checked=NA, style="text-align: center",
           p('Click and drag to zoom in. Double click to zoom out.'))
         }),
         dataTableOutput("info")
      )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  values <- reactiveValues(sequence = "")
  ranges <- reactiveValues(x = c(0,1e7))
  
  proxy = dataTableProxy('info')
  
  observeEvent(input$seqPlot_brush, {
    brush <- input$seqPlot_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
    } else {
      ranges$x <- c(0,1e7)
    }
  })
  observeEvent(input$seqPlot_dblclick, {
    
    ranges$x <- c(0,1e7)
    
  })
  
  validate_input <- function(file_input, text_input) {
    if (is.null(file_input) & text_input == "") {
      "Please upload or enter a sequence!"
    } else if ( text_input != "" & gsub('[ATGC]',"", toupper(text_input)) != ""){
      "Only DNA sequences are allowed. Sequence contains non base characters (ACGT)!"
    } else if ( text_input != "" & nchar(text_input)<15){
      "Please enter a sequence ≥ 15 bp !"
    } else if (!is.null(file_input) & !check_fasta(file_input$datapath)){
      "File is not in fasta format!"
    } else if (!is.null(file_input)){
      if (file.info(file_input$datapath)$size>10000000){
        "File is too large (>10MB)!"
      }
    }
    else {
      NULL
    }
  }
 
  
  
  resultCond <- reactive({
    #print(validate_input(input$uploaded_file, values$sequence))
    #print(check_fasta(input$uploaded_file$datapath))
    #print(input$uploaded_file)
    validate(
      validate_input(input$uploaded_file, values$sequence)
      #need((input$uploaded_file != "" || values$sequence != ""),"Please upload a sequence")# %then%
      #need(length(values$sequence)>3, "A minimum length of 15 bp is required")

    )
    start.time <- Sys.time()
    if(!is.null(input$uploaded_file)){
      iPASS(x = input$uploaded_file$datapath, InputType = "fasta")
    } else {
      iPASS(values$sequence, "text")
    }
  })
  
  alt_scores <- reactive({
      result <- resultCond()
      result$position <- factor(result$position, levels = result$position)
      resu <- result %>% 
        filter(pos>ranges$x[1] &pos<ranges$x[2])
      if (input$checkbox){
        resu
      }
     else {
      resu <- 
        resu %>% select(-alt_sequence, -alt_score)
    }
  })
  #resultCond <- reactive({ parseFasta(input$uploaded_file$datapath) })

  
  
  selected <-  reactiveValues(x = NULL, y = NULL)
  observe({
    if(input$submit > 0) {
      values$sequence <- isolate(input$seq_input)
      
    }
    
  })
  
  observeEvent(input$seqPlot_plot_click, {
     
      print(as.numeric(round(input$seqPlot_plot_click$x)))
      resu <- alt_scores()
      print(which.min(abs(resu$pos -input$seqPlot_plot_click$x)))
      proxy %>% selectRows(which.min(abs(resu$pos -input$seqPlot_plot_click$x)))
      
  })
  
  output$seqPlot <- renderPlot({
      resu <- alt_scores()
      resu$selected <- "a"
      resu$selected[input$info_rows_selected] <- "b"
      if("alt_score" %in% colnames(resu)){
        p <- resu %>% 
          gather("type","score", -position,-pos, -selected, -alt_sequence, -sequence, -aa_position) %>% 
          mutate(type = case_when(type == "score" ~ "wild-type",
                                  TRUE ~ "optimized"),
                   type = factor(type, levels = c("wild-type", "optimized"))) %>% 
          ggplot(aes(pos, score, color = selected, fill=type))+geom_col( size=1, width = 0.8, position = position_dodge())+
          theme(
            #axis.text.x=element_text(angle = -90, vjust = 0.5,
             #                        hjust=0),
            legend.title=element_blank(), 
            legend.position = "top"
          )+
          guides(color=FALSE)+
          xlab("Position")+
          ylab("iPASS score")+geom_hline(yintercept = 1, linetype="dashed", color="darkblue")+ggtitle("iPASS analysis")+
          scale_color_manual(values = c("lightgrey", "darkred"))+scale_fill_brewer(type="qual", palette = "Paired")
      } else{
      p <- resu %>% 
        ggplot(aes(pos, score, color = selected))+geom_col(fill="lightgrey", size=1, width = 0.8)+
        theme(
          #axis.text.x=element_text(angle = -90, vjust = 0.5,
          #                         hjust=0),
          legend.position = "none"
          )+xlab("Position")+
        ylab("iPASS score")+geom_hline(yintercept = 1, linetype="dashed", color="darkblue")+ggtitle("iPASS analysis")+
        scale_color_manual(values = c("lightgrey", "darkred"))
      }
      if (nrow(resu)<31){
        p + scale_x_continuous(breaks = resu$pos, labels= resu$position)
        
      } else {
        p
      }

   })
  output$info <- DT::renderDataTable(
    {
      resu <- alt_scores()
    if ("alt_score" %in% colnames(resu)){
    resu <- resu  %>% 
      dplyr::mutate(Codon = unlist(lapply(position, function(y) strsplit(as.character(y), "\n", fixed=T)[[1]][2])),
                    Position = pos,
                    AA = unlist(lapply(position, function(y) strsplit(as.character(y), "\n", fixed=T)[[1]][3])),
                    alt_score = round(alt_score, digits = 2)
             ) %>% 
      dplyr::select(score, Position, aa_position, Codon, AA, sequence, alt_score, alt_sequence) %>% 
      dplyr::rename(Score = score, 
                    `Position (bp)` = Position, 
                    `Position (AA)` = aa_position,
                    `Replaced AA` = AA,
                    `Replaced Codon` =Codon,
                    Sequence = sequence,
                    `Optimized Sequence` = alt_sequence,
                    `Optimized Score` = alt_score)
    
    datatable(resu, selection = 'single',rownames= FALSE, escape = F,
              options = list(
                columnDefs = list(list(className = 'dt-center', targets = "_all")))) %>% 
      formatRound(1, digits = 2) 
    } else {
      resu <- resu  %>% 
          dplyr::mutate(Codon = unlist(lapply(position, function(y) strsplit(as.character(y), "\n", fixed=T)[[1]][2])),
                        Position = pos,
                        AA = unlist(lapply(position, function(y) strsplit(as.character(y), "\n", fixed=T)[[1]][3]))
          ) %>% 
          dplyr::select(score, Position, aa_position, Codon, AA, sequence) %>% 
          dplyr::rename(Score = score, 
                        `Position (bp)` = Position,
                        `Position (AA)` = aa_position,
                        `Replaced AA` = AA,
                        `Replaced Codon` =Codon,
                        Sequence = sequence)}
      
      datatable(resu, selection = 'single',rownames= FALSE, escape = F,
                options = list(
                  columnDefs = list(list(className = 'dt-center', targets = "_all")))) %>% 
        formatRound("Score", digits = 2)
    
    }, 
    options = list(                # records/page options
                   lengthChange = FALSE,                       # show/hide records per page dropdown
                   bFilter=FALSE),server = TRUE
  )
}

# Run the application 
shinyApp(ui = ui, server = server)




