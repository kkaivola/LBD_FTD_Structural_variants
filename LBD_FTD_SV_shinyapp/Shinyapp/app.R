library(shiny)
library(shinyWidgets)
library(datamods)
  
#set working directory where files are/put files to your working directory

  
LBD <- read.delim("LBD_all_info_unfiltered_with_CADDSV_genehancer.txt", stringsAsFactors = FALSE)
FTD <- read.delim("FTD_all_info_unfiltered_with_CADDSV_genehancer.txt", stringsAsFactors = FALSE)
canonical_transcript_per_gene <- read.delim("canonical_transcript_per_gene.txt", stringsAsFactors = FALSE)
  
ui <- navbarPage("LBD/FTD candidate gene structural variants",
                   tabPanel("Data frame",
                            sidebarLayout(
                              sidebarPanel(
                                width = 3,
                                radioButtons(inputId = "dataset", label = "Phenotype", choices = c("LBD", "FTD")),
                                filter_data_ui("filtering")
                              ),
                              mainPanel(
                                DT::dataTableOutput(outputId = "table")
                              )
                            )
                   ),
                   tabPanel("Exonic SV visualisation",
                            sidebarLayout(
                              sidebarPanel(
                                width = 3,
                                radioButtons(inputId = "pheno", label = "Phenotype", choices = c("LBD", "FTD"), inline = TRUE),
                                pickerInput(inputId = "gene", label = "Gene", choices = unique(sort(LBD$GENE)), multiple = FALSE),
                                textOutput(outputId = "SV_info"),
                                textOutput(outputId = "canonical_transcript"),
                                htmlOutput(outputId = "color_codes")
                              ),
                              mainPanel(
                                imageOutput(outputId =  "img")
                              )
                            )
                   ),
                   tabPanel("Exonic, non-coding and promoter SV visualisation",
                            sidebarLayout(
                              sidebarPanel(
                                width = 3,
                                radioButtons(inputId = "pheno2", label = "Phenotype", choices = c("LBD", "FTD"), inline = TRUE),
                                pickerInput(inputId = "gene2", label = "Gene", choices = unique(sort(LBD$GENE)), multiple = FALSE),
                                textOutput(outputId = "SV_info2"),
                                textOutput(outputId = "canonical_transcript2"),
                                htmlOutput(outputId = "color_codes2")
                              ),
                              mainPanel(
                                imageOutput(outputId =  "img2")
                              )
                            )
                   )
  )
  
  
server <- function(input, output, session) {
    
    output$img <- renderImage({
      list(
        src = file.path("NIH_SV_images/",input$pheno,paste0(input$gene,"_",input$pheno,"_exonic.png")))
    }, deleteFile=FALSE)
    
    output$img2 <- renderImage({
      list(
        src = file.path("NIH_SV_images/",input$pheno2,paste0(input$gene2,"_",input$pheno2,"_candidate_gene_overlap.png")))
    }, deleteFile=FALSE)
    
    data <- reactive({
      get(input$dataset)
    })
    
    output$canonical_transcript <- renderText({ 
      paste0("Canonical transcript: ", canonical_transcript_per_gene[canonical_transcript_per_gene$Gene == input$gene,1])
    })
    
    output$SV_info <- renderText({
      "Number of carriers and within parentheses variant QUAL is presented on variants"
    })
    
    output$SV_info2 <- renderText({
      "Number of carriers and within parentheses variant QUAL is presented on variants"
    })
    
    output$canonical_transcript2 <- renderText({ 
      paste0("Canonical transcript: ", canonical_transcript_per_gene[canonical_transcript_per_gene$Gene == input$gene2,1])
    })
    
    output$color_codes <- renderText({ 
      str1 <- paste0("<font color=\"#5282AF\"><b>","Duplications in blue","</b></font>")
      str2 <- paste0("<font color=\"#D6402D\"><b>","Deletions in red","</b></font>")
      str3 <- paste0("<font color=\"#377246\"><b>","BNDs in dark green","</b></font>")
      str4 <- paste0("<font color=\"#D579E1\"><b>","Insertions in purple","</b></font>")
      str5 <- paste0("<font color=\"#F69A57\"><b>","Inversions in orange","</b></font>")
      str6 <- paste0("<font color=\"#86E496\"><b>","CPX in light green","</b></font>")
      str7 <- paste0("<font color=\"#F6AAFF\"><b>","Mobile element insertions in pink","</b></font>")
      paste0(str1, "<br/>", str2, "<br/>", str3, "<br/>", str4, "<br/>", str5, "<br/>", str6, "<br/>", str7)
    })
    
    output$color_codes2 <- renderText({ 
      str1 <- paste0("<font color=\"#5282AF\"><b>","Duplications in blue","</b></font>")
      str2 <- paste0("<font color=\"#D6402D\"><b>","Deletions in red","</b></font>")
      str3 <- paste0("<font color=\"#377246\"><b>","BNDs in dark green","</b></font>")
      str4 <- paste0("<font color=\"#D579E1\"><b>","Insertions in purple","</b></font>")
      str5 <- paste0("<font color=\"#FG9A57\"><b>","Inversions in orange","</b></font>")
      str6 <- paste0("<font color=\"#86E496\"><b>","CPX in light green","</b></font>")
      str7 <- paste0("<font color=\"#F6AAFF\"><b>","Mobile element insertions in pink","</b></font>")
      paste0(str1, "<br/>", str2, "<br/>", str3, "<br/>", str4, "<br/>", str5, "<br/>", str6, "<br/>", str7)
    })
    
    vars <- reactive({
      if (identical(input$dataset, "LBD")) {
        colnames(LBD)[c(1,6:17,19:22)]
      } else {
        colnames(FTD)[c(1,6:17,19:22)]
      }
    })
    
    res_filter <- filter_data_server(
      id = "filtering",
      data = data,
      name = reactive(input$dataset),
      widget_num = "range",
      widget_date = "slider",
      widget_char = "picker",
      label_na = "Missing",
      drop_ids = FALSE,
      vars = vars
    )
    
    
    output$table <- DT::renderDT({
      res_filter$filtered()
    }, options = list(pageLength = 50))
  }
  
shinyApp(ui, server)

