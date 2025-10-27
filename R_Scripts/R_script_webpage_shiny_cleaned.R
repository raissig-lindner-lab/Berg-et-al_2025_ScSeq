### Code to create the web app for the single cell dataset
library(shiny)
library(tidyverse)
library(ggPlantmap)
library(Cairo)
library(bslib)
library(shinythemes)
library(ggpubr)
library(shinycssloaders)
library(arrow)

### read files for ggPlantmap
combined_data2 <- read.table("www/averaged_expression_per_celltype.csv")
dev.plantmap <- readRDS("www/leaf_dev_plantmap.rds")



### read files with data for labelled UMAP plots
label_info_all <- readRDS("www/All_labelled.rds")
label_info_epi <- readRDS("www/Epidermis_labelled.rds")
label_info_stom <- readRDS("www/Stomata_labelled.rds")



### Location file required for feature plots
location_file <- readRDS("www/Location_file_arrowfiles.rds")







map_expression <- function(data, accession_nr) {
  onegene <- data %>% filter(gene == accession_nr) 
  
  ### Merging with the ggplantmap data using the ggPlantmap.merge() function  
  final.table <- ggPlantmap.merge(dev.plantmap, onegene, id.x = "ROI.name",id.y="name")
  names(final.table)[names(final.table) == 'value'] <- 'Expression'
  
  for(row in 1:length(final.table$Expression)) {
    if(is.na(final.table$Expression[row])) {
      final.table$Expression[row] = 0
    }
  }
  
  ### Plot expression with ggplantmap
  ggPlantmap.heatmap(final.table, value.quant = Expression) + 
          scale_fill_gradient2(low="lightyellow", mid = "red", high="darkblue", midpoint = max(final.table$Expression)/2) +
          labs(title=NULL) +
          theme(plot.title = element_text(hjust = 0.1)) +
          annotate("segment", x = c(1160, 1160, 1160, 1550, 1550), y = c(-1170, -830, -450, -1020, -610), 
                   xend = c(1160, 1160, 1160, 1550, 1550), yend = c(-1020, -670, -300, -870, -450), 
                   arrow = arrow(type = "closed", length = unit(0.03, "npc")))+
          annotate("segment", x = c(450, 120, 700, 700, 700), y = c(-40, -120, -120, -120, -120), 
                   xend = c(450, 230, 130, 470, 750), yend = c(-160, -350, -330, -260, -290), colour = "#707070") +
          annotate("text", x = c(450, 120, 700, 380), y = c(0, -80, -80, -720), 
                   label = c("Epidermis", "Mesophyll", "Vasculature", "Shoot apex"), 
                   colour = "#707070", size = 3.5) +
          annotate("text", x = c(1040, 1000, 990, 990, 1390, 1400, 1400, 1820, 1850, 2250), y = c(-140, -560, -930, -1240, -310, -740, -1105, -350, -960, -650), 
                   label = c("Late\n GC", "Early\nGC", "Dividing\nGMC", "Early\nGMC", "Late\nHC", "Middle\nHC", "Early\nHC", "SC", "Silica cell", "Interstomatal cell"), 
                   colour = "#707070", size = 3.5, angle = 90, lineheight = 1)
    
}





featureplot_website <- function(dataset, gene=NULL) {
  if(gene == "x") {
    ggplot()+geom_point()+theme_classic()
  }
  
  else {
    ### check which files have the data for the selected gene and load the respective files
    file_name_all1 <- location_file[gene, "all_cells1"]
    file_name_all2 <- location_file[gene, "all_cells2"]
    file_name_epi <- location_file[gene, "epidermis"]
    file_name_stom <- location_file[gene, "stomatal_cell_files"]
    
  
    dt1 <- as.data.frame(read_feather(file_name_all1, col_select = c("cell", gene)))
    rownames(dt1) <- dt1$cell
    
    dt2 <- as.data.frame(read_feather(file_name_all2, col_select = c("cell", gene)))
    rownames(dt2) <- dt2$cell
    
    gene_info_all <- rbind(dt1, dt2)
    gene_info_all <- merge(gene_info_all, label_info_all[,1:2], by = c("row.names"))
    
    gene_info_epi <- as.data.frame(read_feather(file_name_epi, col_select = c("cell", gene)))
    rownames(gene_info_epi) <- gene_info_epi$cell
    gene_info_epi <- merge(gene_info_epi, label_info_epi[,1:2], by = c("row.names"))
    
    gene_info_stom <- as.data.frame(read_feather(file_name_stom, col_select = c("cell", gene)))
    rownames(gene_info_stom) <- gene_info_stom$cell
    gene_info_stom <- merge(gene_info_stom, label_info_stom[,1:2], by = c("row.names"))
    
    
    if(dataset == "All cells") {
      feat_data <- gene_info_all[,c("umap_1", "umap_2", gene)]
      
      feat_data <- feat_data[order(feat_data[,3]),]
      
      ggplot(feat_data)+
        geom_point(aes(x=umap_1, y=umap_2, colour = feat_data[,3]))+
        scale_colour_gradient2(low = "lightyellow", mid = "red", high = "darkblue", midpoint = max(feat_data[,3])/2)+
        theme_classic() +
        theme(legend.position = "bottom")+
        labs(colour = "Expression")
    }
    
    else {
      if(dataset == "Epidermis") {
        feat_data <- gene_info_epi[,c("umap_1", "umap_2", gene)]
        
        feat_data <- feat_data[order(feat_data[,3]),]
        
        ggplot(feat_data)+
          geom_point(aes(x=umap_1, y=umap_2, colour = feat_data[,3]))+
          scale_colour_gradient2(low = "lightyellow", mid = "red", high = "darkblue", midpoint = max(feat_data[,3])/2)+
          theme_classic()+
          theme(legend.position = "bottom")+
          labs(colour = "Expression")
      }
      
      else {
        if(dataset == "Stomatal cell files") {
          feat_data <- gene_info_stom[,c("umap_1", "umap_2", gene)]
          
          feat_data <- feat_data[order(feat_data[,3]),]
          
          ggplot(feat_data)+
            geom_point(aes(x=umap_1, y=umap_2, colour = feat_data[,3]))+
            scale_colour_gradient2(low = "lightyellow", mid = "red", high = "darkblue", midpoint = max(feat_data[,3])/2)+
            theme_classic() +
            theme(legend.position = "bottom")+
            labs(colour = "Expression")
        }
      }
    }
  }
}







umapplots_website <- function(dataset) {
  if(dataset == "All cells") {
    new_data <- label_info_all[,c("umap_1", "umap_2", "labels")]
    
    ggplot(new_data, aes(x=umap_1, y=umap_2, colour = labels))+
      geom_point(size = 1)+
      scale_colour_manual(values = c("#e76254", "#72bcd5", "#aadce0", "#376795", "#ffd06f", "#ffe6b7", "#ef8a47"))+
      theme_classic()+
      guides(colour = guide_legend(override.aes = list(shape = 15, size = 5), nrow = 2))+
      theme(legend.position = "bottom")+
      labs(colour = "Tissue")
  }
  
  else {
    if(dataset == "Epidermis") {
      new_data <- label_info_epi[,c("umap_1", "umap_2", "celltypes")]
      
      ggplot(new_data, aes(x=umap_1, y=umap_2, colour = celltypes))+
        geom_point(size = 1)+
        scale_colour_manual(values = c("#574571", "#7fa074", "#c1d1aa", "lightgrey", "#2c4b27", "#b695bc"))+
        theme_classic() +
        guides(colour = guide_legend(override.aes = list(shape = 15, size = 5), nrow = 2))+
        theme(legend.position = "bottom")+
        labs(colour = "Cell types")
    }
    
    else {
      if(dataset == "Stomatal cell files") {
        new_data <- label_info_stom[,c("umap_1", "umap_2", "stages")]
        
        ggplot(new_data, aes(x=umap_1, y=umap_2, colour = factor(stages, levels = c("Unknown", "Interstomatal cells", "Stage 0-1", 
                                                                 "SMCs", "SCs", "Early GMCs", 
                                                                 "Dividing GMCs", "Early GCs", "Late GCs"))))+
          geom_point(size = 1)+
          scale_colour_manual(values = c("#AECB72", "#78C7B8", "#8CC483",
                                         "#D8D97A", "#67AFC2", "#538EB9",
                                         "#356493", "#1D4573", "#0A2E57"))+
          theme_classic() +
          guides(colour = guide_legend(override.aes = list(shape = 15, size = 5), nrow = 3))+
          theme(legend.position = "bottom")+
          labs(colour = "Cell types")
      }
    }
  }
}




















#### User interface ------------------------------------------------------------------------
ui <- page_fixed(
  theme = shinytheme("sandstone"),
  
  titlePanel("Brachypodium leaf single-cell atlas"),
  
  layout_sidebar(
    sidebar = sidebar(
                      textInput("gene", label = "Gene Accession (e.g. BdiBd21-3.1G0240400):", value = "x"),
                      br(),
                      br(),
                      br(),
                      strong("Abbreviated terms:"),
                      "HC = Prickle hair cell",
                      br(),
                      "GMC = Guard mother cell",
                      br(),
                      "GC = Guard cell",
                      br(),
                      "SC = Subsidiary cell"
    ),
    plotOutput("gene_map", width = "800px", height = "400px") %>% withSpinner(color="#0dc5c1")),
  
  page_navbar(
    title = "Subset",
    id = "navbar",
    theme = shinytheme("sandstone"),
    nav_panel(title = "All cells", plotOutput("featureplot_all") %>% withSpinner(color="#0dc5c1")),
    nav_panel(title = "Epidermis", plotOutput("featureplot_epi") %>% withSpinner(color="#0dc5c1")),
    nav_panel(title = "Stomatal cell files", plotOutput("featureplot_stom") %>% withSpinner(color="#0dc5c1"))),
  
  br(), br(), br(), br(), br(), br(), br(),
  
  p("Click the button below to download the plots for your selected gene (may take a few seconds):"),
  downloadButton("download_plots", label = "Download PDF"),
  
  br(), br(), br(),
  
  p("Berg et al.", style = "font-size: 125%;"), 
  
  p(a("Stomatal Biology Lab, Institute of Plant Sciences, University of Bern, Switzerland", href = "https://raissiglab.org/", target = "_blank")),
  
  br(),br(),br(),br(),br()
)











#### Server data -----------------------------------------------------------------------------------
server <- function(input,output) {
  options(shiny.usecairo=T)
  
  ### Gene map
  data <- combined_data2
  
  gene_map_plot <- reactive({
    map_plot <- map_expression(data, input$gene)
    return(map_plot)})
  
  
  output$gene_map <- renderPlot({
    gene_map_plot()
    }
    , res = 96)
  
  all_cells_plots <- reactive({
    all_plot <- ggarrange(featureplot_website(dataset = "All cells", input$gene), umapplots_website(dataset = "All cells"))
    return(all_plot)})
  
  output$featureplot_all <- renderPlot({
    all_cells_plots()
  }, res = 96, width = 900, height = 500)
  
  
  epi_plots <- reactive({
    epi_plot <- ggarrange(featureplot_website(dataset = "Epidermis", input$gene), umapplots_website(dataset = "Epidermis"))
    return(epi_plot)})
  
  output$featureplot_epi <- renderPlot({
    epi_plots()
  }, res = 96, width = 900, height = 500)
  
  
  stom_plots <- reactive({
    stom_plot <- ggarrange(featureplot_website(dataset = "Stomatal cell files", input$gene), umapplots_website(dataset = "Stomatal cell files"))
    return(stom_plot)})
  
  output$featureplot_stom <- renderPlot({
    stom_plots()
  }, res = 96, width = 900, height = 500)
  
  
  plots_for_download <- function(identifier){
    print(ggarrange(gene_map_plot(), all_cells_plots(), epi_plots(), stom_plots(), nrow = 4))
  }
  
  output$download_plots <- downloadHandler(
    filename = function() {
      # Use the selected Gene and Dataset as the suggested file name
      paste0(input$gene, "_", Sys.Date(), ".pdf")
    },
   
    content = function(file) {
      # Write the dataset to the `file` that will be downloaded
      pdf(file, width = 10, height = 20)
      
      plots_for_download()
      
      dev.off()
    }
  )
}



shinyApp(ui, server)
