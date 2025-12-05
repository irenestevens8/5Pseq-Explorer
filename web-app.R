#install.packages(c("shiny","ggplot2", "DT", "plotly", "viridis", "reshape2", "ggrepel", "rsconnect", "shinydashboard", "dplyr", "purrr", "readxl", "zip", "httr", "readxl"))
library(shiny)
library(ggplot2)
library(dplyr)
library(DT)
library(plotly)
library(viridis)
library(reshape2)
library(ggrepel)
library(rsconnect)
library(shinydashboard)
library(purrr)
library(readxl)
library(zip)
library(httr)
library(readxl)

base_url <- "http://data.pelechanolab.com/software/5PSeq_explorer/"
yeast_path <- paste0(base_url, "data/yeast")
candida_path <- paste0(base_url, "data/candida")
bacteria_path <- paste0(base_url, "data/bacteria")
metadata_url <- paste0(base_url, "metadata_v2.xls")
rna_stats_url <- paste0(base_url, "rna_stats_v2.xls")

download_excel <- function(url) {
  tmp <- tempfile(fileext = ".xls")
  resp <- httr::GET(url, httr::write_disk(tmp, overwrite = TRUE))
  if (httr::status_code(resp) != 200) {
    stop(paste("Failed to download:", url))
  }
  return(tmp)
}

metadata <- readxl::read_excel(download_excel(metadata_url)) %>%
  dplyr::mutate(Sample = trimws(as.character(Sample))) %>%
  dplyr::distinct(Sample, .keep_all = TRUE)  

if (!"Sample" %in% names(metadata)) {
  stop("The metadata file must contain a 'Sample' column.")
}

rna_stats_raw <- suppressMessages(readxl::read_excel(download_excel(rna_stats_url)))

if (all(c("Sample", "Type", "Count") %in% names(rna_stats_raw))) {
  rna_stats <- rna_stats_raw
} else {
  rna_stats <- readxl::read_excel(download_excel(rna_stats_url),
                                  col_names = c("Sample", "Type", "Count"))
}

rna_stats <- rna_stats %>%
  dplyr::mutate(
    Sample = trimws(as.character(Sample)),
    Type = trimws(as.character(Type))
  )

missing_in_stats <- setdiff(metadata$Sample, rna_stats$Sample)
extra_in_stats <- setdiff(rna_stats$Sample, metadata$Sample)

if (length(missing_in_stats) > 0) {
  warning(paste("Samples in metadata missing from RNA composition file:",
                paste(missing_in_stats, collapse = ", ")))
}

if (length(extra_in_stats) > 0) {
  message(paste("Additional samples in RNA composition file not found in metadata:",
                paste(extra_in_stats, collapse = ", ")))
}

if (length(missing_in_stats) > 0) {
  rna_stats <- dplyr::bind_rows(
    rna_stats,
    data.frame(Sample = missing_in_stats, Type = NA, Count = NA)
  )
}

rna_stats <- rna_stats %>%
  dplyr::mutate(Sample = factor(Sample, levels = metadata$Sample)) %>%
  dplyr::arrange(Sample)

message("Metadata and RNA stats successfully loaded and aligned.")
message(paste("Total metadata samples:", nrow(metadata)))
message(paste("Total RNA stats entries:", nrow(rna_stats)))



if (!"Sample" %in% colnames(metadata)) {
  stop("The metadata.txt file must contain a 'Sample' column.")
}


get_folder_names <- function(path) {
  list.dirs(path, full.names = FALSE, recursive = FALSE)
}

# updated for URL handling 05/21/25
read_data_summary <- function(base_path, folder_name) {
  file_url <- paste(base_path, folder_name, "protein_coding", "data_summary.txt", sep = "/")
  tryCatch({
    summary_data <- read.delim(url(file_url), header = FALSE, stringsAsFactors = FALSE)
    colnames(summary_data) <- c("Metric", "Value")
    num_of_reads <- summary_data[summary_data$Metric == "NumOfReads", "Value"]
    num_of_map_positions <- summary_data[summary_data$Metric == "NumOfMapPositions", "Value"]
    data.frame(
      Folder = folder_name,
      NumOfReads = as.numeric(num_of_reads),
      NumOfMapPositions = as.numeric(num_of_map_positions)
    )
  }, error = function(e) {
    data.frame(Folder = folder_name, NumOfReads = NA, NumOfMapPositions = NA)
  })
}


# updated for URL handling 05/21/25
read_frame_stats <- function(base_path, folder_name) {
  file_url <- paste(base_path, folder_name, "protein_coding", "frame_stats.txt", sep = "/")
  tryCatch({
    read.delim(url(file_url), header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  }, error = function(e) NULL)
}


# updated for URL handling 05/21/25
read_meta_counts_start <- function(base_path, folder_name, num_of_map_positions) {
  file_url <- paste(base_path, folder_name, "protein_coding", "meta_counts_START.txt", sep = "/")
  tryCatch({
    meta_counts <- read.table(url(file_url), col.names = c("TSS", "Count"))
    meta_counts$CPM <- meta_counts$Count / num_of_map_positions * 1e6
    meta_counts
  }, error = function(e) NULL)
}

# updated for URL handling 05/21/25
read_meta_counts_term <- function(base_path, folder_name, num_of_map_positions) {
  file_url <- paste(base_path, folder_name, "protein_coding", "meta_counts_TERM.txt", sep = "/")
  tryCatch({
    meta_counts <- read.table(url(file_url), col.names = c("TERM", "Count"))
    meta_counts$CPM <- meta_counts$Count / num_of_map_positions * 1e6
    meta_counts
  }, error = function(e) NULL)
}

# updated for URL handling 05/21/25
read_amino_acid_pauses <- function(base_path, folder_name, num_of_map_positions) {
  file_url <- paste(base_path, folder_name, "protein_coding", "amino_acid_pauses.txt", sep = "/")
  tryCatch({
    pauses_data <- read.table(url(file_url), header = TRUE, row.names = 1, check.names = FALSE)
    pauses_data / num_of_map_positions * 1e6
  }, error = function(e) NULL)
}

# added on 10/30/25: Function to read codon_pauses.txt and calculate normalized values
# updated for URL handling 05/23/25
read_codon_pauses <- function(base_path, folder_name, num_of_map_positions) {
  file_url <- paste(base_path, folder_name, "protein_coding", "codon_pauses.txt", sep = "/")
  tryCatch({
    codon_pauses_data <- read.table(url(file_url), header = TRUE, row.names = 1, check.names = FALSE)
    codon_pauses_data <- codon_pauses_data / num_of_map_positions * 1e6  
    codon_pauses_data
  }, error = function(e) NULL)
}

# added 11/1/24:# Function to read and process gene-level frame counts and transcript assembly data for violin plot
# updated for URL handling 05/23/25

read_process_calculate <- function(base_path, folder_name) {
  transcript_url <- paste(base_path, folder_name, "protein_coding", "transcript_assembly.txt", sep = "/")
  frame_counts_url <- paste(base_path, folder_name, "protein_coding", "frame_counts_START.txt", sep = "/")
  
  tryCatch({
    transcript_data <- read.delim(url(transcript_url), header = TRUE, stringsAsFactors = FALSE)
    frame_counts <- read.delim(url(frame_counts_url), header = TRUE, row.names = NULL)
  }, error = function(e) {
    return(NULL)  
  })
  
  transcript_ids <- transcript_data$ID
  
  if (nrow(frame_counts) != length(transcript_ids)) {
    warning(paste("Mismatch between transcript IDs and frame counts in", folder_name))
    return(NULL)
  }
  
  combined_data <- data.frame(
    gene_ID = transcript_ids,
    F0 = frame_counts$F0,
    F1 = frame_counts$F1,
    F2 = frame_counts$F2,
    stringsAsFactors = FALSE
  )
  
  combined_data$Fsum <- combined_data$F0 + combined_data$F1 + combined_data$F2
  combined_data <- combined_data[combined_data$Fsum > 50, ]
  
  combined_data <- combined_data %>%
    mutate(
      F0_to_Fsum = F0/Fsum,
      F1_to_Fsum = F1/Fsum,
      F2_to_Fsum = F2/Fsum,
      F1_to_F0 = F1/(F0 + 0.1),
      F1_to_F2 = F1/(F2 + 0.1),
      F0_to_F1 = F0/(F1 + 0.1),
      F0_to_F2 = F0/(F2 + 0.1),
      F2_to_F1 = F2/(F1 + 0.1),
      F2_to_F0 = F2/(F0 + 0.1)
    )
  
  return(combined_data)
}



# update 10/02/24: frame prefs stacked plot
# update 10/29/24: combined metagene plots

ui <- dashboardPage(
  skin = "black",
  dashboardHeader(
    title = tags$div(
      tags$img(src = "https://pelechanolab.com/wp-content/uploads/2019/12/fivepseq_logoAsset-24@2x-100.jpg",
               height = "40px", style = "margin-right: 10px;"),
      "Explorer"
    ),
    titleWidth = 300
  ),
  dashboardSidebar(
    radioButtons("domain", label = h3("Choose domain"),
                 choices = list("Eukaryota" = "eukaryota", 
                                "Bacteria" = "bacteria"),
                 selected = "eukaryota"),
    conditionalPanel(
      condition = "input.domain == 'eukaryota'",
      selectInput("eukaryota_phylum", "Choose phylum",
                  choices = c("All", "Ascomycota"), selected = "Ascomycota")
    ),
    conditionalPanel(
      condition = "input.domain == 'bacteria'",
      selectInput("bacteria_phylum", "Choose phylum",
                  choices = NULL, selected = "All")
    ),
    selectInput("filter_species", "Species",
                choices = NULL, selected = "All"),
    selectInput("filter_strain", "Strain",
                choices = NULL, selected = "All"),
    radioButtons("replicate_filter", "Replicates Display",
                 choices = c("Individual" = "individual", "Merged" = "merged"),
                 selected = "individual"),
    
    h3("Additional Filters"),
    
    selectInput("filter_keyword", "Keyword",
                choices = NULL, selected = "All"),
    
    selectInput("filter_treatment", "Treatment",
                choices = NULL, selected = "All"),
    
    selectInput("filter_publication", "Publication",
                choices = NULL, selected = "All"),
    
    selectInput("filter_year", "Year",
                choices = NULL, selected = "All"),
    
    selectInput("filter_first_author", "First Author",
                choices = NULL, selected = "All"),
    
    selectInput("filter_lab_pi", "Lab PI",
                choices = NULL, selected = "All")
  ),
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .dataTables_wrapper {
          overflow-x: scroll;
        }
        .dataTables_length, .dataTables_filter, .dataTables_info, .dataTables_paginate {
          margin: 10px;
        }
        .dataTables_scrollHeadInner, .dataTables_scrollBody {
          width: 100% !important;
        }
        .dataTable td, .dataTable th {
          padding: 8px;
        }
        .dataTable {
          border-collapse: collapse;
          width: 100%;
        }
        .dataTable th, .dataTable td {
          border: 1px solid #ddd;
        }
        .dataTable tr:hover {
          background-color: #f5f5f5;
        }
      "))
    ),
    tags$div(style = "display:none;",
             checkboxGroupInput("checkGroup", label = NULL,
                                choices = unique(metadata$Sample),
                                selected = NULL)
    ),
    
    tabsetPanel(
      tabPanel("Welcome",
               fluidRow(
                 column(12,
                        div(style = "background-color: #ffffff;
                                    padding: 40px;
                                    border-radius: 12px;
                                    box-shadow: 0 4px 8px rgba(50,50,93,0.1);
                                    margin: 20px 0;",
                            tags$style(HTML("
                              .welcome-icon { color: #4a90e2; margin-right: 12px; }
                              .emphasis { color: #2c3e50; font-weight: 500; }
                              .contact-link { color: #4a90e2; text-decoration: underline; }
                            ")),
                            h2(icon("door-open", class = "welcome-icon"),
                               "Welcome to 5P-Seq Explorer",
                               style = "color: #2c3e50; margin-bottom: 25px; border-bottom: 2px solid #f0f0f0; padding-bottom: 15px;"),
                            div(style = "max-width: 800px; margin: 0 auto;",
                                p(style = "font-size: 16px; line-height: 1.7; color: #4a4a4a; margin-bottom: 20px;",
                                  "5P-seq is a method for profiling 5' phosphorylated mRNA undergoing decay. ",
                                  "Degradation coincides with translation, as the enzymatic cleavage of the 5' end of mRNAs occurs right behind the last trailing ribosome. ",
                                  "Therefore, information about ribosome dynamics can be obtained by profiling the degradome (5P-seq)."),
                                p(style = "font-size: 16px; line-height: 1.7; color: #4a4a4a; margin-bottom: 20px;",
                                  "Here you can explore ribosome dynamics in 723 datasets from bacteria, Saccharomyces and Candida albicans in reponse to various stresses and drug treatments. Navigate the Tabs (right) to explore ribosome positions at transcription start and termination sites, amino acid and codon level ribosome pauses, and translation frame usage across all genes."),
                                div(style = "background: #f3f8fe; padding: 20px; border-radius: 8px; margin: 25px 0;",
                                    p(style = "font-size: 16px; line-height: 1.7; color: #4a4a4a; margin: 0;",
                                      "Several filters can help you navigate through the available data. Two to highlight are ",
                                      strong(style = "color: #2c3e50;", "'Condition'"), " filter and ",
                                      strong(style = "color: #2c3e50;", "'Species'"), " filter.")
                                ),
                                p(style = "font-size: 16px; line-height: 1.7; color: #4a4a4a; margin-bottom: 20px;",
                                  "This application supports merging biological replicates — select the ",
                                  strong(style = "color: #2c3e50;", "'Merged'"), " button in the sidebar to average replicates."),
                                p(style = "font-size: 16px; line-height: 1.7; color: #4a4a4a; margin-bottom: 20px;",
                                  "To download the data used in generating the plots, navigate to the Download tab."),
                                p(style = "font-size: 16px; line-height: 1.7; color: #4a4a4a; margin-top: 25px;",
                                  "For further inquiries, contact ",
                                  tags$a(href = "irene.stevens@ki.se", "Irene Stevens",
                                         style = "color: #4a90e2; text-decoration: underline; font-weight: 500;"),
                                  ".")
                            )
                        )
                 )
               )
      ),
      
      tabPanel("Metadata",
               fluidRow(
                 box(DTOutput("metadata_table"), title = "Metadata", width = 12, height = "800px")
               )
      ),
      tabPanel("Mapping stats",
               box(plotlyOutput("barplot"), width = 12),
               box(plotlyOutput("rna_stacked_barplot"), width = 12)
      ),
      tabPanel("Frame Stats",
               box(uiOutput("frame_stats_plots"), width = 12)
      ),
      tabPanel("Metagene START Profile",
               fluidRow(
                 column(12, box(plotlyOutput("combined_tss_plot"), width = 12)),
                 column(12, box(uiOutput("metagene_tss_plots"), width = 12))
               )
      ),
      tabPanel("Metagene STOP Profile",
               fluidRow(
                 column(12, box(plotlyOutput("combined_term_plot"), width = 12)),
                 column(12, box(uiOutput("metagene_term_plots"), width = 12))
               )
      ),
      tabPanel("Amino Acid Protection Heatmap",
               box(checkboxInput("normalize_heatmap",
                                 label = "Normalize by row total (Amino Acid)",
                                 value = TRUE), width = 12),
               box(uiOutput("amino_acid_stall_heatmaps"), width = 12)
      ),
      tabPanel("Amino Acid Protection Lineplot",
               fluidRow(
                 column(4, box(selectInput("amino_acid", "Select Amino Acid", choices = NULL),
                               uiOutput("aa_sample_selector"),
                               width = 12)),
                 column(8, box(plotlyOutput("amino_acid_line_plot"), width = 12))
               )
      ),
      tabPanel("Amino Acid Protection Scatterplot",
               fluidRow(
                 column(4, box(numericInput("aa_position", "Select Position", value = -17, min = -30, max = 2),
                               selectInput("aa_sample1", "Select Sample 1", choices = NULL),
                               selectInput("aa_sample2", "Select Sample 2", choices = NULL),
                               width = 12)),
                 column(8, box(plotlyOutput("amino_acid_scatter_plot"), width = 12))
               )
      ),
      tabPanel("Codon Protection Heatmap",
               box(checkboxInput("normalize_codon_heatmap",
                                 label = "Normalize by row total (Codon)",
                                 value = TRUE), width = 12),
               box(uiOutput("codon_stall_heatmaps"), width = 12)
      ),
      tabPanel("Codon Protection Lineplot",
               fluidRow(
                 column(4, box(selectInput("codon", "Select Codon", choices = NULL),
                               uiOutput("sample_selector"),
                               width = 12)),
                 column(8, box(plotlyOutput("codon_stalls_plot"), width = 12))
               )
      ),
      tabPanel("Codon Protection Scatterplot",
               fluidRow(
                 column(4, box(numericInput("codon_position", "Select Position", value = -17, min = -30, max = 2),
                               selectInput("codon_sample1", "Select Sample 1", choices = NULL),
                               selectInput("codon_sample2", "Select Sample 2", choices = NULL),
                               width = 12)),
                 column(8, box(plotlyOutput("codon_scatter_plot"), width = 12))
               )
      ),
      tabPanel("Gene Frame Preferences",
               fluidRow(
                 column(4, box(selectInput("frame_selection", "Select Frame Proportion", choices = list(
                   "F0_to_Fsum" = "F0_to_Fsum",
                   "F1_to_Fsum" = "F1_to_Fsum",
                   "F2_to_Fsum" = "F2_to_Fsum",
                   "F1_to_F0" = "F1_to_F0",
                   "F1_to_F2" = "F1_to_F2",
                   "F0_to_F1" = "F0_to_F1",
                   "F0_to_F2" = "F0_to_F2",
                   "F2_to_F1" = "F2_to_F1",
                   "F2_to_F0" = "F2_to_F0"),
                   selected = "F1_to_Fsum"),
                   width = 12)),
                 column(8,
                        box(uiOutput("violin_plots"),  width = 12),
                        uiOutput("frame_plot_explanation")
                 )
               )
      ),
      tabPanel("Download Data",
               fluidRow(
                 box(
                   title = "Download Selected Data",
                   downloadButton("downloadData", "Download"),
                   p("After selecting samples on Metadata page, download here a ZIP archive of processed files used in the generation of all plots. Each sample's directory includes the following files:"),
                   tags$ul(
                     tags$li(HTML("<strong>data_summary.txt</strong>.<br>This file contains quick mapping statistics.")),
                     tags$li(HTML("<strong>frame_stats.txt</strong>.<br>This file containes summary frame usage statistics.")),
                     tags$li(HTML("<strong>meta_counts_START.txt</strong>.<br>This file contains raw read counts over the START of all genes.")),
                     tags$li(HTML("<strong>meta_counts_TERM.txt</strong>.<br>This file containes raw read counts over the TERMINATION site of all genes.")),
                     tags$li(HTML("<strong>amino_acid_pauses.txt</strong>.<br>This file containes raw reads over amino acids from position -30 to + 8.")),
                     tags$li(HTML("<strong>codon_pauses.txt</strong>.<br>This file contains raw reads over codons from position -30 to +8.")),
                     tags$li(HTML("<strong>frame_counts_START.txt</strong>.<br>This file containes the frame usage for all genes. The row number index corresponds to the transcript assembly row number.")),
                     tags$li(HTML("<strong>transcript_assembly.txt</strong>.<br>This file should be used concatenated with frame_counts_START.txt to obtain the gene frame usage."))
                   ),
                   width = 12
                 )
               )
      )
    )
  )
)

server <- function(input, output, session) {
  
  domainOrganisms <- reactive({
    if (input$domain == "eukaryota") {
      c("yeast", "candida")
    } else {
      "bacteria"
    }
  })
  
  get_base_path <- function(sample_name) {
    sample_organism <- metadata$Organism[metadata$Sample == sample_name]
    
    if (length(sample_organism) == 0 || is.na(sample_organism[1])) {
      warning(paste("No organism found for sample:", sample_name, "- defaulting to yeast"))
      return(yeast_path)
    }
    
    sample_organism <- trimws(tolower(sample_organism[1]))  
    
    switch(sample_organism,
           "yeast"    = yeast_path,
           "candida"  = candida_path,
           "bacteria" = bacteria_path,
           { warning(paste("Unknown organism:", sample_organism, "- defaulting to yeast"))
             yeast_path })
  }
  
  unique(metadata$Organism)
  
  
  observeEvent(input$domain, {
    organisms <- domainOrganisms()
    domain_metadata <- metadata %>% filter(Organism %in% organisms)
    
    if (input$domain == "bacteria") {
      phyla <- sort(unique(domain_metadata$Phylum))
      updateSelectInput(session, "bacteria_phylum", 
                        choices = c("All", phyla), selected = "All")
    }
  })
  
  observe({
    organisms <- domainOrganisms()
    domain_metadata <- metadata %>% filter(Organism %in% organisms)
    
    if (input$domain == "eukaryota" && !is.null(input$eukaryota_phylum) && input$eukaryota_phylum != "All") {
      domain_metadata <- domain_metadata %>% filter(Phylum == input$eukaryota_phylum)
    } else if (input$domain == "bacteria" && !is.null(input$bacteria_phylum) && input$bacteria_phylum != "All") {
      domain_metadata <- domain_metadata %>% filter(Phylum == input$bacteria_phylum)
    }
    
    species_choices <- c("All", sort(unique(domain_metadata$Species)))
    updateSelectInput(session, "filter_species", choices = species_choices, selected = "All")
  })
  
  observe({
    organisms <- domainOrganisms()
    domain_metadata <- metadata %>% filter(Organism %in% organisms)
    
    if (input$domain == "eukaryota" && !is.null(input$eukaryota_phylum) && input$eukaryota_phylum != "All") {
      domain_metadata <- domain_metadata %>% filter(Phylum == input$eukaryota_phylum)
    } else if (input$domain == "bacteria" && !is.null(input$bacteria_phylum) && input$bacteria_phylum != "All") {
      domain_metadata <- domain_metadata %>% filter(Phylum == input$bacteria_phylum)
    }
    
    if (!is.null(input$filter_species) && input$filter_species != "All") {
      domain_metadata <- domain_metadata %>% filter(Species == input$filter_species)
    }
    
    strain_choices <- c("All", sort(unique(domain_metadata$Strain)))
    updateSelectInput(session, "filter_strain", choices = strain_choices, selected = "All")
  })
  
  observeEvent(input$domain, {
    organisms <- domainOrganisms()
    domain_metadata <- metadata %>% filter(Organism %in% organisms)
    
    updateSelectInput(session, "filter_keyword", 
                      choices = c("All", sort(unique(domain_metadata$Keyword))),
                      selected = "All")
    updateSelectInput(session, "filter_treatment", 
                      choices = c("All", sort(unique(domain_metadata$treatment))),
                      selected = "All")
    updateSelectInput(session, "filter_publication",
                      choices = c("All", sort(unique(domain_metadata$Publication))),
                      selected = "All")
    updateSelectInput(session, "filter_year",
                      choices = c("All", sort(unique(domain_metadata$Year))),
                      selected = "All")
    updateSelectInput(session, "filter_first_author",
                      choices = c("All", sort(unique(domain_metadata$First_author))),
                      selected = "All")
    updateSelectInput(session, "filter_lab_pi",
                      choices = c("All", sort(unique(domain_metadata$Lab_PI))),
                      selected = "All")
  })
  
  currentFilteredMetadata <- reactive({
    filteredMetadata()
  })
  
  filteredMetadata <- reactive({
    organisms <- domainOrganisms()
    
    data <- metadata %>% 
      filter(Organism %in% organisms)
    
    if (input$domain == "eukaryota" && !is.null(input$eukaryota_phylum) && input$eukaryota_phylum != "All") {
      data <- data %>% filter(Phylum == input$eukaryota_phylum)
    } else if (input$domain == "bacteria" && !is.null(input$bacteria_phylum) && input$bacteria_phylum != "All") {
      data <- data %>% filter(Phylum == input$bacteria_phylum)
    }
    
    if (!is.null(input$filter_species) && input$filter_species != "All") {
      data <- data %>% filter(Species == input$filter_species)
    }
    
    if (!is.null(input$filter_strain) && input$filter_strain != "All") {
      data <- data %>% filter(Strain == input$filter_strain)
    }
    
    if (!is.null(input$filter_keyword) && input$filter_keyword != "All") {
      data <- data %>% filter(Keyword == input$filter_keyword)
    }
    if (!is.null(input$filter_treatment) && input$filter_treatment != "All") {
      data <- data %>% filter(treatment == input$filter_treatment)
    }
    if (!is.null(input$filter_publication) && input$filter_publication != "All") {
      data <- data %>% filter(Publication == input$filter_publication)
    }
    if (!is.null(input$filter_year) && input$filter_year != "All") {
      data <- data %>% filter(Year == input$filter_year)
    }
    if (!is.null(input$filter_first_author) && input$filter_first_author != "All") {
      data <- data %>% filter(First_author == input$filter_first_author)
    }
    if (!is.null(input$filter_lab_pi) && input$filter_lab_pi != "All") {
      data <- data %>% filter(Lab_PI == input$filter_lab_pi)
    }
    
    data
  })
  
  
  shinyInput <- function(FUN, len, id, ...) {
    inputs <- character(len)
    for(i in seq_len(len)) {
      inputs[i] <- as.character(FUN(paste0(id, i), label = NULL, ...))
    }
    inputs
  }
  
  output$metadata_table <- renderDT({
    currentMetadata <- currentFilteredMetadata()  
    metadata_display <- data.frame(
      Select = shinyInput(checkboxInput, nrow(currentMetadata), 'select_', value = FALSE),
      currentMetadata,
      stringsAsFactors = FALSE
    )
    datatable(
      metadata_display,
      escape = FALSE,  
      options = list(
        autoWidth = TRUE,
        scrollX = TRUE,
        scrollY = "600px", 
        scrollCollapse = TRUE,
        pageLength = 25,    
        lengthMenu = c(10, 25, 50, 100),  
        width = "100%",
        dom = "Blfrtip",
        buttons = c("copy", "csv", "excel", "pdf", "print"),
        filter = "top"
      ),
      callback = JS("
    table.on('change', 'input[type=\"checkbox\"]', function(){
      var selected = [];
      table.$('input[type=\"checkbox\"]').each(function(){
        if(this.checked){
          selected.push(this.id);
        }
      });
      Shiny.setInputValue('checkboxIDs', selected, {priority: 'event'});
    });
  "),
      class = "display nowrap",
      rownames = FALSE
    )
  }, server = FALSE)
  
  
  selected_groups <- reactive({
    if (input$replicate_filter == "merged") {
      unique(metadata$ReplicateGroup[metadata$Sample %in% input$checkGroup])
    } else {
      input$checkGroup
    }
  })
  
  
  observe({
    if (is.null(input$checkboxIDs) || length(input$checkboxIDs) == 0) {
      updateCheckboxGroupInput(session, "checkGroup", selected = character(0))
    } else {
      current_meta <- currentFilteredMetadata()
      rows <- as.numeric(gsub("select_", "", input$checkboxIDs))
      selected_samples <- current_meta$Sample[rows]
      updateCheckboxGroupInput(session, "checkGroup", selected = selected_samples)
    }
  })
  
  
  output$barplot <- renderPlotly({
    req(input$checkGroup) 
    
    data_list <- lapply(input$checkGroup, function(folder_name) {
      base_path <- get_base_path(folder_name)  
      read_data_summary(base_path, folder_name)
    })
    
    plot_data <- do.call(rbind, data_list)
    
    if (nrow(plot_data) > 0 && !all(is.na(plot_data$NumOfReads))) {
      p <- ggplot(plot_data, aes(x = Folder)) +
        geom_bar(aes(y = NumOfReads, text = paste("Total reads:", NumOfReads)), 
                 stat = "identity", position = "dodge", alpha = 0.7, fill = "#FFDDD6") +
        geom_bar(aes(y = NumOfMapPositions, text = paste("Mapped reads:", NumOfMapPositions)), 
                 stat = "identity", position = "dodge", alpha = 0.7, fill = "#FF876F") +
        theme_minimal() +
        labs(title = "Mapping summary", x = "Sample", y = "Read counts") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      ggplotly(p, tooltip = "text")
    } else {
      ggplotly(ggplot() +
                 labs(title = "No data available", x = "", y = "") +
                 theme_minimal())
    }
  })
  
  
  output$stacked_barplot <- renderPlotly({
    req(input$checkGroup)
    
    frame_data_list <- lapply(input$checkGroup, function(folder_name) {
      base_path <- get_base_path(folder_name)  
      stats_data <- read_frame_stats(base_path, folder_name)
      if (!is.null(stats_data)) {
        data.frame(
          Sample = folder_name,
          Frame = c("F0", "F1", "F2"),
          Percentage = as.numeric(stats_data["f_perc", ])
        )  
      }  
    })  
    
    frame_data <- do.call(rbind, frame_data_list)
    
    stacked_plot <- ggplot(frame_data, aes(x = Sample, y = Percentage, fill = Frame)) +
      geom_bar(stat = "identity", position = "stack") +
      theme_minimal() +
      labs(title = "Frame Percentage Across Samples", 
           x = "Sample", 
           y = "Frame Percentage") +
      scale_fill_manual(values = c("#414487FF", "#22A884FF","#FFDDD6")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggplotly(stacked_plot)
  })
  
  
  output$frame_stats_plots <- renderUI({
    req(input$checkGroup)
    
    plot_output_list <- lapply(input$checkGroup, function(folder_name) {
      base_path <- get_base_path(folder_name)  
      stats_data <- read_frame_stats(base_path, folder_name)
      if (!is.null(stats_data)) {
        plot_id_perc <- paste0("plot_perc_", folder_name)
        
        df <- data.frame(
          Frame = c("F0", "F1", "F2"), 
          f_perc = as.numeric(stats_data["f_perc", ])
        )  
        
        p <- ggplot(df, aes(x = Frame, y = f_perc)) +
          geom_bar(stat = "identity", fill = "#414487FF", alpha = 0.7) +
          theme_minimal() +
          labs(title = paste("Frame Percentage for", folder_name),
               x = "Frame", 
               y = "% Reads") +
          theme(axis.text = element_text(size = 10),
                plot.title = element_text(size = 11))
        
        output[[plot_id_perc]] <- renderPlotly(ggplotly(p))
        plotlyOutput(plot_id_perc)
      }  
    })  
    
    tagList(
      plotlyOutput("stacked_barplot"),
      hr(style = "border-top: 2px solid #bbb; margin: 20px 0;"),
      tags$h4("Individual Frame Statistics", 
              style = "font-size: 16px; margin-bottom: 15px;"),
      do.call(tagList, plot_output_list)
    )
  })
  
  
  ## Added on 10/29/24: combined TSS Plot
  output$combined_tss_plot <- renderPlotly({
    req(input$checkGroup)
    
    combined_data <- tryCatch({
      data_list <- lapply(input$checkGroup, function(folder_name) {
        base_path <- get_base_path(folder_name)
        summary_data <- read_data_summary(base_path, folder_name)
        meta_counts <- read_meta_counts_start(base_path, folder_name, summary_data$NumOfMapPositions)
        
        if (is.null(meta_counts) || nrow(meta_counts) == 0) return(NULL)
        meta_counts$Sample <- folder_name
        meta_counts
      })
      
      combined_data <- do.call(rbind, Filter(Negate(is.null), data_list))
      combined_data
    }, error = function(e) NULL)
    
    if (is.null(combined_data) || nrow(combined_data) == 0) {
      return(plotly_empty() %>% layout(title = "No data available for selected samples"))
    }
    
    combined_plot <- ggplot(combined_data, aes(x = TSS, y = CPM, color = Sample, text = paste("Sample:", Sample))) +
      geom_line(linewidth = 0.5) +
      scale_color_viridis(discrete = TRUE) +
      theme_bw() +
      labs(title = "Combined Metagene TSS Profile", x = "TSS", y = "Counts Per Million")
    
    ggplotly(combined_plot, tooltip = "text")
  })
  
  
  
  output$combined_tss_plot <- renderPlotly({
    req(input$checkGroup)
    
    combined_data <- tryCatch({
      if (input$replicate_filter == "merged") {
        selected_groups <- unique(metadata$ReplicateGroup[metadata$Sample %in% input$checkGroup])
        
        group_results <- lapply(selected_groups, function(group) {
          group_samples <- metadata$Sample[metadata$ReplicateGroup == group]
          
          group_data <- lapply(group_samples, function(folder_name) {
            tryCatch({
              base_path <- get_base_path(folder_name)
              summary_data <- read_data_summary(base_path, folder_name)
              if (is.null(summary_data)) return(NULL)
              
              meta_counts <- read_meta_counts_start(base_path, folder_name, summary_data$NumOfMapPositions)
              if (!is.null(meta_counts) && nrow(meta_counts) > 0) {
                meta_counts$Sample <- folder_name
                return(meta_counts)
              }
              NULL
            }, error = function(e) NULL)
          }) %>% bind_rows()
          
          if (!is.null(group_data) && nrow(group_data) > 0) {
            group_data %>%
              group_by(TSS) %>%
              summarise(CPM = mean(CPM, na.rm = TRUE), .groups = "drop") %>%
              mutate(Sample = group)
          } else NULL
        })
        
        bind_rows(group_results)
      } else {
        sample_results <- lapply(input$checkGroup, function(folder_name) {
          tryCatch({
            base_path <- get_base_path(folder_name)
            summary_data <- read_data_summary(base_path, folder_name)
            if (is.null(summary_data)) return(NULL)
            
            meta_counts <- read_meta_counts_start(base_path, folder_name, summary_data$NumOfMapPositions)
            if (!is.null(meta_counts) && nrow(meta_counts) > 0) {
              meta_counts$Sample <- folder_name
              return(meta_counts)
            }
            NULL
          }, error = function(e) NULL)
        })
        bind_rows(sample_results)
      }
    }, error = function(e) NULL)
    
    if (is.null(combined_data) || nrow(combined_data) == 0) {
      return(plotly_empty() %>% layout(title = "No data available for selected samples"))
    }
    
    combined_plot <- ggplot(combined_data, aes(x = TSS, y = CPM, color = Sample)) +
      geom_line(linewidth = 0.5) +
      scale_color_viridis(discrete = TRUE) +
      theme_bw() +
      labs(
        title = "Combined Metagene TSS Profile",
        x = "TSS",
        y = "Counts Per Million"
      )
    
    ggplotly(combined_plot)
  })
  
  
  output$metagene_tss_plots <- renderUI({
    req(input$checkGroup)
    
    if (input$replicate_filter == "merged") return(NULL)
    
    plot_output_list <- lapply(input$checkGroup, function(folder_name) {
      tryCatch({
        base_path <- get_base_path(folder_name)  
        summary_data <- read_data_summary(base_path, folder_name)
        meta_counts <- read_meta_counts_start(base_path, folder_name, summary_data$NumOfMapPositions)
        if (!is.null(meta_counts) && nrow(meta_counts) > 0) {
          plot_id <- paste0("plot_tss_", folder_name)
          
          output[[plot_id]] <- renderPlotly({
            p <- ggplot(meta_counts, aes(x = TSS, y = CPM)) +
              geom_line(color = "#414487FF") +
              theme_bw() +
              labs(title = paste("Metagene TSS profile for", folder_name),
                   x = "TSS", 
                   y = "Counts Per Million")
            ggplotly(p)
          })
          
          plotlyOutput(plot_id)
        } else {
          NULL
        }
      }, error = function(e) NULL)
    })
    
    plot_output_list <- Filter(Negate(is.null), plot_output_list)
    
    if (length(plot_output_list) == 0) {
      return(tags$p("No data available for selected samples"))
    }
    
    do.call(tagList, plot_output_list)
  })
  
  
  output$combined_term_plot <- renderPlotly({
    req(input$checkGroup)
    
    combined_data <- tryCatch({
      if (input$replicate_filter == "merged") {
        selected_groups <- unique(metadata$ReplicateGroup[metadata$Sample %in% input$checkGroup])
        
        group_results <- lapply(selected_groups, function(group) {
          group_samples <- metadata$Sample[metadata$ReplicateGroup == group]
          
          group_data <- lapply(group_samples, function(folder_name) {
            tryCatch({
              base_path <- get_base_path(folder_name)
              summary_data <- read_data_summary(base_path, folder_name)
              if (is.null(summary_data)) return(NULL)
              
              meta_counts <- read_meta_counts_term(base_path, folder_name, summary_data$NumOfMapPositions)
              if (!is.null(meta_counts) && nrow(meta_counts) > 0) {
                meta_counts$Sample <- folder_name
                return(meta_counts)
              }
              NULL
            }, error = function(e) NULL)
          }) %>% bind_rows()
          
          if (!is.null(group_data) && nrow(group_data) > 0) {
            group_data %>%
              group_by(TERM) %>%
              summarise(CPM = mean(CPM, na.rm = TRUE), .groups = "drop") %>%
              mutate(Sample = group)
          } else NULL
        })
        
        bind_rows(group_results)
      } else {
        sample_results <- lapply(input$checkGroup, function(folder_name) {
          tryCatch({
            base_path <- get_base_path(folder_name)
            summary_data <- read_data_summary(base_path, folder_name)
            if (is.null(summary_data)) return(NULL)
            
            meta_counts <- read_meta_counts_term(base_path, folder_name, summary_data$NumOfMapPositions)
            if (!is.null(meta_counts) && nrow(meta_counts) > 0) {
              meta_counts$Sample <- folder_name
              return(meta_counts)
            }
            NULL
          }, error = function(e) NULL)
        })
        bind_rows(sample_results)
      }
    }, error = function(e) NULL)
    
    if (is.null(combined_data) || nrow(combined_data) == 0) {
      return(plotly_empty() %>% layout(title = "No data available for selected samples"))
    }
    
    combined_plot <- ggplot(combined_data, aes(x = TERM, y = CPM, color = Sample)) +
      geom_line(linewidth = 0.5) +
      scale_color_viridis(discrete = TRUE) +
      theme_bw() +
      labs(
        title = "Combined Metagene TERM Profile",
        x = "TERM",
        y = "Counts Per Million"
      )
    
    ggplotly(combined_plot)
  })
  
  output$metagene_term_plots <- renderUI({
    req(input$checkGroup)
    
    if (input$replicate_filter == "merged") return(NULL)
    
    plot_output_list <- lapply(input$checkGroup, function(folder_name) {
      base_path <- get_base_path(folder_name)  
      summary_data <- read_data_summary(base_path, folder_name)
      meta_counts <- read_meta_counts_term(base_path, folder_name, summary_data$NumOfMapPositions)
      if (!is.null(meta_counts)) {
        plot_id <- paste0("plot_term_", folder_name)
        
        output[[plot_id]] <- renderPlotly({
          p <- ggplot(meta_counts, aes(x = TERM, y = CPM)) +
            geom_line(color = "#414487FF") +
            theme_bw() +
            labs(title = paste("Metagene TERM profile for", folder_name),
                 x = "TERM", 
                 y = "Counts Per Million")
          ggplotly(p)
        })
        
        plotlyOutput(plot_id)
      }
    })
    
    do.call(tagList, plot_output_list)
  })
  
  
  output$amino_acid_stall_heatmaps <- renderUI({
    req(input$checkGroup)
    
    if (input$replicate_filter == "merged") {
      selected_groups <- unique(metadata$ReplicateGroup[metadata$Sample %in% input$checkGroup])
      
      plot_output_list <- lapply(selected_groups, function(group) {
        group_samples <- metadata$Sample[metadata$ReplicateGroup == group]
        
        group_data <- lapply(group_samples, function(folder_name) {
          base_path <- get_base_path(folder_name)
          summary_data <- read_data_summary(base_path, folder_name)
          if (!is.null(summary_data)) {
            pauses_data <- read_amino_acid_pauses(base_path, folder_name, summary_data$NumOfMapPositions)
            if (!is.null(pauses_data)) {
              amino_acids <- rownames(pauses_data)
              pauses_numeric <- as.data.frame(sapply(pauses_data, as.numeric))
              
              if (input$normalize_heatmap) {
                pauses_numeric <- t(apply(pauses_numeric, 1, function(x) x / sum(x, na.rm = TRUE)))
              }
              
              pauses_numeric <- as.data.frame(pauses_numeric)
              rownames(pauses_numeric) <- amino_acids
              pauses_numeric$Amino_Acid <- rownames(pauses_numeric)
              
              melt(pauses_numeric, id.vars = "Amino_Acid", variable.name = "Position") %>%
                mutate(Position = as.numeric(as.character(Position)),
                       Sample = folder_name)
            }
          }
        }) %>% bind_rows()
        
        averaged_data <- group_data %>%
          filter(Position >= -20 & Position <= -3) %>%
          group_by(Amino_Acid, Position) %>%
          summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
        
        averaged_data$Position <- factor(averaged_data$Position, levels = -20:-3)
        
        n_aa <- length(unique(averaged_data$Amino_Acid))
        plot_height <- max(600, n_aa * 25)  
        
        p <- ggplot(averaged_data, aes(x = Position, y = Amino_Acid, fill = value)) +
          geom_tile() +
          scale_fill_viridis() +
          theme_minimal() +
          labs(title = paste("Amino Acid Stalls Heatmap for", group), 
               x = "Position", 
               y = "Amino Acid") +
          theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            axis.text.y = element_text(size = 9) 
          )
        
        plot_id <- paste0("plot_heatmap_", group)
        output[[plot_id]] <- renderPlotly({
          ggplotly(p, height = plot_height) %>%
            layout(
              margin = list(l = 100),  
              yaxis = list(tickfont = list(size = 10))
            )
        })
        plotlyOutput(plot_id, height = paste0(plot_height, "px"))
      })
      
    } else {
      plot_output_list <- lapply(input$checkGroup, function(folder_name) {
        base_path <- get_base_path(folder_name)
        summary_data <- read_data_summary(base_path, folder_name)
        if (!is.null(summary_data)) {
          pauses_data <- read_amino_acid_pauses(base_path, folder_name, summary_data$NumOfMapPositions)
          if (!is.null(pauses_data)) {
            amino_acids <- rownames(pauses_data)
            pauses_numeric <- as.data.frame(sapply(pauses_data, as.numeric))
            
            if (input$normalize_heatmap) {
              pauses_numeric <- t(apply(pauses_numeric, 1, function(x) x / sum(x, na.rm = TRUE)))
            }
            
            pauses_numeric <- as.data.frame(pauses_numeric)
            rownames(pauses_numeric) <- amino_acids
            pauses_numeric$Amino_Acid <- rownames(pauses_numeric)
            
            melted_data <- melt(pauses_numeric, id.vars = "Amino_Acid", variable.name = "Position")
            melted_data$Position <- as.numeric(as.character(melted_data$Position))
            
            melted_filtered <- melted_data[melted_data$Position >= -20 & melted_data$Position <= -3, ]
            
            melted_filtered$Position <- factor(melted_filtered$Position, levels = -20:-3)
            
            n_aa <- length(unique(melted_filtered$Amino_Acid))
            plot_height <- max(600, n_aa * 25)  
            
            p <- ggplot(melted_filtered, aes(x = Position, y = Amino_Acid, fill = value)) +
              geom_tile() +
              scale_fill_viridis() +
              theme_minimal() +
              labs(title = paste("Amino Acid Stalls Heatmap for", folder_name), 
                   x = "Position", 
                   y = "Amino Acid") +
              theme(
                axis.text.x = element_text(angle = 90, hjust = 1),
                axis.text.y = element_text(size = 9)  
              )
            
            plot_id <- paste0("plot_heatmap_", folder_name)
            output[[plot_id]] <- renderPlotly({
              ggplotly(p, height = plot_height) %>%
                layout(
                  margin = list(l = 100), 
                  yaxis = list(tickfont = list(size = 10))
                )
            })
            plotlyOutput(plot_id, height = paste0(plot_height, "px"))
          }
        }
      })
    }
    do.call(tagList, plot_output_list)
  })
  
  
  observeEvent(input$checkGroup, {
    req(input$checkGroup)
    available_samples <- input$checkGroup
    updateSelectInput(session, "sample1", choices = available_samples, selected = NULL)
    updateSelectInput(session, "sample2", choices = available_samples, selected = NULL)
    updateSelectInput(session, "aa_sample1", choices = available_samples, selected = NULL)
    updateSelectInput(session, "aa_sample2", choices = available_samples, selected = NULL)
    updateSelectInput(session, "codon_sample1", choices = available_samples, selected = NULL)
    updateSelectInput(session, "codon_sample2", choices = available_samples, selected = NULL)
  })
  
  output$amino_acid_pauses_plot <- renderPlotly({
    req(input$sample1, input$sample2, input$position)
    position <- as.character(input$position)
    
    sample1_data <- tryCatch({
      base_path1 <- get_base_path(input$sample1)  
      summary_data <- read_data_summary(base_path1, input$sample1)
      read_amino_acid_pauses(base_path1, input$sample1, summary_data$NumOfMapPositions)
    }, error = function(e) NULL)
    
    sample2_data <- tryCatch({
      base_path2 <- get_base_path(input$sample2)  
      summary_data <- read_data_summary(base_path2, input$sample2)
      read_amino_acid_pauses(base_path2, input$sample2, summary_data$NumOfMapPositions)
    }, error = function(e) NULL)
    
    if (!is.null(sample1_data) && 
        !is.null(sample2_data) &&
        position %in% colnames(sample1_data) && 
        position %in% colnames(sample2_data)) {
      
      combined_data <- data.frame(
        Amino_Acid = rownames(sample1_data),
        Sample1 = sample1_data[[position]],
        Sample2 = sample2_data[[position]]
      )
      
      p <- ggplot(combined_data, aes(x = Sample1, y = Sample2, label = Amino_Acid)) +
        geom_point(color = "#414487FF", size = 2.5) +
        geom_text_repel(aes(label = Amino_Acid), 
                        size = 3.5,
                        max.overlaps = 20,
                        box.padding = 0.5,
                        segment.color = "grey50") +
        theme_minimal(base_size = 12) +
        labs(
          title = paste("Amino Acid Pauses Comparison at Position", position),
          x = paste("Sample 1:", input$sample1),
          y = paste("Sample 2:", input$sample2)
        ) +
        geom_abline(linetype = "dashed", color = "grey50") +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))
      
      ggplotly(p, tooltip = c("x", "y", "label")) %>%
        layout(hoverlabel = list(bgcolor = "white"))
    } else {
      plotly_empty() %>%
        layout(title = "Could not load data for selected samples/position")
    }
  })
  
  output$violin_plots <- renderUI({
    req(input$checkGroup, input$frame_selection)
    
    plot_output_list <- tryCatch({
      if (input$replicate_filter == "merged") {
        selected_groups <- unique(metadata$ReplicateGroup[metadata$Sample %in% input$checkGroup])
        
        lapply(selected_groups, function(group) {
          samples_in_group <- metadata$Sample[metadata$ReplicateGroup == group]
          group_data_list <- lapply(samples_in_group, function(sample) {
            tryCatch({
              base_path <- get_base_path(sample)  
              data <- read_process_calculate(base_path, sample)  
              data
            }, error = function(e) NULL)
          }) %>% purrr::discard(is.null)
          
          if (length(group_data_list) == 0) return(NULL)
          
          group_data <- dplyr::bind_rows(group_data_list)
          n_transcripts <- nrow(group_data)
          
          if (n_transcripts == 0) return(NULL)
          
          p <- ggplot(group_data, aes(x = "", y = .data[[input$frame_selection]])) +
            geom_violin(fill = "#365c8d", trim = FALSE, width = 0.8) +
            geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
            annotate("text", x = 0.5, y = Inf, 
                     label = paste("N =", format(n_transcripts, big.mark = ",")), 
                     vjust = 1.5, size = 3, color = "gray30") +
            labs(
              title = paste0(group, "\n"),
              y = gsub("_", " ", input$frame_selection),
              x = NULL
            ) +
            theme_minimal() +
            theme(
              plot.title = element_text(size = 9, face = "bold", hjust = 0.5, 
                                        margin = margin(b = 5)),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              panel.grid.major.x = element_blank(),
              plot.margin = unit(c(2,5,2,5), "mm")
            )
          
          plot_id <- paste0("violin_", group)
          output[[plot_id]] <- renderPlotly({
            ggplotly(p, tooltip = "y") %>% 
              layout(
                showlegend = FALSE,
                margin = list(t = 40),
                annotations = list(
                  x = 0.5, y = 1.05,
                  text = paste0("<b>", group, "</b>"),
                  showarrow = FALSE,
                  xref = "paper",
                  yref = "paper",
                  font = list(size = 8)
                )
              )
          })
          box(plotlyOutput(plot_id), width = 8)
        })
        
      } else {
        lapply(input$checkGroup, function(sample) {
          tryCatch({
            base_path <- get_base_path(sample) 
            data <- read_process_calculate(base_path, sample)  
            if (is.null(data) || nrow(data) == 0) return(NULL)
            
            n_transcripts <- nrow(data)
            p <- ggplot(data, aes(x = "", y = .data[[input$frame_selection]])) +
              geom_violin(fill = "#FF876F", trim = FALSE, width = 0.8) +
              geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.5) +
              annotate("text", x = 0.5, y = Inf, 
                       label = paste("N =", format(n_transcripts, big.mark = ",")), 
                       vjust = 1.5, size = 3, color = "gray30") +
              labs(
                y = gsub("_", " ", input$frame_selection),
                x = NULL
              ) +
              theme_minimal() +
              theme(
                plot.title = element_text(size = 9, hjust = 0.5, margin = margin(b = 5)),
                axis.text.x = element_blank(),
                axis.title.x = element_blank(),
                panel.grid.major.x = element_blank(),
                plot.margin = unit(c(2,5,2,5), "mm")
              )
            
            plot_id <- paste0("violin_", sample)
            output[[plot_id]] <- renderPlotly({
              ggplotly(p, tooltip = "y") %>% 
                layout(
                  showlegend = FALSE,
                  margin = list(t = 40),
                  annotations = list(
                    x = 0.5, y = 1.05,
                    text = paste0("<b>", sample, "</b>"),
                    showarrow = FALSE,
                    xref = "paper",
                    yref = "paper",
                    font = list(size = 8)
                  )
                )
            })
            box(plotlyOutput(plot_id), width = 8)
          }, error = function(e) NULL)
        })
      }
    }, error = function(e) NULL)
    
    plot_output_list <- Filter(Negate(is.null), plot_output_list)
    
    if (length(plot_output_list) == 0) {
      return(tags$p("No data available for selected samples"))
    }
    
    tagList(
      fluidRow(
        lapply(plot_output_list, function(plot) {
          column(6, plot)
        })
      )
    )
  })
  
  output$frame_plot_explanation <- renderUI({
    HTML("<div style='margin-top: 20px; padding: 10px; background-color: #f8f9fa; border-radius: 5px;'>
       <strong>Explanation:</strong> Violin plots display proportions of frame counts across all genes (min ≥50 reads). 
       In eukaryotes, F1/Fsum is the normal codon protection phenotype.
       </div>")
  })
  
  output$rna_stacked_barplot <- renderPlotly({
    req(input$checkGroup)
    
    filtered_rna <- rna_stats %>%
      filter(Sample %in% input$checkGroup) %>%
      group_by(Sample) %>%
      mutate(
        Total = sum(Count),
        Percentage = Count / Total * 100
      ) %>%
      ungroup()
    
    p <- ggplot(
      filtered_rna, 
      aes(
        x = Sample, 
        y = Count, 
        fill = Type,
        text = paste(
          "Type:", Type, "<br>",
          "Count:", Count, "<br>",
          "Percentage:", round(Percentage, 1), "%"
        )
      )
    ) +
      geom_bar(stat = "identity", position = "fill") +
      scale_y_continuous(labels = scales::percent_format()) +
      scale_fill_viridis(discrete = TRUE) +
      theme_minimal() +
      labs(
        title = "RNA Composition", 
        x = "Sample", 
        y = "Read Counts"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggplotly(p, tooltip = "text")
  })
  
  observe({
    req(input$checkGroup)
    
    if (input$replicate_filter == "merged") {
      first_group <- unique(metadata$ReplicateGroup[metadata$Sample %in% input$checkGroup])[1]
      group_samples <- metadata$Sample[metadata$ReplicateGroup == first_group]
      first_sample <- group_samples[1]
    } else {
      first_sample <- input$checkGroup[1]
    }
    
    base_path <- get_base_path(first_sample)  
    summary_data <- read_data_summary(base_path, first_sample)
    codon_pauses_data <- read_codon_pauses(base_path, first_sample, summary_data$NumOfMapPositions)
    
    if (!is.null(codon_pauses_data)) {
      updateSelectInput(session, "codon", choices = rownames(codon_pauses_data))
    }
  })
  
  output$sample_selector <- renderUI({
    req(input$checkGroup)
    
    if (input$replicate_filter == "merged") {
      choices <- unique(metadata$ReplicateGroup[metadata$Sample %in% input$checkGroup])
      checkboxGroupInput("selected_samples", "Select Groups", 
                         choices = choices, selected = choices)
    } else {
      checkboxGroupInput("selected_samples", "Select Samples", 
                         choices = input$checkGroup, selected = input$checkGroup)
    }
  })
  
  output$codon_stalls_plot <- renderPlotly({
    req(input$codon, input$selected_samples)
    
    if (input$replicate_filter == "merged") {
      selected_samples <- metadata$Sample[metadata$ReplicateGroup %in% input$selected_samples]
    } else {
      selected_samples <- input$selected_samples
    }
    
    codon_data <- do.call(rbind, lapply(selected_samples, function(folder_name) {
      base_path <- get_base_path(folder_name)  
      summary_data <- read_data_summary(base_path, folder_name)
      if (!is.null(summary_data)) {
        pauses_data <- read_codon_pauses(base_path, folder_name, summary_data$NumOfMapPositions)
        if (!is.null(pauses_data)) {
          codon_row <- pauses_data[input$codon, , drop = FALSE]
          codon_row$OriginalSample <- folder_name
          return(codon_row)
        }
      }
    }))
    
    melted_data <- codon_data %>%
      melt(id.vars = "OriginalSample", variable.name = "Position", value.name = "CPM") %>%
      mutate(Position = as.numeric(as.character(Position)))
    
    if (input$replicate_filter == "merged") {
      melted_data <- melted_data %>%
        left_join(metadata %>% select(OriginalSample = Sample, ReplicateGroup), by = "OriginalSample") %>%
        group_by(ReplicateGroup, Position) %>%
        summarise(CPM = mean(CPM, na.rm = TRUE), .groups = "drop") %>%
        rename(Sample = ReplicateGroup)
    } else {
      melted_data <- melted_data %>% rename(Sample = OriginalSample)
    }
    
    p <- ggplot(melted_data, aes(x = Position, y = CPM, color = Sample)) +
      geom_line(size = 0.5) +
      theme_minimal() +
      labs(title = paste("Codon Stalls for:", input$codon),
           x = "Position Relative to Codon", 
           y = "Counts Per Million") +
      scale_color_viridis(discrete = TRUE) +
      theme(legend.position = "bottom")
    
    ggplotly(p) %>% layout(legend = list(orientation = "h", y = -0.2))
  })
  
  observe({
    req(input$checkGroup)
    
    if (input$replicate_filter == "merged") {
      first_group <- unique(metadata$ReplicateGroup[metadata$Sample %in% input$checkGroup])[1]
      group_samples <- metadata$Sample[metadata$ReplicateGroup == first_group]
      first_sample <- group_samples[1]
    } else {
      first_sample <- input$checkGroup[1]
    }
    
    base_path <- get_base_path(first_sample)  
    summary_data <- read_data_summary(base_path, first_sample)
    aa_pauses_data <- read_amino_acid_pauses(base_path, first_sample, summary_data$NumOfMapPositions)
    
    if (!is.null(aa_pauses_data)) {
      updateSelectInput(session, "amino_acid", choices = rownames(aa_pauses_data))
    }
  })
  
  output$aa_sample_selector <- renderUI({
    req(input$checkGroup)
    
    if (input$replicate_filter == "merged") {
      choices <- unique(metadata$ReplicateGroup[metadata$Sample %in% input$checkGroup])
      checkboxGroupInput("selected_aa_samples", "Select Groups", 
                         choices = choices, selected = choices)
    } else {
      checkboxGroupInput("selected_aa_samples", "Select Samples", 
                         choices = input$checkGroup, selected = input$checkGroup)
    }
  })
  
  output$amino_acid_line_plot <- renderPlotly({
    req(input$amino_acid, input$selected_aa_samples)
    
    if (input$replicate_filter == "merged") {
      selected_samples <- metadata$Sample[metadata$ReplicateGroup %in% input$selected_aa_samples]
    } else {
      selected_samples <- input$selected_aa_samples
    }
    
    aa_data <- do.call(rbind, lapply(selected_samples, function(folder_name) {
      base_path <- get_base_path(folder_name)  
      summary_data <- read_data_summary(base_path, folder_name)
      if (!is.null(summary_data)) {
        pauses_data <- read_amino_acid_pauses(base_path, folder_name, summary_data$NumOfMapPositions)
        if (!is.null(pauses_data)) {
          aa_row <- pauses_data[input$amino_acid, , drop = FALSE]
          aa_row$OriginalSample <- folder_name
          return(aa_row)
        }
      }
    }))
    
    melted_data <- aa_data %>%
      melt(id.vars = "OriginalSample", variable.name = "Position", value.name = "CPM") %>%
      mutate(Position = as.numeric(as.character(Position)))
    
    if (input$replicate_filter == "merged") {
      melted_data <- melted_data %>%
        left_join(metadata %>% select(OriginalSample = Sample, ReplicateGroup), by = "OriginalSample") %>%
        group_by(ReplicateGroup, Position) %>%
        summarise(CPM = mean(CPM, na.rm = TRUE), .groups = "drop") %>%
        rename(Sample = ReplicateGroup)
    } else {
      melted_data <- melted_data %>% rename(Sample = OriginalSample)
    }
    
    p <- ggplot(melted_data, aes(x = Position, y = CPM, color = Sample)) +
      geom_line(size = 0.5) +
      theme_minimal() +
      labs(title = paste("Amino Acid Stalls for:", input$amino_acid),
           x = "Position Relative to Amino Acid", 
           y = "Counts Per Million") +
      scale_color_viridis(discrete = TRUE) +
      theme(legend.position = "bottom")
    
    ggplotly(p) %>% layout(legend = list(orientation = "h", y = -0.2))
  })
  
  output$amino_acid_scatter_plot <- renderPlotly({
    req(input$aa_sample1, input$aa_sample2, input$aa_position)
    position <- as.character(input$aa_position)
    
    sample1_data <- tryCatch({
      base_path1 <- get_base_path(input$aa_sample1)  
      summary_data <- read_data_summary(base_path1, input$aa_sample1)
      read_amino_acid_pauses(base_path1, input$aa_sample1, summary_data$NumOfMapPositions)
    }, error = function(e) NULL)
    
    sample2_data <- tryCatch({
      base_path2 <- get_base_path(input$aa_sample2)  
      summary_data <- read_data_summary(base_path2, input$aa_sample2)
      read_amino_acid_pauses(base_path2, input$aa_sample2, summary_data$NumOfMapPositions)
    }, error = function(e) NULL)
    
    if (!is.null(sample1_data) && 
        !is.null(sample2_data) &&
        position %in% colnames(sample1_data) && 
        position %in% colnames(sample2_data)) {
      
      combined_data <- data.frame(
        Amino_Acid = rownames(sample1_data),
        Sample1 = sample1_data[[position]],
        Sample2 = sample2_data[[position]]
      )
      
      p <- ggplot(combined_data, aes(x = Sample1, y = Sample2, label = Amino_Acid)) +
        geom_point(color = "#FF876F", size = 2.5) +
        geom_text_repel(aes(label = Amino_Acid), 
                        size = 3.5,
                        max.overlaps = 20,
                        box.padding = 0.5,
                        segment.color = "grey50") +
        theme_minimal(base_size = 12) +
        labs(
          title = paste("Amino Acid Pauses Comparison at Position", position),
          x = paste("Sample 1:", input$aa_sample1, "(CPM)"),
          y = paste("Sample 2:", input$aa_sample2, "(CPM)")
        ) +
        geom_abline(linetype = "dashed", color = "grey50") +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))
      
      ggplotly(p, tooltip = c("x", "y", "label")) %>%
        layout(hoverlabel = list(bgcolor = "white"))
    } else {
      plotly_empty() %>%
        layout(title = "Could not load data for selected samples/position")
    }
  })
  
  output$codon_stall_heatmaps <- renderUI({
    req(input$checkGroup)
    
    if (input$replicate_filter == "merged") {
      selected_groups <- unique(metadata$ReplicateGroup[metadata$Sample %in% input$checkGroup])
      
      plot_output_list <- lapply(selected_groups, function(group) {
        group_samples <- metadata$Sample[metadata$ReplicateGroup == group]
        
        group_data <- lapply(group_samples, function(folder_name) {
          base_path <- get_base_path(folder_name)
          summary_data <- read_data_summary(base_path, folder_name)
          if (!is.null(summary_data)) {
            pauses_data <- read_codon_pauses(base_path, folder_name, summary_data$NumOfMapPositions)
            if (!is.null(pauses_data)) {
              codons <- rownames(pauses_data)
              pauses_numeric <- as.data.frame(sapply(pauses_data, as.numeric))
              
              if (input$normalize_codon_heatmap) {
                pauses_numeric <- t(apply(pauses_numeric, 1, function(x) x / sum(x, na.rm = TRUE)))
              }
              
              pauses_numeric <- as.data.frame(pauses_numeric)
              rownames(pauses_numeric) <- codons
              pauses_numeric$Codon <- rownames(pauses_numeric)
              
              melt(pauses_numeric, id.vars = "Codon", variable.name = "Position") %>%
                mutate(Position = as.numeric(as.character(Position)),
                       Sample = folder_name)
            }
          }
        }) %>% bind_rows()
        
        averaged_data <- group_data %>%
          filter(Position >= -20 & Position <= -3) %>%
          group_by(Codon, Position) %>%
          summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
        
        averaged_data$Position <- factor(averaged_data$Position, levels = -20:-3)
        
        n_codons <- length(unique(averaged_data$Codon))
        plot_height <- max(600, n_codons * 25)  
        
        p <- ggplot(averaged_data, aes(x = Position, y = Codon, fill = value)) +
          geom_tile() +
          scale_fill_viridis() +
          theme_minimal() +
          labs(title = paste("Codon Stalls Heatmap for", group), 
               x = "Position", 
               y = "Codon") +
          theme(
            axis.text.x = element_text(angle = 90, hjust = 1),
            axis.text.y = element_text(size = 9)  
          )
        
        plot_id <- paste0("plot_codon_heatmap_", group)
        output[[plot_id]] <- renderPlotly({
          ggplotly(p, height = plot_height) %>%
            layout(
              margin = list(l = 100),  
              yaxis = list(tickfont = list(size = 10))
            )
        })
        plotlyOutput(plot_id, height = paste0(plot_height, "px"))
      })
      
    } else {
      plot_output_list <- lapply(input$checkGroup, function(folder_name) {
        base_path <- get_base_path(folder_name)
        summary_data <- read_data_summary(base_path, folder_name)
        if (!is.null(summary_data)) {
          pauses_data <- read_codon_pauses(base_path, folder_name, summary_data$NumOfMapPositions)
          if (!is.null(pauses_data)) {
            codons <- rownames(pauses_data)
            pauses_numeric <- as.data.frame(sapply(pauses_data, as.numeric))
            
            if (input$normalize_codon_heatmap) {
              pauses_numeric <- t(apply(pauses_numeric, 1, function(x) x / sum(x, na.rm = TRUE)))
            }
            
            pauses_numeric <- as.data.frame(pauses_numeric)
            rownames(pauses_numeric) <- codons
            pauses_numeric$Codon <- rownames(pauses_numeric)
            
            melted_data <- melt(pauses_numeric, id.vars = "Codon", variable.name = "Position")
            melted_data$Position <- as.numeric(as.character(melted_data$Position))
            
            melted_filtered <- melted_data[melted_data$Position >= -20 & melted_data$Position <= -3, ]
            
            melted_filtered$Position <- factor(melted_filtered$Position, levels = -20:-3)
            
            n_codons <- length(unique(melted_filtered$Codon))
            plot_height <- max(600, n_codons * 25)  
            
            p <- ggplot(melted_filtered, aes(x = Position, y = Codon, fill = value)) +
              geom_tile() +
              scale_fill_viridis() +
              theme_minimal() +
              labs(title = paste("Codon Stalls Heatmap for", folder_name), 
                   x = "Position", 
                   y = "Codon") +
              theme(
                axis.text.x = element_text(angle = 90, hjust = 1),
                axis.text.y = element_text(size = 9)  
              )
            
            plot_id <- paste0("plot_codon_heatmap_", folder_name)
            output[[plot_id]] <- renderPlotly({
              ggplotly(p, height = plot_height) %>%
                layout(
                  margin = list(l = 100), 
                  yaxis = list(tickfont = list(size = 10))
                )
            })
            plotlyOutput(plot_id, height = paste0(plot_height, "px"))
          }
        }
      })
    }
    do.call(tagList, plot_output_list)
  })
  
  output$codon_scatter_plot <- renderPlotly({
    req(input$codon_sample1, input$codon_sample2, input$codon_position)
    position <- as.character(input$codon_position)
    
    sample1_data <- tryCatch({
      base_path1 <- get_base_path(input$codon_sample1)
      summary_data <- read_data_summary(base_path1, input$codon_sample1)
      read_codon_pauses(base_path1, input$codon_sample1, summary_data$NumOfMapPositions)
    }, error = function(e) NULL)
    
    sample2_data <- tryCatch({
      base_path2 <- get_base_path(input$codon_sample2)
      summary_data <- read_data_summary(base_path2, input$codon_sample2)
      read_codon_pauses(base_path2, input$codon_sample2, summary_data$NumOfMapPositions)
    }, error = function(e) NULL)
    
    if (!is.null(sample1_data) && 
        !is.null(sample2_data) &&
        position %in% colnames(sample1_data) && 
        position %in% colnames(sample2_data)) {
      
      combined_data <- data.frame(
        Codon = rownames(sample1_data),
        Sample1 = sample1_data[[position]],
        Sample2 = sample2_data[[position]]
      )
      
      p <- ggplot(combined_data, aes(x = Sample1, y = Sample2, label = Codon)) +
        geom_point(color = "#414487FF", size = 2.5) +
        geom_text_repel(aes(label = Codon), 
                        size = 3.5,
                        max.overlaps = 20,
                        box.padding = 0.5,
                        segment.color = "grey50") +
        theme_minimal(base_size = 12) +
        labs(
          title = paste("Codon Pauses Comparison at Position", position),
          x = paste("Sample 1:", input$codon_sample1, "(CPM)"), 
          y = paste("Sample 2:", input$codon_sample2, "(CPM)")
        ) +
        geom_abline(linetype = "dashed", color = "grey50") +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))
      
      ggplotly(p, tooltip = c("x", "y", "label")) %>%
        layout(hoverlabel = list(bgcolor = "white"))
    } else {
      plotly_empty() %>%
        layout(title = "Could not load data for selected samples/position")
    }
  })
  
  # added april 8, 2025: Download Data
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("5PSeq_Data_", Sys.Date(), ".zip", sep = "")
    },
    content = function(file) {
      req(input$checkGroup)
      samples <- input$checkGroup
      
      showModal(modalDialog("Preparing download...", footer = NULL))
      on.exit(removeModal())
      
      temp_dir <- tempfile()
      dir.create(temp_dir)
      root_dir <- file.path(temp_dir, "5P-seq_data")
      dir.create(root_dir)
      
      required_files <- c(
        "data_summary.txt",
        "frame_stats.txt",
        "meta_counts_START.txt",
        "meta_counts_TERM.txt",
        "amino_acid_pauses.txt",
        "codon_pauses.txt",
        "frame_counts_START.txt",
        "transcript_assembly.txt"
      )
      
      failed_downloads <- character()
      
      tryCatch({
        withProgress(
          message = 'Downloading files',
          value = 0,
          max = length(samples) * length(required_files),
          {
            for (sample in samples) {
              base_path <- get_base_path(sample)  
              sample_dir <- file.path(root_dir, sample, "protein_coding")
              dir.create(sample_dir, recursive = TRUE, showWarnings = FALSE)
              
              for (fname in required_files) {
                incProgress(1, detail = paste(sample, fname))
                
                file_url <- paste(base_path, sample, "protein_coding", fname, sep = "/")
                dest_path <- file.path(sample_dir, fname)
                
                tryCatch({
                  response <- GET(file_url, timeout(30))
                  if (status_code(response) == 200) {
                    writeBin(content(response, "raw"), dest_path)
                  } else {
                    failed_downloads <<- c(failed_downloads, file_url)
                  }
                }, error = function(e) {
                  failed_downloads <<- c(failed_downloads, file_url)
                })
              }
            }
          }
        )
        
        if (length(list.files(root_dir, recursive = TRUE)) == 0) {
          stop("No files could be downloaded")
        }
        
        old_wd <- setwd(temp_dir)
        on.exit(setwd(old_wd), add = TRUE)
        
        zip::zip(
          zipfile = file,
          files = "5P-seq_data",
          recurse = TRUE
        )
        
        if (length(failed_downloads) > 0) {
          showNotification(
            paste("Download completed with", length(failed_downloads), "missing files"),
            type = "warning", duration = 10
          )
        }
      },
      error = function(e) {
        showNotification(
          paste("Download failed:", e$message, 
                if (length(failed_downloads) > 0) paste("\nFailed files:", length(failed_downloads))),
          type = "error", duration = 10
        )
      },
      finally = {
        unlink(temp_dir, recursive = TRUE)
      })
    },
    contentType = "application/zip"
  )
  
}

shinyApp(ui = ui, server = server)
