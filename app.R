#### P. knowlesi Orthologue Expression Viewer ####

# Author: Katlijn De Meulenaere

# Script to generate a Shiny App for investigation of P. knowlesi (A1-H.1; De Meulenaere et al., 2026) gene expression patterns with P. vivax (smru1 clinical isolate; Zhu et al., 2016) and P. falciparum (Subudhi et al., 2020) gene expression patterns, over the IDC (asexual blood stages).
# In the tab 'Orthologue Expression Viewer Pk-Pv', a Pv/Pk gene is entered, the orthologue(s) are searched, and an expression plot for the input gene + orthologue(s) is given together with similarity statistics of the expression patterns.
# In the tab 'Orthologue Expression Viewer Pk-Pf', a Pf/Pk gene is entered, the orthologue(s) are searched, and an expression plot for the input gene + orthologue(s) is given together with similarity statistics of the expression patterns.
# In the tab 'Any Gene Expression Viewer', up to 4 Pk/Pv/Pf genes are entered, and one expression plot is given for all input genes together.

# To launch the App, do:
# library(shiny)
# runApp("path/Pk_comparator") #app.R is stored in /Pk_comparator directory.


#### Load packages ####

library(shiny) #v1.13.0
library(ggplot2) #v4.0.2
library(gridExtra) #v2.3
library(cowplot) #v1.2.0
# R version 4.5.2


#### Import dataframes ####

Pk_expr <- read.delim("./datasets/Pk_expr.txt", check.names = FALSE)
Pv_expr <- read.delim("./datasets/Pv_expr.txt", check.names = FALSE)
Pf_expr <- read.delim("./datasets/Pf_expr.txt", check.names = FALSE)
#Expression of all gene ID's that passed filtering (expression: TPM, log2 with a pseudocount of 1, Z-scored)). Genes without orthologue are present.

sim_PkPv <- read.delim("./datasets/sim_PkPv.txt")
sim_PkPf <- read.delim("./datasets/sim_PkPf.txt")
#orthologous gene ID's (that passed filtering) with their similarity stats: dtw, time shift, outcome (sim/dissim)

ortho_PkPv <- read.delim("./datasets/ortho_PkPv.txt")
ortho_PkPf <- read.delim("./datasets/ortho_PkPf.txt")
#orthologous gene ID's (no filtering applied)

Pk_unfiltered <- read.delim("./datasets/Pk_unfiltered.txt")
Pv_unfiltered <- read.delim("./datasets/Pv_unfiltered.txt")
Pf_unfiltered <- read.delim("./datasets/Pf_unfiltered.txt")
#dataframe with all Pv/Pk genes, and whether they pass the filtering rules.


#later, for web version:

#Pk_expr <- read.delim("./datasets/Pk_expr.txt", check.names = FALSE)
#Pv_expr <- read.delim("./datasets/Pv_expr.txt", check.names = FALSE)
#Pf_expr <- read.delim("./datasets/Pf_expr.txt", check.names = FALSE)
#Expression of all gene ID's that passed filtering (expression: TPM, log2 with a pseudocount of 1, Z-scored)). Genes without orthologue are present.

#sim_PkPv <- read.delim("./datasets/sim_PkPv.txt")
#sim_PkPf <- read.delim("./datasets/sim_PkPf.txt")
#orthologous gene ID's (that passed filtering) with their similarity stats: dtw, time shift, outcome (sim/dissim)

#ortho_PkPv <- read.delim("./datasets/ortho_PkPv.txt")
#ortho_PkPf <- read.delim("./datasets/ortho_PkPf.txt")
#orthologous gene ID's (no filtering applied)

#Pk_unfiltered <- read.delim("./datasets/Pk_unfiltered.txt")
#Pv_unfiltered <- read.delim("./datasets/Pv_unfiltered.txt")
#Pf_unfiltered <- read.delim("./datasets/Pf_unfiltered.txt")
#dataframe with all Pv/Pk/Pf genes, and whether they pass the filtering rules.


#### Helper functions ####

# A. Any gene tab: Extract data for input genes. Input genes are not necessarily orthologues.
get_gene_data <- function(gene_ids) {
  #gene_ids is a vector containing 1 or more gene IDs (Pk or Pv), e.g. gene_ids<-c("PKNH_0931500","PVP01_0109300")
  
  results <- list()
  #fill list elements with empty values:
  for (gene_id in gene_ids) {
    warnings <- character(0)
    df_expr <- NULL
    sampled_points <- NULL
    max_val <- NA
    
    if (grepl("^PKNH_", gene_id)) {
      # === If Pk gene ===
      
      if (!(gene_id %in% rownames(Pk_unfiltered))) {
        warnings <- c(warnings, paste0("Unrecognised gene ID: ", gene_id, " (typo/new gene)."))
        #same gene names as PkPf_unfiltered
        
      } else {
        
        status <- Pk_unfiltered[gene_id, ]
        
        if (status$FCPass == "low") {
          warnings <- c(warnings, paste0("Gene ", gene_id, " has too low transcriptional variation over IDC (log2 fold change between highest and lowest normalised expression level <0.5)."))
        }
        if (status$sumPass == "low") {
          warnings <- c(warnings, paste0("Gene ", gene_id, " has too low expression (sum of normalised expression levels of the 5 sampled time points <0.5)."))
        }
        
        if (length(warnings) == 0) {
          
          #get results list input:
          
          expr <- Pk_expr[gene_id, , drop = FALSE]
          expr_vec <- as.numeric(expr)
          reltime <- as.numeric(colnames(expr))
          df_expr <- data.frame(reltime = reltime, values = expr_vec, gene = gene_id)
          
          max_val <- reltime[which.max(expr_vec)]
          
          # sampled IDC points (5)
          sampled_points <- df_expr[df_expr$reltime %in% c(0.185, 0.520, 0.740, 0.890, 1.000), ]
        }
        
      }
    } else if (grepl("^PVP01_", gene_id)) {
      # === If Pv gene ===
      
      if (!(gene_id %in% rownames(Pv_unfiltered))) {
        warnings <- c(warnings, paste0("Unrecognised gene ID: ", gene_id, " (typo/new gene)."))
        
      } else {
        
        status <- Pv_unfiltered[gene_id, ]
        
        if (status$FCPass == "low") {
          warnings <- c(warnings, paste0("Gene ", gene_id, " has too low transcriptional variation over IDC (log2 fold change between highest and lowest normalised expression level <0.5)."))
        }
        if (status$sumPass == "low") {
          warnings <- c(warnings, paste0("Gene ", gene_id, " has too low expression (sum of normalised expression levels of the 7 sampled time points <0.5)."))
        }
        
        if (length(warnings) == 0) {
          
          #get results list input:
          
          expr <- Pv_expr[gene_id, , drop = FALSE]
          expr_vec <- as.numeric(expr)
          reltime <- as.numeric(colnames(expr))
          df_expr <- data.frame(reltime = reltime, values = expr_vec, gene = gene_id)
          
          max_val <- reltime[which.max(expr_vec)]
          
          # sampled IDC points (7)
          sampled_points <- df_expr[df_expr$reltime %in% c(0.125, 0.375, 0.500, 0.625, 0.750, 0.875, 1.000), ]
        }
      }
      
    } else if (grepl("^PF3D7_", gene_id)) {
      # === If Pf gene ===
      
      if (!(gene_id %in% rownames(Pf_unfiltered))) {
        warnings <- c(warnings, paste0("Unrecognised gene ID: ", gene_id, " (typo/new gene)."))
        
      } else {
        
        status <- Pf_unfiltered[gene_id, ]
        
        if (status$FCPass == "low") {
          warnings <- c(warnings, paste0("Gene ", gene_id, " has too low transcriptional variation over IDC (log2 fold change between highest and lowest normalised expression level <0.5)."))
        }
        if (status$sumPass == "low") {
          warnings <- c(warnings, paste0("Gene ", gene_id, " has too low expression (sum of normalised expression levels across sampled time points <0.5)."))
        }
        
        if (length(warnings) == 0) {
          
          #get results list input:
          
          expr <- Pf_expr[gene_id, , drop = FALSE]
          expr_vec <- as.numeric(expr)
          reltime <- as.numeric(colnames(expr))
          df_expr <- data.frame(reltime = reltime, values = expr_vec, gene = gene_id)
          
          max_val <- reltime[which.max(expr_vec)]
          
          # sampled IDC points(25)
          sampled_points <- df_expr[df_expr$reltime %in% c(0.000,0.040,0.085,0.125,0.165,0.210,0.250,0.290,0.335,0.375,0.415,0.460,0.500,0.540,0.585,
                                                           0.625,0.665,0.710,0.750,0.790,0.835,0.875,0.915,0.960,1.000), ]
          
        }
      }
      
    } else {
      warnings <- c(warnings, paste0("Unrecognised gene ID: ", gene_id, " (must start with PVP01_..., PKNH_..., or PF3D7_...)."))
    }
    
    results[[gene_id]] <- list(
      #Insert the below objects as the elements of the list named gene_id.
      expression = df_expr,
      sampled_points = sampled_points,
      max = max_val,
      gene = gene_id,
      warning = if (length(warnings) > 0) paste(warnings, collapse = "\n") else NULL
    )
  }
  return(results)
  
  #structure 'results': nested list. Each element is an input gene. Each element itself is a sublist.
  #-gene 1 list:
  #    -expression (df)
  #    -sampled points (df)
  #    -max (num)
  #    -gene (chr)
  #    -warning (chr)
  #-gene 2 list: same structure, only if >1 gene was entered
  #etc
}


# B. Any gene tab: Generate plot
make_plot <- function(results, add_citation = FALSE) {
  # nested results list made in function A.
  # add_citation default is FALSE.
  
  # only make plot for genes without warning
  valid <- Filter(function(x) is.null(x$warning), results) #go through list results, and only keep genes with no warning; throw out multi_warn warning if >1 orthologues (see C.).
  if (length(valid) == 0) return(NULL)                     #if every gene has a warning, no plot will be rendered. Doesn't continue.
  
  # bind expression/sampled_points data
  df_expr <- do.call(rbind, lapply(valid, `[[`, "expression"))        #extract 'expression' dataframe from 'results' list, stack it row-wise
  df_sampled  <- do.call(rbind, lapply(valid, `[[`, "sampled_points"))  #extract 'sampled_points' dataframe from 'results' list, stack it row-wise
  df_max      <- data.frame(reltime = sapply(valid, `[[`, "max"),
                            gene = sapply(valid, `[[`, "gene"))         #extract 'max' and 'gene' from 'results' list per gene_id, into a 2-column dataframe (cols reltime and gene)
  
  # Offset overlapping maxima, such that vertical dotted lines in plot are not on top of each other
  if (nrow(df_max) > 1) {
    dup_groups <- split(df_max, df_max$reltime) #splits df_max into a list of different data.frames, one per unique value in $reltime.
    df_max$offset <- 0
    for (grp in dup_groups) {
      if (nrow(grp) > 1) { #meaning there are genes with the same maximum
        # spread them symmetrically around the true reltime
        n <- nrow(grp)
        shifts <- seq(-0.002, 0.002, length.out = n)  # if 2 genes with same maxima, they will be plotted at max-0.002 and max+0.002, if 3 genes, then at max-0.002, max and max+0.002, etc
        df_max$offset[df_max$gene %in% grp$gene] <- shifts
      }
    }
    df_max$reltime_offset <- df_max$reltime + df_max$offset
  } else {
    df_max$reltime_offset <- df_max$reltime
  }
  
  #plot
  p<-ggplot(df_expr, aes(x = reltime, y = values, group = gene)) +
    geom_line(aes(color = gene), linewidth = 0.8) +
    geom_point(data = df_sampled, aes(x = reltime, y = values, color = gene), size = 2) +
    geom_vline(data = df_max, aes(xintercept = reltime_offset, color = gene, linetype = "peak expression time point"),
               linewidth = 1, show.legend = c(color = FALSE)) +
    scale_x_continuous(limits = c(0, 1.01)) + #>1 for the potential offsets
    scale_color_manual(
      values = c("#D81B60", "#1E88E5", "#FFC107", "#004D40")) +
    scale_linetype_manual(
      name = NULL,
      values = c("peak expression time point" = "dotted")) +
    guides(color = guide_legend(
      order = 1,
      override.aes = list(
        linetype = "solid",
        linewidth = 0.8)),
      linetype = guide_legend(
        order = 2,
        override.aes = list(
          color = "black",
          linewidth = 1))) +
    xlab("Relative time") +
    ylab("Gene expression\nnormalised and Z-scored") +
    theme_bw() +
    theme(
      axis.text = element_text(size = 12, color = "black"),
      axis.title = element_text(size = 14, color = "black"),
      legend.title = element_blank(),
      legend.text = element_text(size = 14, color = "black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.key.height = unit(0.8, "cm")
    )
  
  # Optionally add citation (in case figure is downloaded)
  if (add_citation) {
    citation <- cowplot::ggdraw() +
      theme(plot.background = element_rect(fill = "white", color = NA)) +
      cowplot::draw_label(
        "De Meulenaere et al., 2026, doi: 10.64898/2026.02.09.704036.",
        x = 0, y = 0.5, hjust = 0, vjust = 0.5,
        size = 10, color = "gray40"
      )
    p <- cowplot::plot_grid(p, citation, ncol = 1, rel_heights = c(1, 0.1))
  }
  
  return(p)
}

# C. orthologue tab Pk-Pv: Extract data for a single Pk/Pv input gene and its Pv/Pk orthologue(s).
get_ortho_data_Pv <- function(gene_id) {
  # gene_id_PkPv : single Pv or Pk gene ID (string)
  
  # check existence of the input gene:
  if (grepl("^PKNH_", gene_id)) {
    
    if (!(gene_id %in% rownames(Pk_unfiltered))) {
      return(list(list(warning = "Unrecognised input gene ID (typo/new gene)."))) #create a nested list that just contains the warning (as this is the output format of the function)
      #same gene names as PkPf_unfiltered
    } #doesn't continue
    
    # get orthologue(s) in Pv
    ortho_hits <- ortho_PkPv$Pv_GeneID[ortho_PkPv$Pk_GeneID == gene_id]
    
    species <- "Pk"
    
  } else if (grepl("^PVP01_", gene_id)) {
    
    if (!(gene_id %in% rownames(Pv_unfiltered))) {
      return(list(list(warning = "Unrecognised input gene ID (typo/new gene)."))) #create a nested list that just contains the warning (as this is the output format of the function)
    } #doesn't continue
    
    # get orthologue(s) in Pk
    ortho_hits <- ortho_PkPv$Pk_GeneID[ortho_PkPv$Pv_GeneID == gene_id]
    
    species <- "Pv"
    
  } else {
    return(list(list(warning = "Unrecognised input gene ID (must start with PVP01_... or PKNH_...)."))) #create a nested list that just contains the warning (as this is the output format of the function)
  } #doesn't continue
  
  # call existing helper function to get 'results' list (and filtering warnings) for input gene
  results_input <- get_gene_data(gene_id)
  
  #warning no orthologues:
  if (length(ortho_hits) == 0) {
    ortho_warn <- list(list(warning = "No orthologues found for the input gene")) #create a nested list that just contains the warning (like this it follows the format of results_input and they can be merged)
    return(c(ortho_warn, results_input)) #also still show input gene results list, because it might have warnings that we want to display.
  } #doesn't continue
  
  # call existing helper function to get 'results' list (and filtering warnings) for orthologous gene(s)
  results_ortho <- get_gene_data(ortho_hits)
  
  # Add similarity stats to orthologue results
  for (i in seq_along(results_ortho)) {
    results_i <- results_ortho[[i]]
    ortho_id <- results_i$gene
    
    if (species == "Pk") {
      sim_retrieve <- sim_PkPv[sim_PkPv$Pk_GeneID == gene_id & sim_PkPv$Pv_GeneID == ortho_id, ]
    }
    if (species == "Pv") {
      sim_retrieve <- sim_PkPv[sim_PkPv$Pv_GeneID == gene_id & sim_PkPv$Pk_GeneID == ortho_id, ] #empty (0 rows) if the gene is not in 'sim_PkPv' (when it's not passing filtering)
    }
    
    if (nrow(sim_retrieve) == 1) {
      results_i$dtw      <- sim_retrieve$dtw
      results_i$distance <- sim_retrieve$distance
      results_i$outcome  <- ifelse(sim_retrieve$outcome == "sim", "similar", "dissimilar") # sim_PkPv$outcome only contains 'sim' or 'dissim'
    } else {
      # no similarity row retrieved (ortho_id not in 'sim_PkPv', because it is filtered out)
      results_i$dtw      <- NA
      results_i$distance <- NA
      results_i$outcome  <- NA
    }
    
    results_ortho[[i]] <- results_i
  }
  
  results <- c(results_input, results_ortho)
  
  # If multiple orthologues found, add a top-level warning element
  if (length(ortho_hits) > 1) {
    multi_warn <- paste0("The input gene has ", length(ortho_hits), " orthologues: ", paste(ortho_hits, collapse = ", "))
    results <- c(list(list(warning = multi_warn)), results)
  }
  
  return(results)
  
  #structure 'results': nested list. Each element is a gene (input gene or orthologue). Each element itself is a sublist. The multi_warn warning can be added at the start when present, in an empty list element.
  #-list without name
  #    -warning: multi_warn warning in case there are >1 orthologues
  #-input gene list:
  #    -expression (df)
  #    -sampled points (df)
  #    -max (num)
  #    -gene (chr)
  #    -warning (chr)
  #- orthologue gene list:
  #    -expression (df)
  #    -sampled points (df)
  #    -max (num)
  #    -gene (chr)
  #    -warning (chr)
  #    -dtw (num)
  #    -distance (num)
  #    -outcome (chr)
  #- optionally additional orthologue gene lists (if >1 orthologue)
}

# D. orthologue tab Pk-Pf: Extract data for a single Pk/Pf input gene and its Pf/Pk orthologue(s).
get_ortho_data_Pf <- function(gene_id) {
  # gene_id : single Pf or Pk gene ID (string)
  
  # check existence of the input gene:
  if (grepl("^PKNH_", gene_id)) {
    
    if (!(gene_id %in% rownames(Pk_unfiltered))) {
      return(list(list(warning = "Unrecognised input gene ID (typo/new gene)."))) #create a nested list that just contains the warning (as this is the output format of the function)
      #same gene names as PkPf_unfiltered
    } #doesn't continue
    
    # get orthologue(s) in Pf
    ortho_hits <- ortho_PkPf$Pf_GeneID[ortho_PkPf$Pk_GeneID == gene_id]
    
    species <- "Pk"
    
  } else if (grepl("^PF3D7_", gene_id)) {
    
    if (!(gene_id %in% rownames(Pf_unfiltered))) {
      return(list(list(warning = "Unrecognised input gene ID (typo/new gene)."))) #create a nested list that just contains the warning (as this is the output format of the function)
    } #doesn't continue
    
    # get orthologue(s) in Pk
    ortho_hits <- ortho_PkPf$Pk_GeneID[ortho_PkPf$Pf_GeneID == gene_id]
    
    species <- "Pf"
    
  } else {
    return(list(list(warning = "Unrecognised input gene ID (must start with PF3D7_... or PKNH_...)."))) #create a nested list that just contains the warning (as this is the output format of the function)
  } #doesn't continue
  
  # call existing helper function to get 'results' list (and filtering warnings) for input gene
  results_input <- get_gene_data(gene_id)
  
  #warning no orthologues:
  if (length(ortho_hits) == 0) {
    ortho_warn <- list(list(warning = "No orthologues found for the input gene")) #create a nested list that just contains the warning (like this it follows the format of results_input and they can be merged)
    return(c(ortho_warn, results_input)) #also still show input gene results list, because it might have warnings that we want to display.
  } #doesn't continue
  
  # call existing helper function to get 'results' list (and filtering warnings) for orthologous gene(s)
  results_ortho <- get_gene_data(ortho_hits)
  
  # Add similarity stats to orthologue results
  for (i in seq_along(results_ortho)) {
    results_i <- results_ortho[[i]]
    ortho_id <- results_i$gene
    
    if (species == "Pk") {
      sim_retrieve <- sim_PkPf[sim_PkPf$Pk_GeneID == gene_id & sim_PkPf$Pf_GeneID == ortho_id, ]
    }
    if (species == "Pf") {
      sim_retrieve <- sim_PkPf[sim_PkPf$Pf_GeneID == gene_id & sim_PkPf$Pk_GeneID == ortho_id, ] #empty (0 rows) if the gene is not in 'sim_PkPf' (when it's not passing filtering)
    }
    
    if (nrow(sim_retrieve) == 1) {
      results_i$dtw      <- sim_retrieve$dtw
      results_i$distance <- sim_retrieve$distance
      results_i$outcome  <- ifelse(sim_retrieve$outcome == "sim", "similar", "dissimilar") # sim_PkPf$outcome only contains 'sim' or 'dissim'
    } else {
      # no similarity row retrieved (ortho_id not in 'sim_PkPf', because it is filtered out)
      results_i$dtw      <- NA
      results_i$distance <- NA
      results_i$outcome  <- NA
    }
    
    results_ortho[[i]] <- results_i
  }
  
  results <- c(results_input, results_ortho)
  
  # If multiple orthologues found, add a top-level warning element
  if (length(ortho_hits) > 1) {
    multi_warn <- paste0("The input gene has ", length(ortho_hits), " orthologues: ", paste(ortho_hits, collapse = ", "))
    results <- c(list(list(warning = multi_warn)), results)
  }
  
  return(results)
  
  #structure 'results': nested list. For structure, see at end of function get_ortho_data_Pv
}


# E. orthologue tab Pk-Pv/Pk-Pf: generate plot for every orthologous pair separately (>1 orthologous pair if the single input gene had multiple orthologues)
make_ortho_pair_plots <- function(results, gene_id, add_citation = FALSE) {
  # makes orthologous pairs, makes a plot for each pair, and returns a list of ggplot objects (one per input–orthologue pair)
  # uses results nested list made in function C.
  # add_citation default is FALSE.
  
  # keep only entries without warning
  valid <- Filter(function(x) is.null(x$warning), results) #go through list results, and only keep genes with no warning; throw out multi_warn warning if >1 orthologues (see C.).
  
  # if input gene was removed by warnings, nothing to plot
  if (!(gene_id %in% names(valid))) return(NULL)
  
  # split into input gene and orthologue(s)
  input_res <- valid[[gene_id]]
  ortho_res <- valid[names(valid) != gene_id]
  
  # if no orthologues remain (all filtered out), nothing to plot
  if (length(ortho_res) == 0) return(NULL)
  
  #now pair input gene with (each) orthologue, and generate plots
  plots <- list()
  for (ortho_id in names(ortho_res)) {
    pair <- list(input_res, ortho_res[[ortho_id]])
    p <- make_plot(pair, add_citation = add_citation) #what is entered when calling the make_ortho_pair_plots function (F or T) is passed here, and taken into the make_plot function.
    plots[[paste(gene_id, ortho_id, sep = "_vs_")]] <- p
  }
  
  # if no plots were created, return NULL
  if (length(plots) == 0) return(NULL)
  
  return(plots)
  #list of ggplot objects (one per input–orthologue pair)
}



#### 1. User interface object ####
UI <- fluidPage(
  # Add global CSS
  tags$head(
    tags$style(HTML("
      body { padding-bottom: 60px; }
      h4.section-title {
        font-size: 1.2em;
        font-weight: 600;
        margin-top: 25px;
        margin-bottom: 10px;
      }
    "))
  ),
  
  titlePanel(HTML("<i>P. knowlesi</i> Orthologue Expression Viewer")),
  
  tabsetPanel(
    tabPanel("Orthologue Expression Viewer Pk-Pv",
             sidebarLayout(
               sidebarPanel(
                 textInput("gene_id_PkPv", HTML("Visualise the expression of a <i>P. knowlesi</i> or <i>P. vivax</i> gene and its orthologue(s) over the asexual blood cycle.<br><br>
                                           Enter a <i>P. vivax</i> (PVP01_) or <i>P. knowlesi</i> (PKNH_) gene ID:"), 
                           "PVP01_0109300"), #default value: PvCRT
                 actionButton("go_ortho_PkPv", "Generate plot"),
                 br(), br(),
                 verbatimTextOutput("warningText_ortho_PkPv"),
                 uiOutput("summaryBox_PkPv")
               ),
               mainPanel(plotOutput("genePlot_ortho_PkPv"))
             )
    ),
    
    tabPanel("Orthologue Expression Viewer Pk-Pf",
             sidebarLayout(
               sidebarPanel(
                 textInput("gene_id_PkPf", HTML("Visualise the expression of a <i>P. knowlesi</i> or <i>P. falciparum</i> gene and its orthologue(s) over the asexual blood cycle.<br><br>
                                           Enter a <i>P. falciparum</i> (PF3D7_) or <i>P. knowlesi</i> (PKNH_) gene ID:"), 
                           "PF3D7_0709000"), #default value: PfCRT
                 actionButton("go_ortho_PkPf", "Generate plot"),
                 br(), br(),
                 verbatimTextOutput("warningText_ortho_PkPf"),
                 uiOutput("summaryBox_PkPf")
               ),
               mainPanel(plotOutput("genePlot_ortho_PkPf"))
             )
    ),
    
    tabPanel("Any Gene Expression Viewer",
             sidebarLayout(
               sidebarPanel(
                 textAreaInput("multi_genes", HTML("Visualise the expression of <i>P. knowlesi</i>, <i>P. vivax</i> and/or <i>P. falciparum</i> gene(s) over the asexual blood cycle.<br><br>
                                                  Enter up to 4 <i>P. knowlesi</i> (PKNH_), <i>P. vivax</i> (PVP01_) or <i>P. falciparum</i> (PF3D7_) gene IDs:<br>
                                                  <span style='font-size: 0.9em; color: gray;'>(comma-separated, combination of <i>P. vivax</i>, <i>P. knowlesi</i> and <i>P. falciparum</i> genes is allowed)"
                 ),
                 "PVP01_0623800,PKNH_1230100,PVP01_0934200,PF3D7_0731500", #default values: PvDBP, PkRON2, PvAMA1, PfEBA175
                 rows=3, resize="none"),
                 actionButton("go_multi", "Generate plot"),
                 br(), br(),
                 verbatimTextOutput("warningText_multi"),
                 uiOutput("download_ui")
               ),
               mainPanel(plotOutput("genePlot_multi"))
             )
    ),
    
    tabPanel("Information",
             fluidRow(
               column(10, offset = 1,
                      HTML("
                    <br>
                    <p>This tab contains additional information on how data was obtained and processed.</p>

                    <p>The code for this Shiny App can be found on this Github repository: https://github.com/kdemeulenaere/Pk_Orthologue_Expression_Viewer.</p>

                    <p><b>For more methodological details, see:</b><br>
                    De Meulenaere et al., 2026, doi: 10.64898/2026.02.09.704036.</p>

                    <h4 class='section-title'>Sample and data origins:</h4>
                    <ul>
                      <li><b><i>P. knowlesi</i>:</b> De Meulenaere et al., 2026, doi: 10.64898/2026.02.09.704036. Line A1-H.1 was synchronised and collected at 5 IDC time points (5, 14, 20, 24, 27 hpi). This line was grown <i>in vitro</i> in human erythrocytes, and does not produce gametocytes. The 5 hpi time point can show minor schizont contamination.</li>
                      <li><b><i>P. vivax</i>:</b> Zhu et al., 2016, Scientific reports. doi: 10.1038/srep20498. Clinical isolate smru1 contained mainly ring stages, was then matured <i>ex vivo</i>, and collected at 7 IDC time points (6, 18, 24, 30, 36, 42, 48 hpi). Since smru1 originates from a patient, it contains gametocytes.</li>
                      <li><b><i>P. falciparum</i>:</b> Subudhi et al., 2020, Nature Communications. doi: 10.1038/s41467-020-16593-y. Strain II3 (derived from 3D7) was synchronised and collected at 25 IDC time points (0-48 hpi, collected every 2 hours). The 2 replicates per time point were averaged after normalisation. This line does produce gametocytes.</li>
                    </ul>

                    <h4 class='section-title'>Plot:</h4>
                    <ul>
                      <li><b>Data:</b> Transcriptome data was collected at 5 (<i>P. knowlesi</i>), 7 (<i>P. vivax</i>), or 25 (<i>P. falciparum</i>) time points of the IDC. These time points are shown as dots. To visualise the full time range, expression was estimated by smoothing (LOWESS, only done for <i>P. falciparum</i>) followed by interpolation between the time points (PCHIP, done for all species). This interpolated data is shown as a line. No data is shown between time point 0 and the first sampled time point because no interpolation is possible here.</li>
                      <li><b>X-axis:</b> Expressed in relative time, since <i>P. knowlesi</i> A1-H.1 has a 27-hour IDC, while <i>P. vivax</i> and <i>P. falciparum</i> have a 48-hour IDC. Full IDC length was normalised to a 0–1 scale (<i>P. knowlesi:</i> /27; <i>P. vivax</i> and <i>P. falciparum:</i> /48).
                      <div style='display: flex; gap: 70px; margin-top: 8px; margin-bottom: 8px; align-items: flex-start;'>
  
                            <!-- Pk -->
                            <div>
                            <table style='border-collapse: collapse; font-size: 0.9em; color: gray; margin-top: 4px;'>
                            <tr>
                            <th colspan='2' style='text-align:left; padding:4px 12px; color: gray;'>
                            <i>P. knowlesi</i>
                            </th>
                            </tr>
                            <tr>
                            <th style='text-align: left; padding: 4px 12px;'>IDC time (hpi)</th>
                            <th style='text-align: left; padding: 4px 12px;'>Relative time</th>
                            </tr>
                            <tr><td style='padding: 4px 12px;'>5</td><td style='padding: 4px 12px;'>0.185</td></tr>
                            <tr><td style='padding: 4px 12px;'>14</td><td style='padding: 4px 12px;'>0.520</td></tr>
                            <tr><td style='padding: 4px 12px;'>20</td><td style='padding: 4px 12px;'>0.740</td></tr>
                            <tr><td style='padding: 4px 12px;'>24</td><td style='padding: 4px 12px;'>0.890</td></tr>
                            <tr><td style='padding: 4px 12px;'>27</td><td style='padding: 4px 12px;'>1.000</td></tr>
                            </table>
                            </div>
                            
                            <!-- Pv -->
                            <div>
                            <table style='border-collapse: collapse; font-size: 0.9em; color: gray; margin-top: 4px;'>
                            <tr>
                            <th colspan='2' style='text-align:left; padding:4px 12px; color: gray;'>
                            <i>P. vivax</i>
                            </th>
                            </tr>
                            <tr>
                            <th style='text-align: left; padding: 4px 12px;'>IDC time (hpi)</th>
                            <th style='text-align: left; padding: 4px 12px;'>Relative time</th>
                            </tr>
                            <tr><td style='padding: 4px 12px;'>6</td><td style='padding: 4px 12px;'>0.125</td></tr>
                            <tr><td style='padding: 4px 12px;'>18</td><td style='padding: 4px 12px;'>0.375</td></tr>
                            <tr><td style='padding: 4px 12px;'>24</td><td style='padding: 4px 12px;'>0.500</td></tr>
                            <tr><td style='padding: 4px 12px;'>30</td><td style='padding: 4px 12px;'>0.625</td></tr>
                            <tr><td style='padding: 4px 12px;'>36</td><td style='padding: 4px 12px;'>0.750</td></tr>
                            <tr><td style='padding: 4px 12px;'>42</td><td style='padding: 4px 12px;'>0.875</td></tr>
                            <tr><td style='padding: 4px 12px;'>48</td><td style='padding: 4px 12px;'>1.000</td></tr>
                            </table>
                            </div>
                            
                            <!-- Pf -->
                            <div>
                            <div style='font-size: 0.9em; color: gray; padding:4px 12px;'>
                            <b><i>P. falciparum</i></b>
                            </div>
                            
                            <div style='display: flex; gap: 20px;'>
                            
                            <!-- Pf col1 -->
                            <table style='border-collapse: collapse; font-size: 0.9em; color: gray; margin-top: 4px;'>
                            <tr>
                            <th style='text-align: left; padding: 4px 8px;'>IDC time (hpi)</th>
                            <th style='text-align: left; padding: 4px 8px;'>Relative time</th>
                            </tr>
                            <tr><td style='padding: 4px 8px;'>0</td><td>0.000</td></tr>
                            <tr><td style='padding: 4px 8px;'>2</td><td>0.040</td></tr>
                            <tr><td style='padding: 4px 8px;'>4</td><td>0.085</td></tr>
                            <tr><td style='padding: 4px 8px;'>6</td><td>0.125</td></tr>
                            <tr><td style='padding: 4px 8px;'>8</td><td>0.165</td></tr>
                            <tr><td style='padding: 4px 8px;'>10</td><td>0.210</td></tr>
                            <tr><td style='padding: 4px 8px;'>12</td><td>0.250</td></tr>
                            <tr><td style='padding: 4px 8px;'>14</td><td>0.290</td></tr>
                            <tr><td style='padding: 4px 8px;'>16</td><td>0.335</td></tr>
                            </table>
                            
                            <!-- Pf col2 -->
                            <table style='border-collapse: collapse; font-size: 0.9em; color: gray; margin-top: 4px;'>
                            <tr>
                            <th style='text-align: left; padding: 4px 8px;'>IDC time (hpi)</th>
                            <th style='text-align: left; padding: 4px 8px;'>Relative time</th>
                            </tr>
                            <tr><td style='padding: 4px 8px;'>18</td><td>0.375</td></tr>
                            <tr><td style='padding: 4px 8px;'>20</td><td>0.415</td></tr>
                            <tr><td style='padding: 4px 8px;'>22</td><td>0.460</td></tr>
                            <tr><td style='padding: 4px 8px;'>24</td><td>0.500</td></tr>
                            <tr><td style='padding: 4px 8px;'>26</td><td>0.540</td></tr>
                            <tr><td style='padding: 4px 8px;'>28</td><td>0.585</td></tr>
                            <tr><td style='padding: 4px 8px;'>30</td><td>0.625</td></tr>
                            <tr><td style='padding: 4px 8px;'>32</td><td>0.665</td></tr>
                            </table>
                            
                            <!-- Pf col3 -->
                            <table style='border-collapse: collapse; font-size: 0.9em; color: gray; margin-top: 4px;'>
                            <tr>
                            <th style='text-align: left; padding: 4px 8px;'>IDC time (hpi)</th>
                            <th style='text-align: left; padding: 4px 8px;'>Relative time</th>
                            </tr>
                            <tr><td style='padding: 4px 8px;'>34</td><td>0.710</td></tr>
                            <tr><td style='padding: 4px 8px;'>36</td><td>0.750</td></tr>
                            <tr><td style='padding: 4px 8px;'>38</td><td>0.790</td></tr>
                            <tr><td style='padding: 4px 8px;'>40</td><td>0.835</td></tr>
                            <tr><td style='padding: 4px 8px;'>42</td><td>0.875</td></tr>
                            <tr><td style='padding: 4px 8px;'>44</td><td>0.915</td></tr>
                            <tr><td style='padding: 4px 8px;'>46</td><td>0.960</td></tr>
                            <tr><td style='padding: 4px 8px;'>48</td><td>1.000</td></tr>
                            </table>
                          </div>
                        </div>
                        </li>
                      <li><b>Y-axis:</b> Transcript counts were TPM-normalised, log<sub>2</sub>-transformed (with a pseudocount of 1), then Z-scored (different genes can be shown on the same y-axis range).</li>
                      <small style='color: gray;'>Transcriptional <i>timing</i> can be compared between genes, expression <i>levels</i> have no meaning due to Z-scoring. <i>P. knowlesi</i> expression levels (TPM-normalised) can be found in Supplementary File S4 of De Meulenaere et al. (2026).</small></li>
                      <li><b>Dotted vertical line:</b> Indicates the time point of peak expression (maximal gene expression).<br>
                      <small style='color: gray;'>If multiple genes have the same peak expression time point, the dotted vertical lines are slightly shifted to prevent an overlap.</small></li>
                    </ul>

                    <h4 class='section-title'>Orthologue summary (in Orthologue Expression Viewer):</h4>
                    <ul>
                      <li><b>Orthologue ID:</b> A protein BLAST from <i>P. vivax</i> to <i>P. knowlesi</i> or from <i>P. falciparum</i> to <i>P. knowlesi</i> was carried out under the default settings. BLAST-resulting protein pairs with a % identity >=45% and with >=50% of the query protein aligned to its hit protein, were considered orthologues.<br>
                      <small style='color: gray;'>For ApiAP2 transcription factors, ApiAP2 genes with shared PlasmoDB gene names, a PlasmoDB OrthoMCL synteny hit or an indicated orthology in Jeninga et al. (2019, doi: 10.3390/pathogens8020047) were considered orthologues as well, in case there was no hit among the BLAST-based orthologue list.</small></li>
                      <li><b>Dynamic time warping score</b> Dynamic time warping (DTW; Sakoe-Chiba window type and temporal shift of 0.2 relative time units allowed) was applied to the expression levels over the time interval that was sampled in both species (relative time 0.185 – 1). The lower the score, the more similar the orthologue expression profiles are.</li>
                      <li><b>Time shift between maxima:</b> The difference (in relative time units) between the peak expression time points (maximal expression) of both orthologues. Circularity of the IDC is taken into account (for example, the time difference between relative times 0.9 and 0.1 would be 0.2).</li>
                      <li><b>Classification:</b> Expression patterns of orthologous gene pairs were classified as similar or dissimilar, based on the DTW score and time shift between the maxima. A full explanation can be found in Figure 9 of De Meulenaere et al. (2026).</li>
                    </ul>

                    <h4 class='section-title'>Warnings:</h4>
                    <ul>
                      <li><b>Gene has too low transcriptional variation over IDC (log<sub>2</sub> fold change between highest and lowest normalised expression level &lt; 0.5)</b>: genes with low transcriptional variation are not suitable for Z-scoring (used for the y-axis and orthologue comparison statistics), as this artificially amplifies minor fluctuations in expression levels. Transcriptional variation was quantified as the log<sub>2</sub> fold change between the minimal and maximal expression level (pseudocount of 1).</li>
                      <li><b>Gene has too low expression (sum of normalised expression levels of sampled time points &lt; 0.5)</b>: genes with virtually no expression are not shown.</li>
                      <li><b>The input gene has multiple orthologues:</b> when multiple orthologues are found for a single input gene ID, a plot and orthologue summary is generated for each of them.</li>
                      <li><b>Unrecognised gene ID (must start with PVP01_..., PKNH_... or PF3D7_...):</b> the gene ID is not from the correct reference genome, or the prefix was entered incorrectly (<i>e.g.</i> should be uppercase).</li>
                      <li><b>Unrecognised gene ID (typo/new gene):</b> the gene ID was not detected in the dataset (typo, or a newer gene ID absent in the used PlasmoDB reference genomes: strain H v2 release 66, PVP01 v2 release 66, 3D7 v3 release 67).</li>
                      <li><b>You entered more than 4 gene IDs. Only the first 4 are processed and plotted</b>: the 'Any Gene tab' input box only accepts up to 4 gene ID's.</li>
                    </ul>

                    <h4 class='section-title'>Abbreviations:</h4>
                    <ul>
                      <li><b>DTW:</b> dynamic time warping</li>
                      <li><b>hpi:</b> hours post-invasion</li>
                      <li><b>IDC:</b> intraerythrocytic developmental cycle</li>
                      <li><b>LOWESS:</b> Locally Weighted Scatterplot Smoothing</li>
                      <li><b>PCHIP:</b> Piecewise Cubic Hermite Interpolating Polynomial</li>
                      
                    </ul>
                  "),
                      
                      br(), br()
               )
             )
    )
  ),
  
  # Sticky footer citation + app version
  div(
    id = "app-footer",
    div(
      HTML("<span style='margin-right: 4em;'>v2</span>If you use this tool, please cite: De Meulenaere et al., 2026, doi: 10.64898/2026.02.09.704036."),
      style = "font-size: 0.9em; color: gray; text-align: left;
               position: fixed; bottom: 0; left: 0; width: 100%;
               background: white; padding: 5px 20px; z-index: 1000;
               border-top: 1px solid #ccc;"
    )
  )
)

#### 2. Server function ####
Server <- function(input, output) {
  
  # ==== Orthologue Pk-Pv tab ====
  
  # Retrieve data for plotting and summary box
  data_ortho_PkPv <- eventReactive(input$go_ortho_PkPv, {
    get_ortho_data_Pv(input$gene_id_PkPv)
  }) #run if user presses go_ortho_PkPv button
  
  # Plot
  output$genePlot_ortho_PkPv <- renderPlot({
    dlist <- data_ortho_PkPv()
    if (is.null(dlist)) return(NULL) #if no output from get_ortho_data_Pv, return NULL and do not plot (non-existing gene entered)
    
    plots_list <- make_ortho_pair_plots(dlist, input$gene_id_PkPv)
    if (is.null(plots_list)) return(NULL) #if nothing could be plotted and make_ortho_pair_plots returned NULL.
    
    # plots_list is a list of ggplot objects; arrange them vertically
    gridExtra::grid.arrange(grobs = plots_list, ncol = 1)
  }, height = function() {
    dlist <- data_ortho_PkPv()
    if (is.null(dlist)) return(600)
    
    # calculate number of valid orthologs (has dtw in list indicative of an ortholog, no warning)
    n_plots <- sum(vapply(dlist, function(x) is.null(x$warning) && !is.null(x$dtw), logical(1))) #counts if TRUE
    
    if (n_plots == 0) {400} else {
      max(600, 300 * n_plots)  # 300px per plot, at least 600px
    }
  })
  
  # Similarity stats + download button
  output$summaryBox_PkPv <- renderUI({
    dlist <- data_ortho_PkPv()
    if (is.null(dlist)) return(NULL)
    
    # keep only entries without warnings, and keep only orthologues (which have a dtw element in the list)
    valid_ortho <- Filter(function(x) is.null(x$warning) && !is.null(x$dtw), dlist)
    if (length(valid_ortho) == 0) return(NULL)
    
    plots_list <- make_ortho_pair_plots(dlist, input$gene_id_PkPv, add_citation = TRUE)
    if (is.null(plots_list)) return(NULL)
    
    tagList(lapply(seq_along(valid_ortho), function(i) {
      dd <- valid_ortho[[i]]
      download_id <- paste0("download_plot_PkPv_", i)  # unique ID
      
      # Create download handler dynamically
      output[[download_id]] <- downloadHandler(
        filename = function() paste0("Expression_",input$gene_id_PkPv,"_", dd$gene, ".png"),
        content = function(file) {
          plot_name <- paste(input$gene_id_PkPv, dd$gene, sep = "_vs_")
          p <- plots_list[[plot_name]]
          ggsave(file, p, width = 10, height = 6, dpi=300)
        }
      )
      
      tags$div(
        style = "border: 1px solid #ccc; margin: 10px; padding: 10px;",
        tags$h4("Orthologue summary"),
        tags$p(strong("Orthologue ID: "), dd$gene),
        tags$p(strong("Dynamic time warping score: "), round(dd$dtw, 3)),
        tags$p(strong("Time shift between maxima: "), dd$distance),
        tags$p(strong("Classification: "), dd$outcome),
        downloadButton(download_id, "Download plot")
      )
    }))
  })
  
  # Warnings
  output$warningText_ortho_PkPv <- renderText({
    dlist <- data_ortho_PkPv()
    if (is.null(dlist)) return(NULL)
    
    # extract warnings for each gene
    ws <- vapply(dlist, function(x) if (!is.null(x$warning)) as.character(x$warning) else NA_character_, FUN.VALUE = character(1))
    ws <- ws[!is.na(ws)]
    
    if (length(ws) == 0) return(NULL)
    paste(ws, collapse = "\n\n")
  })
  
  # ==== Orthologue Pk-Pf tab ====
  
  # Retrieve data for plotting and summary box
  data_ortho_PkPf <- eventReactive(input$go_ortho_PkPf, {
    get_ortho_data_Pf(input$gene_id_PkPf)
  }) #run if user presses go_ortho_PkPf button
  
  # Plot
  output$genePlot_ortho_PkPf <- renderPlot({
    dlist <- data_ortho_PkPf()
    if (is.null(dlist)) return(NULL) #if no output from get_ortho_data_Pv, return NULL and do not plot (non-existing gene entered)
    
    plots_list <- make_ortho_pair_plots(dlist, input$gene_id_PkPf)
    if (is.null(plots_list)) return(NULL) #if nothing could be plotted and make_ortho_pair_plots returned NULL.
    
    # plots_list is a list of ggplot objects; arrange them vertically
    gridExtra::grid.arrange(grobs = plots_list, ncol = 1)
  }, height = function() {
    dlist <- data_ortho_PkPf()
    if (is.null(dlist)) return(600)
    
    # calculate number of valid orthologs (has dtw in list indicative of an ortholog, no warning)
    n_plots <- sum(vapply(dlist, function(x) is.null(x$warning) && !is.null(x$dtw), logical(1))) #counts if TRUE
    
    if (n_plots == 0) {400} else {
      max(600, 300 * n_plots)  # 300px per plot, at least 600px
    }
  })
  
  # Similarity stats + download button
  output$summaryBox_PkPf <- renderUI({
    dlist <- data_ortho_PkPf()
    if (is.null(dlist)) return(NULL)
    
    # keep only entries without warnings, and keep only orthologues (which have a dtw element in the list)
    valid_ortho <- Filter(function(x) is.null(x$warning) && !is.null(x$dtw), dlist)
    if (length(valid_ortho) == 0) return(NULL)
    
    plots_list <- make_ortho_pair_plots(dlist, input$gene_id_PkPf, add_citation = TRUE)
    if (is.null(plots_list)) return(NULL)
    
    tagList(lapply(seq_along(valid_ortho), function(i) {
      dd <- valid_ortho[[i]]
      download_id <- paste0("download_plot_PkPf_", i)  # unique ID
      
      # Create download handler dynamically
      output[[download_id]] <- downloadHandler(
        filename = function() paste0("Expression_",input$gene_id_PkPf,"_", dd$gene, ".png"),
        content = function(file) {
          plot_name <- paste(input$gene_id_PkPf, dd$gene, sep = "_vs_")
          p <- plots_list[[plot_name]]
          ggsave(file, p, width = 10, height = 6, dpi=300)
        }
      )
      
      tags$div(
        style = "border: 1px solid #ccc; margin: 10px; padding: 10px;",
        tags$h4("Orthologue summary"),
        tags$p(strong("Orthologue ID: "), dd$gene),
        tags$p(strong("Dynamic time warping score: "), round(dd$dtw, 3)),
        tags$p(strong("Time shift between maxima: "), dd$distance),
        tags$p(strong("Classification: "), dd$outcome),
        downloadButton(download_id, "Download plot")
      )
    }))
  })
  
  # Warnings
  output$warningText_ortho_PkPf <- renderText({
    dlist <- data_ortho_PkPf()
    if (is.null(dlist)) return(NULL)
    
    # extract warnings for each gene
    ws <- vapply(dlist, function(x) if (!is.null(x$warning)) as.character(x$warning) else NA_character_, FUN.VALUE = character(1))
    ws <- ws[!is.na(ws)]
    
    if (length(ws) == 0) return(NULL)
    paste(ws, collapse = "\n\n")
  })
  
  # ==== Any gene tab ====
  
  # Retrieve data for plotting
  data_multi <- eventReactive(input$go_multi, {
    ids <- strsplit(input$multi_genes, ",")[[1]]
    ids <- trimws(ids)  # remove spaces if they were entered
    ids <- ids[ids != ""]  # drop empties
    TooMany_warning <- NULL
    if (length(ids) > 4) {
      TooMany_warning <- "You entered more than 4 gene IDs. Only the first 4 are processed and plotted."
    }
    ids <- ids[1:min(4, length(ids))]  # max 4 genes can be entered. If more, only the first 4 are taken.
    res <- get_gene_data(ids)
    if (!is.null(TooMany_warning)) {
      # attach warning as an extra element in the list
      res[[".TooMany_warning"]] <- list(warning = TooMany_warning)
    }
    res
  })
  
  # Retrieve plot
  plot_multi <- reactive({
    dlist <- data_multi()
    if (is.null(dlist)) return(NULL)
    
    p <- make_plot(dlist)
    if (is.null(p)) return(NULL)
    p
  })
  
  # Plot
  output$genePlot_multi <- renderPlot({
    plot_multi()
  }, height = 600)
  
  # Warnings
  output$warningText_multi <- renderText({
    dlist <- data_multi()
    if (is.null(dlist)) return(NULL)
    
    # extract warnings for each gene
    ws <- vapply(dlist, function(x) if (!is.null(x$warning)) as.character(x$warning) else NA_character_, FUN.VALUE = character(1))
    ws <- ws[!is.na(ws)]
    
    if (length(ws) == 0) return(NULL)
    paste(ws, collapse = "\n\n")
  })
  
  # Download button
  output$download_ui <- renderUI({
    req(plot_multi())  # only proceed and show download button if plot_multi() is not NULL
    downloadButton("download_plot_multi", "Download plot") #only shows the button if plot exists (otherwise returns NULL)
  })
  
  # Download handler
  output$download_plot_multi <- downloadHandler(
    filename = function() {
      paste0("multi_gene_expression.png")
    },
    content = function(file) {
      dlist <- data_multi()
      if (is.null(dlist)) return(NULL)
      p <- make_plot(dlist, add_citation = TRUE)
      ggsave(file, p, width = 10, height = 6, dpi = 300)
    }
  )
  
}


#### 3. Call to shinyApp function ####
shinyApp(ui = UI, server = Server)




