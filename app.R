#### P. knowlesi - P. vivax Expression Viewer ####

# Author: Katlijn De Meulenaere

# Script to generate a Shiny App for investigation of P. knowlesi (A1-H.1; De Meulenaere et al., 2025) and P. vivax (smru1 clinical isolate; Zhu et al., 2016) gene expression patterns over the IDC (asexual blood stages).
# In the first tab of the App, a Pv/Pk gene is entered, the orthologue(s) are searched, and an expression plot for the input gene + orthologue(s) is given together with similarity statistics of the expression patterns.
# In the second tab of the App, up to 4 Pk/Pv genes are entered, and one expression plot is given for all input genes together.

# To launch the App, do:
# library(shiny)
# runApp("path/Pk_Pv_comparator") #app.R is stored in /Pk_Pv_comparator directory.


#### Load packages ####

library(shiny) #v1.9.1
library(ggplot2) #v3.5.2
library(gridExtra) #v2.3
library(cowplot) #v1.1.3
# R version 4.3.0

#### Import dataframes ####

Pk_expr <- read.delim("./datasets/Pk_expr.txt", check.names = FALSE)
Pv_expr <- read.delim("./datasets/Pv_expr.txt", check.names = FALSE)
#Expression of all gene ID's that passed filtering (expression: TPM, log2 with a pseudocount of 1, Z-scored)). Genes without orthologue are present.

sim <- read.delim("./datasets/sim.txt")
#orthologous gene ID's (that passed filtering) with their similarity stats: dtw, time distance, outcome (sim/dissim)

ortho <- read.delim("./datasets/ortho.txt")
#orthologous gene ID's (no filtering applied)

Pk_unfiltered <- read.delim("./datasets/Pk_unfiltered.txt")
Pv_unfiltered <- read.delim("./datasets/Pv_unfiltered.txt")
#dataframe with all Pv/Pk genes, and whether they pass the filtering rules.


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
          sampled_points <- df_expr[df_expr$reltime %in% c(0.185, 0.52, 0.74, 0.89, 1), ]
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
          sampled_points <- df_expr[df_expr$reltime %in% c(0.125, 0.375, 0.5, 0.625, 0.75, 0.875, 1), ]
        }
      }
      
    } else {
      warnings <- c(warnings, paste0("Unrecognised gene ID: ", gene_id, " (must start with PVP01_... or PKNH_...)."))
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
    geom_vline(data = df_max, aes(xintercept = reltime_offset, color = gene, linetype = "maximal expression of gene"),
               linewidth = 1, show.legend = c(color = FALSE)) +
    scale_x_continuous(limits = c(0, 1.01)) + #>1 for the potential offsets
    scale_color_manual(
      values = c("#D81B60", "#1E88E5", "#FFC107", "#004D40")) +
    scale_linetype_manual(
      name = NULL,
      values = c("maximal expression of gene" = "dotted")) +
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
      cowplot::draw_label(
        "De Meulenaere et al., 2025, XXX. doi: XXX.",
        x = 0, y = 0.5, hjust = 0, vjust = 0.5,
        size = 10, color = "gray40"
      )
    p <- cowplot::plot_grid(p, citation, ncol = 1, rel_heights = c(1, 0.1))
  }
  
  return(p)
}

# C. orthologue tab: Extract data for a single input gene and its orthologue(s).
get_ortho_data <- function(gene_id) {
  # gene_id : single Pv or Pk gene ID (string)
  
  # check existence of the input gene:
  if (grepl("^PKNH_", gene_id)) {
    
    if (!(gene_id %in% rownames(Pk_unfiltered))) {
      return(list(list(warning = "Unrecognised input gene ID (typo/new gene)."))) #create a nested list that just contains the warning (as this is the output format of the function)
    } #doesn't continue
    
    # get orthologue(s) in Pv
    ortho_hits <- ortho$Pv_GeneID[ortho$Pk_GeneID == gene_id]
    
    species <- "Pk"
    
  } else if (grepl("^PVP01_", gene_id)) {
    
    if (!(gene_id %in% rownames(Pv_unfiltered))) {
      return(list(list(warning = "Unrecognised input gene ID (typo/new gene)."))) #create a nested list that just contains the warning (as this is the output format of the function)
    } #doesn't continue
    
    # get orthologue(s) in Pk
    ortho_hits <- ortho$Pk_GeneID[ortho$Pv_GeneID == gene_id]
    
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
      sim_retrieve <- sim[sim$Pk_GeneID == gene_id & sim$Pv_GeneID == ortho_id, ]
    }
    if (species == "Pv") {
      sim_retrieve <- sim[sim$Pv_GeneID == gene_id & sim$Pk_GeneID == ortho_id, ] #empty (0 rows) if the gene is not in 'sim' (when it's not passing filtering)
    }
    
    if (nrow(sim_retrieve) == 1) {
      results_i$dtw      <- sim_retrieve$dtw
      results_i$distance <- sim_retrieve$distance
      results_i$outcome  <- ifelse(sim_retrieve$outcome == "sim", "similar", "dissimilar") # sim$outcome only contains 'sim' or 'dissim'
    } else {
      # no similarity row retrieved (ortho_id not in 'sim', because it is filtered out)
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


# D. orthologue tab: generate plot for every orthologous pair separately (>1 orthologous pair if the single input gene had multiple orthologues)
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
  
  titlePanel(HTML("<i>P. knowlesi</i> - <i>P. vivax</i> Expression Viewer")),
  
  tabsetPanel(
    tabPanel("Orthologue Expression Viewer",
             sidebarLayout(
               sidebarPanel(
                 textInput("gene_id", HTML("Visualise the expression of a <i>P. knowlesi</i> or <i>P. vivax</i> gene and its orthologue(s) over the asexual blood cycle.<br><br>
                                           Enter a <i>P. vivax</i> (PVP01_) or <i>P. knowlesi</i> (PKNH_) gene ID:"), 
                           "PVP01_0109300"), #default value: PvCRT
                 actionButton("go_ortho", "Generate plot"),
                 br(), br(),
                 verbatimTextOutput("warningText_ortho"),
                 uiOutput("summaryBox")
               ),
               mainPanel(plotOutput("genePlot_ortho"))
             )
    ),
    
    tabPanel("Any Gene Expression Viewer",
             sidebarLayout(
               sidebarPanel(
                 textAreaInput("multi_genes", HTML("Visualise the expression of <i>P. knowlesi</i> and/or <i>P. vivax</i> gene(s) over the asexual blood cycle.<br><br>
                                                  Enter up to 4 <i>P. knowlesi</i> (PKNH_) or <i>P. vivax</i> (PVP01_) gene IDs:<br>
                                                  <span style='font-size: 0.9em; color: gray;'>(comma-separated, combination of <i>P. vivax</i> and <i>P. knowlesi</i> genes is allowed)"
                                                   ),
                               "PVP01_0623800,PKNH_1230100,PVP01_0934200", #default values: PvDBP, PkRON2, PvAMA1
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

                    <p>The code for this Shiny App can be found on: XXX ENTER GITHUB LINK XXX</p>

                    <p><b>For more methodological details, see:</b><br>
                    De Meulenaere et al., 2025, XXX. doi: XXX.</p>

                    <h4 class='section-title'>Sample and data origins:</h4>
                    <ul>
                      <li><b><i>P. knowlesi</i>:</b> De Meulenaere et al., 2025, XXX. doi: XXX. Line A1-H.1 was synchronised and collected at 5 IDC time points (5, 14, 20, 24, 27 hpi). This line was grown <i>in vitro</i> in human erythrocytes, and does not produce gametocytes. The 5 hpi time point can show minor schizont contamination.</li>
                      <li><b><i>P. vivax</i>:</b> Zhu et al., 2016, Scientific reports. doi: 10.1038/srep20498. Clinical isolate smru1 contained mainly ring stages, was then matured <i>ex vivo</i>, and collected at 7 IDC time points (6, 18, 24, 30, 36, 42, 48 hpi). Since smru1 originates from a patient, it contains gametocytes.</li>
                    </ul>

                    <h4 class='section-title'>Plot:</h4>
                    <ul>
                      <li><b>Data:</b> Transcriptome data was collected at 5 (<i>P. knowlesi</i>) or 7 (<i>P. vivax</i>) time points of the IDC. These time points are shown as dots. In between those collected time points, expression was estimated by interpolation (PCHIP). This interpolated data is shown as a line. No data is shown between time point 0 and the first sampled time point because no interpolation is possible here.</li>
                      <li><b>X-axis:</b> Expressed in relative time, since <i>P. knowlesi</i> A1-H.1 has a 27-hour IDC and <i>P. vivax</i> a 48-hour IDC. Full IDC length was normalised to a range of 0–1 (<i>P. knowlesi:</i> /27; <i>P. vivax:</i> /48).
                      <div style='display: flex; gap: 30px; margin-top: 8px; margin-bottom: 8px;'>
                          <div>
                            <table style='border-collapse: collapse; font-size: 0.9em; color: gray; margin-top: 4px;'>
                              <tr>
                                <th colspan='2' style='text-align:left; padding:4px 12px; color: gray;'>
                                  <i>P. knowlesi</i></th>
                              </tr>
                              <tr>
                                <th style='text-align: left; padding: 4px 12px;'>Developmental time (hpi)</th>
                                <th style='text-align: left; padding: 4px 12px;'>Relative time</th>
                              </tr>
                              <tr><td style='padding: 4px 12px;'>5</td><td style='padding: 4px 12px;'>0.185</td></tr>
                              <tr><td style='padding: 4px 12px;'>14</td><td style='padding: 4px 12px;'>0.52</td></tr>
                              <tr><td style='padding: 4px 12px;'>20</td><td style='padding: 4px 12px;'>0.74</td></tr>
                              <tr><td style='padding: 4px 12px;'>24</td><td style='padding: 4px 12px;'>0.89</td></tr>
                              <tr><td style='padding: 4px 12px;'>27</td><td style='padding: 4px 12px;'>1</td></tr>
                            </table>
                          </div>
                          <div>
                            <table style='border-collapse: collapse; font-size: 0.9em; color: gray; margin-top: 4px;'>
                              <tr>
                                <th colspan='2' style='text-align:left; padding:4px 12px; color: gray;'>
                                <i>P. vivax</i>
                              </tr>
                              <tr>
                                <th style='text-align: left; padding: 4px 12px;'>Developmental time (hpi)</th>
                                <th style='text-align: left; padding: 4px 12px;'>Relative time</th>
                              </tr>
                              <tr><td style='padding: 4px 12px;'>6</td><td style='padding: 4px 12px;'>0.125</td></tr>
                              <tr><td style='padding: 4px 12px;'>18</td><td style='padding: 4px 12px;'>0.375</td></tr>
                              <tr><td style='padding: 4px 12px;'>24</td><td style='padding: 4px 12px;'>0.5</td></tr>
                              <tr><td style='padding: 4px 12px;'>30</td><td style='padding: 4px 12px;'>0.625</td></tr>
                              <tr><td style='padding: 4px 12px;'>36</td><td style='padding: 4px 12px;'>0.75</td></tr>
                              <tr><td style='padding: 4px 12px;'>42</td><td style='padding: 4px 12px;'>0.875</td></tr>
                              <tr><td style='padding: 4px 12px;'>48</td><td style='padding: 4px 12px;'>1</td></tr>
                            </table>
                          </div>
                        </div>
                        </li>
                      <li><b>Y-axis:</b> Transcript counts were TPM-normalised, log<sub>2</sub>-transformed (with a pseudocount of 1), then Z-scored (to compare between <i>P. knowlesi</i> and <i>P. vivax</i>).</li>
                      <li><b>Dotted vertical line:</b> Indicates the time point of maximal gene expression.<br>
                      <small style='color: gray;'>If multiple genes have the same time point of maximal gene expression, the dotted vertical lines are slightly shifted to prevent an overlap.</small></li>
                    </ul>

                    <h4 class='section-title'>Orthologue summary (in Orthologue Expression Viewer):</h4>
                    <ul>
                      <li><b>Orthologue ID:</b> A protein BLAST from <i>P. vivax</i> to <i>P. knowlesi</i> was carried out under the default settings. BLAST-resulting protein pairs with a % identity >=45% and with >=50% of the query protein aligned to its hit protein, were considered orthologues.<br>
                      <small style='color: gray;'>For ApiAP2 transcription factors, ApiAP2 genes with shared PlasmoDB gene names, a PlasmoDB OrthoMCL synteny hit or an indicated orthology in Jeninga et al. (2019, doi: 10.3390/pathogens8020047) were considered orthologues as well, in case there was no hit among the BLAST-based orthologue list.</small></li>
                      <li><b>Dynamic time warping score</b> Dynamic time warping (DTW; Sakoe-Chiba window type and phase shift of 0.2 relative time units allowed) was applied to the expression levels over the time interval that was sampled in both species (relative time 0.185 – 1). The lower the score, the more similar the orthologue expression profiles are.</li>
                      <li><b>Time distance between maxima:</b> The distance (in relative time units) between the time points of maximal expression of both orthologues. Circularity of the IDC is taken into account (for example, the time distance between relative times 0.9 and 0.1 would be 0.2).</li>
                      <li><b>Classification:</b> Expression patterns of orthologous gene pairs were classified as similar or dissimilar, based on the DTW score and time distance between the maxima. A full explanation can be found in Figure 3 of De Meulenaere et al. (2025).</li>
                    </ul>

                    <h4 class='section-title'>Warnings:</h4>
                    <ul>
                      <li><b>Gene has too low transcriptional variation over IDC (log<sub>2</sub> fold change between highest and lowest normalised expression level &lt; 0.5)</b>: genes with low transcriptional variation are not suitable for Z-scoring (used for the y-axis and orthologue comparison statistics), as this artificially amplifies minor fluctuations in expression levels. Transcriptional variation was quantified as the fold change between the minimal and maximal expression levels (pseudocount of 1).</li>
                      <li><b>Gene has too low expression (sum of normalised expression levels of sampled time points &lt; 0.5)</b>: genes with virtually no expression are not shown.</li>
                      <li><b>The input gene has multiple orthologues:</b> when multiple orthologues are found for a single input gene ID, a plot and orthologue summary is generated for each of them.</li>
                      <li><b>Unrecognised gene ID (must start with PVP01_... or PKNH_...):</b> the gene ID is not from the correct reference genome, or the prefix was entered incorrectly (<i>e.g.</i> should be uppercase).</li>
                      <li><b>Unrecognised gene ID (typo/new gene):</b> the gene ID was not detected in the dataset (typo, or a newer gene ID absent in strain H annotation v2, PlasmoDB release 66).</li>
                      <li><b>You entered more than 4 gene IDs. Only the first 4 are processed and plotted</b>: the 'Any Gene tab' input box only accepts up to 4 gene ID's.</li>
                    </ul>

                    <h4 class='section-title'>Abbreviations:</h4>
                    <ul>
                      <li><b>IDC:</b> intraerythrocytic developmental cycle</li>
                      <li><b>hpi:</b> hours post-invasion</li>
                      <li><b>DTW:</b> dynamic time warping</li>
                      <li><b>PCHIP:</b> Piecewise Cubic Hermite Interpolating Polynomial</li>
                    </ul>
                  "),
                      
                      br(), br()
               )
             )
    )
  ),
  
  # Sticky footer citation
  div(
    id = "app-footer",
    div(
      "If you use this tool, please cite: De Meulenaere et al., 2025, XXX. doi: XXX.",
      style = "font-size: 0.9em; color: gray; text-align: center;
               position: fixed; bottom: 0; left: 0; width: 100%;
               background: white; padding: 5px 0; z-index: 1000;
               border-top: 1px solid #ccc;"
    )
  )
)

#### 2. Server function ####
Server <- function(input, output) {
  
  # ==== Orthologue tab ====
  
  # Retrieve data for plotting and summary box
  data_ortho <- eventReactive(input$go_ortho, {
    get_ortho_data(input$gene_id)
  }) #run if user presses go_ortho button
  
  # Plot
  output$genePlot_ortho <- renderPlot({
    dlist <- data_ortho()
    if (is.null(dlist)) return(NULL) #if no output from get_ortho_data, return NULL and do not plot (non-existing gene entered)
    
    plots_list <- make_ortho_pair_plots(dlist, input$gene_id)
    if (is.null(plots_list)) return(NULL) #if nothing could be plotted and make_ortho_pair_plots returned NULL.
    
    # plots_list is a list of ggplot objects; arrange them vertically
    gridExtra::grid.arrange(grobs = plots_list, ncol = 1)
  }, height = function() {
    dlist <- data_ortho()
    if (is.null(dlist)) return(600)
    
    # calculate number of valid orthologs (has dtw in list indicative of an ortholog, no warning)
    n_plots <- sum(vapply(dlist, function(x) is.null(x$warning) && !is.null(x$dtw), logical(1))) #counts if TRUE
    
    if (n_plots == 0) {400} else {
      max(600, 300 * n_plots)  # 300px per plot, at least 600px
    }
  })
  
  # Similarity stats + download button
  output$summaryBox <- renderUI({
    dlist <- data_ortho()
    if (is.null(dlist)) return(NULL)
    
    # keep only entries without warnings, and keep only orthologues (which have a dtw element in the list)
    valid_ortho <- Filter(function(x) is.null(x$warning) && !is.null(x$dtw), dlist)
    if (length(valid_ortho) == 0) return(NULL)
    
    plots_list <- make_ortho_pair_plots(dlist, input$gene_id, add_citation = TRUE)
    if (is.null(plots_list)) return(NULL)
    
    tagList(lapply(seq_along(valid_ortho), function(i) {
      dd <- valid_ortho[[i]]
      download_id <- paste0("download_plot_", i)  # unique ID
      
      # Create download handler dynamically
      output[[download_id]] <- downloadHandler(
        filename = function() paste0("Expression_",input$gene_id,"_", dd$gene, ".png"),
        content = function(file) {
          plot_name <- paste(input$gene_id, dd$gene, sep = "_vs_")
          p <- plots_list[[plot_name]]
          ggsave(file, p, width = 10, height = 6, dpi=300)
        }
      )
      
      tags$div(
        style = "border: 1px solid #ccc; margin: 10px; padding: 10px;",
        tags$h4("Orthologue summary"),
        tags$p(strong("Orthologue ID: "), dd$gene),
        tags$p(strong("Dynamic time warping score: "), round(dd$dtw, 3)),
        tags$p(strong("Time distance between maxima: "), dd$distance),
        tags$p(strong("Classification: "), dd$outcome),
        downloadButton(download_id, "Download plot")
      )
    }))
  })
  
  # Warnings
  output$warningText_ortho <- renderText({
    dlist <- data_ortho()
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




