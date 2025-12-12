library(tidyverse)
library(ggplot2)
library(googledrive)
library(googlesheets4)
library(data.table)

# Increase connection buffer size for large files
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 100)

##---------------------------------------------------------
## FUNCTION: plot_expression_and_amreads
##---------------------------------------------------------

plot_expression_and_amreads <- function(query_genes,
                                        title_prefix,
                                        trt_filter = c("M", "C"),
                                        time_filter = c("T1", "T2")) {
  tryCatch({
    ## Step 1: Load annotation table
    message("Loading annotation table...")
    DE_ASE_table <- read_csv("data/DE_ASE_table.csv", progress = FALSE)
    
    ## Step 2: Filter selected genes
    message("Filtering selected genes...")
    selected_genes <- DE_ASE_table %>%
      filter(geneID %in% query_genes | Working.Symbol %in% query_genes)
    
    if (nrow(selected_genes) == 0) {
      stop("No genes found matching the query. Please check your gene IDs/symbols.")
    }
    
    ## Step 3: Load normalized counts
    message("Loading normalized counts...")
    norm_counts_DE <- read_csv("data/norm_counts_DE.csv", progress = FALSE) %>%
      select(-any_of("...1")) %>%      # safer than select(-"...1")
      filter(geneID %in% selected_genes$geneID) %>%
      right_join(selected_genes, by = "geneID")
    
    if (nrow(norm_counts_DE) == 0) {
      stop("No rows in norm_counts_DE after filtering for selected genes.")
    }
    
    ## Step 4: Load Total AM reads file
    message("Loading AM reads data...")
    Total_AMreads_df <- read_csv("data/Tota_AMreads_ineachsample.csv",
                                 progress = FALSE)
    
    ## Step 5: Reshape and annotate metadata
    message("Processing data for plotting...")
    long <- norm_counts_DE %>%
      pivot_longer(
        cols = matches("_M_|_C_"),
        names_to = "Sample",
        values_to = "norm_counts"
      ) %>%
      mutate(label = Sample) %>%
      left_join(Total_AMreads_df, by = "label") %>%
      separate(Sample,
               into = c("genotype", "trt", "time", "rep"),
               sep = "_",
               remove = FALSE) %>%
      mutate(
        generation = case_when(
          genotype == "Oh43xB73" ~ "F1",
          genotype %in% c("Oh43", "B73") ~ "F0",
          TRUE ~ NA_character_
        ),
        genotype = factor(genotype,
                          levels = c("Oh43", "B73", "Oh43xB73")),
        generation = factor(generation,
                            levels = c("F0", "F1"))
      ) %>%
      filter(!str_detect(genotype, "W22"))
    
    # Apply treatment/time filters
    long <- long %>%
      filter(trt %in% trt_filter, time %in% time_filter)
    
    if (nrow(long) == 0) {
      stop("No data left after filtering by treatment/time. Try relaxing filters or checking gene IDs.")
    }
    
    ## Step 6: Summary statistics
    message("Calculating summary statistics...")
    long_summary <- long %>%
      group_by(genotype, time, geneID, trt, generation,
               Working.Symbol, ASE_sign) %>%
      summarise(
        mean_norm_counts = mean(norm_counts, na.rm = TRUE),
        se_norm_counts = sd(norm_counts, na.rm = TRUE) / sqrt(n()),
        mean_amreads = mean(Total_AMreads, na.rm = TRUE),
        se_amreads = sd(Total_AMreads, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
      )
    
    if (nrow(long_summary) == 0) {
      stop("Summary table is empty. Something went wrong downstream of filtering.")
    }
    
    ## Step 7: Order genes by ASE sign
    ordered_genes <- long_summary %>%
      distinct(Working.Symbol, ASE_sign) %>%
      group_by(Working.Symbol) %>%
      slice(1) %>%
      ungroup() %>%
      mutate(ASE_sign = factor(ASE_sign, levels = c("up", "down"))) %>%
      arrange(ASE_sign, Working.Symbol)
    
    long_summary$Working.Symbol <- factor(
      long_summary$Working.Symbol,
      levels = ordered_genes$Working.Symbol
    )
    long$Working.Symbol <- factor(
      long$Working.Symbol,
      levels = ordered_genes$Working.Symbol
    )
    
    ## Step 8: Define plot dimensions
    num_genes <- length(unique(long$Working.Symbol))
    height <- max(5, min(15, 0.5 * num_genes))
    width  <- 5
    
    ## Step 9: Expression plot by treatment
    message("Generating expression plot...")
    p1 <- ggplot() +
      geom_point(
        data = long,
        aes(x = trt, y = norm_counts, color = genotype,
            group = interaction(geneID, genotype)),
        alpha = 0.1
      ) +
      geom_line(
        data = long_summary,
        aes(x = trt, y = mean_norm_counts, color = genotype,
            group = interaction(geneID, genotype)),
        linewidth = 1
      ) +
      geom_point(
        data = long_summary,
        aes(x = trt, y = mean_norm_counts, color = genotype),
        size = 3
      ) +
      geom_errorbar(
        data = long_summary,
        aes(
          x = trt,
          ymin = mean_norm_counts - se_norm_counts,
          ymax = mean_norm_counts + se_norm_counts,
          color = genotype
        ),
        width = 0.1
      ) +
      facet_grid(time ~ Working.Symbol, scales = "free_y") +
      scale_color_brewer(palette = "Set2") +
      labs(x = "Treatment", y = "Normalized Counts ± SE") +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "top",
        strip.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 13, face = "bold"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        axis.line = element_line(color = "black")
      )
    
    ## Step 10: Norm counts vs Total AM reads
    message("Generating AM reads plot...")
    p2 <- ggplot() +
      geom_point(
        data = long,
        aes(x = norm_counts, y = Total_AMreads, color = genotype),
        alpha = 0.1
      ) +
      geom_point(
        data = long_summary,
        aes(x = mean_norm_counts, y = mean_amreads, color = genotype),
        size = 3
      ) +
      geom_errorbar(
        data = long_summary,
        aes(
          x = mean_norm_counts,
          ymin = mean_amreads - se_amreads,
          ymax = mean_amreads + se_amreads,
          color = genotype
        ),
        width = 0.1
      ) +
      geom_errorbarh(
        data = long_summary,
        aes(
          y = mean_amreads,
          xmin = mean_norm_counts - se_norm_counts,
          xmax = mean_norm_counts + se_norm_counts,
          color = genotype
        ),
        height = 0.1
      ) +
      facet_grid(time + trt ~ Working.Symbol, scales = "free") +
      scale_color_brewer(palette = "Set2") +
      labs(x = "Normalized Counts ± SE", y = "Total AM Reads ± SE") +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "top",
        strip.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 13, face = "bold"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        axis.line = element_line(color = "black")
      )
    
    # Ensure /plots folder exists
    dir.create("plots", showWarnings = FALSE)
    
    # Save plots
    message("Saving expression plot...")
    de_expr_path <- file.path("plots", paste0(title_prefix, "_expression_plot.png"))
    ggsave(de_expr_path, p1, width = width, height = height, dpi = 600)
    
    message("Saving AM reads plot...")
    de_amread_path <- file.path("plots", paste0(title_prefix, "_norm_vs_amreads.png"))
    ggsave(de_amread_path, p2, width = width, height = height, dpi = 600)
    
    # Return paths to app
    list(
      expr_plot    = de_expr_path,
      amreads_plot = de_amread_path
    )
    
  }, error = function(e) {
    message("Error in plot_expression_and_amreads: ", e$message)
    stop(e)
  })
}

##---------------------------------------------------------
## FUNCTION: plot_ase_expression
##---------------------------------------------------------

plot_ase_expression <- function(query_genes,
                                title_prefix,
                                trt_filter = c("M", "C"),
                                time_filter = c("T1", "T2")) {
  tryCatch({
    ## Step 1: Load gene annotation table
    message("Loading annotation table...")
    DE_ASE_table <- read_csv("data/DE_ASE_table.csv", progress = FALSE)
    
    selected_genes <- DE_ASE_table %>%
      filter(geneID %in% query_genes | Working.Symbol %in% query_genes)
    
    if (nrow(selected_genes) == 0) {
      stop("No genes found matching the query. Please check your gene IDs/symbols.")
    }
    
    ## Step 2: Load ASE-normalized counts
    message("Loading ASE counts...")
    norm_counts_ASE <- read_csv("data/norm_counts_ASE.csv", progress = FALSE) %>%
      select(-any_of("...1")) %>%
      filter(geneID %in% selected_genes$geneID) %>%
      right_join(selected_genes, by = "geneID")
    
    if (nrow(norm_counts_ASE) == 0) {
      stop("No rows in norm_counts_ASE after filtering for selected genes.")
    }
    
    ## Step 3: Load Total AM reads file
    message("Loading AM reads data...")
    Total_AMreads_df <- read_csv("data/Tota_AMreads_ineachsample.csv",
                                 progress = FALSE)
    
    ## Step 4: Reshape data and annotate metadata
    message("Processing data for ASE plotting...")
    long_ASE <- norm_counts_ASE %>%
      pivot_longer(
        cols = matches("_M_|_C_"),
        names_to = "Sample",
        values_to = "norm_counts"
      ) %>%
      mutate(label = gsub("E_|_O|_B", "", Sample)) %>%
      left_join(Total_AMreads_df, by = "label") %>%
      separate(
        Sample,
        into = c("delete", "genotype", "trt", "time", "rep", "allele"),
        sep = "_",
        remove = FALSE
      ) %>%
      mutate(
        generation = case_when(
          genotype == "Oh43xB73" ~ "F1",
          genotype %in% c("Oh43", "B73") ~ "F0",
          TRUE ~ NA_character_
        ),
        allele = case_when(
          allele == "B" ~ "B73",
          allele == "O" ~ "Oh43",
          TRUE ~ allele
        ),
        genotype_1 = case_when(
          genotype == "Oh43xB73" & allele == "B73"   ~ "F1_B73allele",
          genotype == "Oh43xB73" & allele == "Oh43"  ~ "F1_Oh43allele",
          genotype == "Oh43"                         ~ "Oh43",
          genotype == "B73"                          ~ "B73",
          TRUE                                       ~ NA_character_
        ),
        generation = factor(generation, levels = c("F0", "F1")),
        genotype_1 = factor(
          genotype_1,
          levels = c("Oh43", "B73", "F1_Oh43allele", "F1_B73allele")
        )
      ) %>%
      filter(!str_detect(genotype, "W22"))
    
    # Apply filters
    long_ASE <- long_ASE %>%
      filter(trt %in% trt_filter, time %in% time_filter)
    
    if (nrow(long_ASE) == 0) {
      stop("No ASE data left after filtering by treatment/time.")
    }
    
    ## Step 5: Summary statistics
    message("Calculating ASE summary statistics...")
    long_ASE_summary <- long_ASE %>%
      group_by(allele, time, geneID, trt, generation,
               Working.Symbol, ASE_sign) %>%
      summarise(
        mean_norm_counts = mean(norm_counts, na.rm = TRUE),
        se_norm_counts = sd(norm_counts, na.rm = TRUE) / sqrt(n()),
        .groups = "drop"
      )
    
    if (nrow(long_ASE_summary) == 0) {
      stop("ASE summary table is empty after grouping.")
    }
    
    ## Step 6: Order genes by ASE sign
    ordered_genes <- long_ASE_summary %>%
      distinct(Working.Symbol, ASE_sign) %>%
      group_by(Working.Symbol) %>%
      slice(1) %>%
      ungroup() %>%
      mutate(ASE_sign = factor(ASE_sign, levels = c("up", "down"))) %>%
      arrange(ASE_sign, Working.Symbol)
    
    long_ASE_summary$Working.Symbol <- factor(
      long_ASE_summary$Working.Symbol,
      levels = ordered_genes$Working.Symbol
    )
    long_ASE$Working.Symbol <- factor(
      long_ASE$Working.Symbol,
      levels = ordered_genes$Working.Symbol
    )
    
    ## Step 7: Define plot dimensions
    num_genes <- length(unique(long_ASE$Working.Symbol))
    height <- max(5, min(15, 0.5 * num_genes))
    width  <- 5
    
    ## Step 8: Generate ASE plot
    message("Generating ASE plot...")
    p <- ggplot() +
      geom_point(
        data = long_ASE,
        aes(x = generation, y = norm_counts,
            color = allele,
            group = interaction(geneID, allele)),
        alpha = 0.1
      ) +
      geom_line(
        data = long_ASE_summary,
        aes(x = generation, y = mean_norm_counts,
            color = allele,
            group = interaction(geneID, allele)),
        linewidth = 1
      ) +
      geom_point(
        data = long_ASE_summary,
        aes(x = generation, y = mean_norm_counts, color = allele),
        size = 3
      ) +
      geom_errorbar(
        data = long_ASE_summary,
        aes(
          x = generation,
          ymin = mean_norm_counts - se_norm_counts,
          ymax = mean_norm_counts + se_norm_counts,
          color = allele
        ),
        width = 0.1
      ) +
      facet_grid(time + trt ~ Working.Symbol, scales = "free_y") +
      scale_color_brewer(palette = "Set2") +
      labs(x = "Generation", y = "Normalized Counts ± SE") +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "top",
        strip.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 13, face = "bold"),
        panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        axis.line = element_line(color = "black")
      )
    
    # Ensure /plots folder exists
    dir.create("plots", showWarnings = FALSE)
    
    message("Saving ASE plot...")
    ase_path <- file.path("plots", paste0(title_prefix, "_ASE_expression_plot.png"))
    ggsave(ase_path, p, width = width, height = height, dpi = 600)
    
    return(TRUE)
    
  }, error = function(e) {
    message("Error in plot_ase_expression: ", e$message)
    stop(e)
  })
}