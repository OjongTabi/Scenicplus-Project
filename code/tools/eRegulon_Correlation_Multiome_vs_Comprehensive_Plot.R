#!/usr/bin/env Rscript

library(ggplot2)
library(readr)
library(ggrepel)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)  # Add this at the top


# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Ensure two arguments are provided
if (length(args) < 2) {
  stop("Incorrect inputs")
}

file <- args[1]
label <- args[2]

# Function to create scatter plots
create_scatter_plot <- function(file, label) {
  # Load data
  data <- read_csv(file)
  
  data <- data %>%
  mutate(Annotation = case_when(
    Multiome_Corr != 0 & Comprehensive_Corr == 0 ~ "Multiome_Corr",
    Multiome_Corr == 0 & Comprehensive_Corr != 0 ~ "Comprehensive_Corr",
    Multiome_Corr != 0 & Comprehensive_Corr != 0 ~ "Common",
    TRUE ~ "None"
  ))


  # Create scatter plot without annotations
  plot_no_annotations <- ggplot(data, aes(x = Multiome_Corr, y = Comprehensive_Corr, color = Annotation)) +
    geom_point(size = 1, alpha = 0.8) +  # Use color aesthetic for points
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    labs(title = paste(label, "(No Annotations)"), x = "Multiome Correlation", y = "Comprehensive Correlation") +
    theme_minimal() +
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.1)) +  # Set x-axis range and steps
    scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.1)) +  # Set y-axis range and steps
    coord_fixed(ratio = 1) + # Ensure same scale for both axes
    scale_color_manual(values = c("Multiome_Corr" = "green", 
                                  "Comprehensive_Corr" = "brown", 
                                  "Common" = "blue")) +  # Assign custom colors
    theme(legend.position = "right")  # Adjust legend position
  
  # Create scatter plot with annotations using ggrepel
  plot_with_annotations <- ggplot(data, aes(x = Multiome_Corr, y = Comprehensive_Corr, label = eRegulon, color = Annotation)) +
    geom_point(size = 1, alpha = 0.8) +  # Use color aesthetic for points
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    geom_text_repel(max.overlaps = nrow(data),  # Ensure all data points are annotated
                    size = 1.75, 
                    fontface = "bold", 
                    box.padding = 0.1,
                    point.padding = 0.1,
                    force = 1,  # Adjust force to spread labels further apart
                    segment.size = 0.5,  # Set thinner lines between text and points
                    color = "black") +  # Bold text and prevent overcrowding
    labs(title = paste(label, "(With Annotations)"), x = "Multiome Correlation", y = "Comprehensive Correlation") +
    theme_minimal() +
    scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.1)) +  # Set x-axis range and steps
    scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.1)) +  # Set y-axis range and steps
    coord_fixed(ratio = 1) + # Ensure same scale for both axes
    scale_color_manual(values = c("Multiome_Corr" = "green", 
                                  "Comprehensive_Corr" = "brown", 
                                  "Common" = "blue")) +  # Assign custom colors
    theme(legend.position = "right")  # Adjust legend position
  
  # Return the plots
  return(list(plot_no_annotations, plot_with_annotations))
}

# Generate a single plot
plots <- create_scatter_plot(file, label)

combined_plot <- plots[[1]] + plots[[2]]  

# Save the combined plot
ggsave("corr_plot.png", combined_plot, width = 14, height = 20)

