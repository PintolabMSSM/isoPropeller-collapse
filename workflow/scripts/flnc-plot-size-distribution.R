#!/usr/bin/env Rscript

# 13.05.2025 14:36:49 EDT

#############
# ARGUMENTS #
#############

# Argparse library
library(argparse)

# Initialize an ArgumentParser object
parser <- ArgumentParser(description = "Plot FLNC size distributions")

# Defining named arguments
parser$add_argument("-b", "--bin_size", 
                    type     = "integer",
                    default  = 200,
                    help = "Bin size in nucleotides. Default: 200")
parser$add_argument("-m", "--method", 
                    type     = "character",
                    default  = "median",
                    help = "The metric for summarizing the read counts in each bin. One of sum, average, or median. Default: sum")
parser$add_argument("-o", "--output", 
                    type     = "character",
                    default  = "flnc-read-length-distribution.svg",
                    help     = "Output file name")
parser$add_argument("-t", "--plot_type", 
                    type     = "character",
                    default  = "density",
                    help     = "Type of plot: density, boxplot, or violin. Default: density")
parser$add_argument("-g", "--mapping_file",
                    type    = "character",
                    default = NULL,
                    help    = "Optional CSV/TSV file mapping Files to Groups. Columns: File,Group")


# Defining positional arguments
parser$add_argument("files", 
                    type     = "character", 
                    nargs    = "*",
                    help     = "Extra positional arguments")

# Parsing the arguments
args <- parser$parse_args()

# Check if positional arguments are missing (fallback validation)
if (length(args$files) == 0) {
   cat("\nError: At least one positional argument is required.\n\n")
   parser$print_help()
   quit(status = 1)
}

#############
# LIBRARIES #
#############

library(ggplot2)
library(patchwork)
library(dplyr)
library(data.table)

#############
# FUNCTIONS #
#############

binned_histogram_data <- function(file, bin_size, method = "sum") {

   # Read and summarize data
   data           = fread(file, select = "insertlen")
   data           = as.data.frame(table(data$insertlen))
   colnames(data) = c("Number","Frequency")

   # Ensure data is numeric
   data$Number <- as.numeric(as.character(data$Number))

   # Create bins using floor and direct calculation
   data$Bin <- floor(data$Number / bin_size) * bin_size

   # Aggregate the data by bin
   aggregated_data <- data %>%
     group_by(Bin) %>%
     summarize(
       Frequency = case_when(
         method == "sum" ~ sum(Frequency, na.rm = TRUE),
         method == "average" ~ mean(Frequency, na.rm = TRUE),
         method == "median" ~ median(Frequency, na.rm = TRUE),
         TRUE ~ sum(Frequency, na.rm = TRUE) # Default to sum
       )
     ) %>%
   arrange(Bin)

   # Normalize the data
   aggregated_data$Frequency = aggregated_data$Frequency / sum(aggregated_data$Frequency)

   # Add file name and return
   aggregated_data$File = basename(file)
   return(aggregated_data)
}

sliding_window_binned_histogram_data <- function(file, bin_size, step_size, method = "sum") {

    # Read and summarize data
    data           = fread(file, select = "insertlen")
    data           = as.data.frame(table(data$insertlen))
    colnames(data) = c("Number","Frequency")

    # Ensure data is numeric
    data$Number <- as.numeric(as.character(data$Number))

    # Initialize an empty list for sliding window results
    sliding_window_data <- list()

    # Define the range of the window
    min_num <- min(data$Number)
    max_num <- max(data$Number)

    # Apply sliding window method
    for (start in seq(min_num, max_num - bin_size, by = step_size)) {
        end <- start + bin_size

        # Filter data within the window
        window_data <- data[data$Number >= start & data$Number < end, ]

        # Check if the window is empty
        if (nrow(window_data) == 0) {
          aggregated_value <- 0  # Set frequency to 0 for empty window
        } else {
          # Calculate the aggregated value
          aggregated_value <- case_when(
            method == "sum" ~ sum(window_data$Frequency, na.rm = TRUE),
            method == "average" ~ mean(window_data$Frequency, na.rm = TRUE),
            method == "median" ~ median(window_data$Frequency, na.rm = TRUE),
            TRUE ~ sum(window_data$Frequency, na.rm = TRUE)
          )
        }

        # Store the result
        sliding_window_data[[length(sliding_window_data) + 1]] <- data.frame(
          Bin = (start + end) / 2,  # Midpoint of the bin
          Frequency = aggregated_value
        )
    }

    # Combine the sliding window data
    aggregated_data <- do.call(rbind, sliding_window_data)

    # Normalize the data
    aggregated_data$Frequency = aggregated_data$Frequency / sum(aggregated_data$Frequency)

    # Add file name and return
    aggregated_data$File = basename(file)
    return(aggregated_data)
}


boxplot_data <- function(file) {
   # Read only the insertlen column
   data <- fread(file, select = "insertlen")

   # Summarize as a table
   tab           <- as.data.frame(table(data$insertlen))
   colnames(tab) <- c("InsertLen", "Frequency")
   tab$InsertLen <- as.numeric(as.character(tab$InsertLen))
   tab$Frequency <- as.numeric(tab$Frequency)
   tab$File      <- basename(file)

   return(tab)
}


map_files_to_groups <- function(plot_data, mapping, order_by_group = TRUE) {  
   # Merge mapping into plot_data by File name
   merged_data <- merge(plot_data, mapping, by.x = "File", by.y = "File", all.x = TRUE)
   
   # Handle unmapped Files
   merged_data$Group[is.na(merged_data$Group)] <- "Unmapped"
   merged_data$Group <- factor(merged_data$Group)
   
   # Order File factor by group then File, if requested
   if (order_by_group) {
      file_order       <- mapping[order(Group, File), File]

      # Include unmapped files at the end
      extras           <- setdiff(unique(plot_data$File), file_order)
      file_order       <- c(file_order, sort(extras))
      merged_data$File <- factor(merged_data$File, levels = file_order)
   }
   
   return(merged_data)
}


########
# MAIN #
########

# Gather data
plot_data = data.frame()

# Load mapping
if (!is.null(args$mapping_file)) {
   mapping <- fread(args$mapping_file)
   
   if (!all(c("File", "Group") %in% colnames(mapping))) {
      stop("Mapping must contain 'File' and 'Group' columns.")
   }
}

# Make density plots
if (args$plot_type %in% c("density")) {
   for ( file in args$files ){
      cat(paste0("Processing: ", file,"\n"))
      dset = sliding_window_binned_histogram_data(file, bin_size=args$bin_size, step_size=10, method=args$method)
      plot_data = rbind(plot_data, dset)
   }

   # Plot the binned histogram
   p <- ggplot(plot_data, aes(x = Bin, y = Frequency, color=File)) +
      geom_line(linewidth=1) +
      labs(
         title = paste0("FLNC histogram (", args$method, " binning)"),
         x = "Number (Base Pairs)",
         y = "Frequency"
      ) +
      xlim(0,7500) +
      theme_minimal()
   
   ggsave(p, file=args$output, width=20, height=8)
}

# Make boxplots and violin plots
if (args$plot_type %in% c("boxplot", "violin")) {
   for (file in args$files) {
      cat(paste0("Processing for box/violin: ", file, "\n"))
      dset <- boxplot_data(file)  # summarized, unexpanded
      plot_data <- rbind(plot_data, dset)
   }

   # Add optional mapping to groups
   if (!is.null(args$mapping_file)) {
      plot_data <- map_files_to_groups(plot_data, mapping)
   }

   # Check if 'Group' column exists and use it for fill
   use_group <- "Group" %in% colnames(plot_data)
   fill_aes  <- if (use_group) aes(fill = Group) else aes(fill = File)

   # Make boxplot using weights rather than simply parsing all values
   p <- ggplot(plot_data, aes(x = File, y = InsertLen, weight = Frequency)) + fill_aes 

   # Select the plot type
   if (args$plot_type == "boxplot") {
      p <- p + geom_boxplot(outlier.shape = NA)
   } else {
      p <- p + 
           geom_violin(trim = FALSE, scale = "count") +
           geom_boxplot(width = 0.1, outlier.shape = NA)
   }

   # Finalize plot
   p <- p +
      ylim(0,5000) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
   # Save to file
   ggsave(p, file=args$output, width=20, height=8)
}


# Clean up the SVG format for inkscape by removing 'textlength' and 'lengthAdjust' attributes\n",
if (grepl("\\.svg$", args$output, ignore.case = TRUE)){
   svg_content <- readLines(args$output)
   svg_content <- gsub(' textLength=["\'][^"\']+["\']', '', svg_content)
   svg_content <- gsub(' lengthAdjust=["\'][^"\']+["\']', '', svg_content)
   writeLines(svg_content, args$output)
}

