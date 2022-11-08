#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(tidyverse))

count_observed_haplotypes <- function(data, gene, allelic_group, allele) {
  
  target_haplotype <- paste(gene, allelic_group, sep="*")
  
  counted_data <- data %>% 
    separate(
      `observed_haplotypes`,
      into=c("Gene","Allelic_Group","Allele"),
      sep=':|\\*',
      fill='right',
      remove = F
      ) %>%
    filter(Gene == gene & Allelic_Group == allelic_group)

  if ( !is.na(allele) ) {
    counted_data <- filter(counted_data, Allele == allele)
    target_haplotype <- paste(target_haplotype, allele, sep=":")
  }

  result <- counted_data %>%
    select(observed_haplotypes) %>%
    mutate("target_haplotype" = target_haplotype) %>%
    group_by(target_haplotype) %>%
    summarise(
      n=n(),
      .groups='keep'
    )
  
  ## If haplotype is not present, return a count of 0
  if ( nrow(result) == 0 ) {
    result <- tibble("target_haplotype" = target_haplotype, "n"=0)
  }

  return(result)
}

##### MAIN ######
parser <- OptionParser(usage = "usage: %prog (-i input.txt) (-t targets.txt) [-o output.txt]", description="A helper tool to count the number of HLA haplotypes matching the target alleles from arlequin format and print them out in Rock format.")
parser <- add_option(parser, c("-i", "--input"),
  type = "character",
  action = "store", dest = "input_fn",
  help = "The input arlequin file."
)
parser <- add_option(parser, c("-t", "--targets"),
  type = "character",
  action = "store", dest = "targets_fn",
  help = "A file containing the desired target haplotypes, one per line."
)
parser <- add_option(parser, c("-o", "--output"),
  type = "character",
  action = "store", dest = "output_fn",
  help = "The desired output Rock-format file."
)

args <- parse_args(parser)

## validate inputs
if ( !file.exists(args$input_fn) ) write("Input file does not exist!", file = stderr())
if ( !file.exists(args$targets_fn) ) write("Input file does not exist!", file = stderr())

## If no output name is given, print to stdout.
if ( is.null(args$output_fn) ) {
  output_file <- stdout()
} else {
  output_file <- args$output_fn
}

targets <- read_tsv(args$targets_fn, col_names=c("Targets"), col_types='c', show_col_types = FALSE) %>%
  separate(`Targets`, into=c("Gene","Allelic_Group","Allele"), sep=':|\\*', fill='right')

data <- read_tsv(args$input_fn, col_names=c("Ind", "Pop"), show_col_types = FALSE) %>%
  select(-Ind, -Pop) %>%
  pivot_longer(cols=everything(), names_to='filler', values_to="observed_haplotypes") %>%
  select(-filler)

target_haplotype_counts = targets %>% 
  purrr::pmap_dfr(~count_observed_haplotypes(data, ..1, ..2, ..3)) %>%
  pivot_wider(id_cols=everything(), names_from=target_haplotype, values_from=n)

write_tsv(target_haplotype_counts, file = output_file)
