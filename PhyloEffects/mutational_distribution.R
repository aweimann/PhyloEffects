library(tidyverse)
library("scales")
library(ggrepel)
library(forcats)
library(fuzzyjoin)
library(stringr)
library(argparse)

#library(ggtree)
#library(ape)

theme_set(theme_classic(base_size = 18))


reverselog_trans <- function(base = exp(1)) {
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    trans_new(paste0("reverselog-", format(base)), trans, inv, 
              log_breaks(base = base), 
              domain = c(1e-100, Inf))
}


scientific_10 <- function(x) {
      ifelse(x == 0, "0",
      parse(text=gsub("e[+]", " %*% 10^", scales::scientific_format()(x))))
}

#manhattan plot function
manhattan <- function(g, by){    
g  %>% ggplot(aes(position, pval, label = ifelse(!is.na(gene), gene, locus_tag))) +#, color = !!sym(by))) +
    geom_hline(aes(yintercept = (g %>% filter(padj < 0.05) %>% ungroup() %>% summarize(max(pval)))[[1]])) +
    geom_point(pch = 21) +
    geom_text_repel(size = 6) +
    scale_x_continuous(label = scientific_10) +
    scale_y_continuous(trans = reverselog_trans(10),
        breaks = trans_breaks("log10", function(x) 10^x),
        labels = trans_format("log10", math_format(.x)))+
    # scale_colour_manual(values = c("grey", "black"), guide = 'none') +
    scale_size_manual(guide = 'none',  values = c(sig = 3, `non-sig` = 1), 1 ) +
    ylab(expression(log[10]~group("(", p-value, ")")))+
    xlab("Genomic position")
#https://slowkow.com/notes/ggplot2-qqplot/
}

#qqplot function
gg_qqplot <- function(ps, ci = 0.95) {
  n  <- length(ps)
  df <- data.frame(
    observed = -log10(sort(ps)),
    expected = -log10(ppoints(n)),
    clower   = -log10(qbeta(p = (1 - ci) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper   = -log10(qbeta(p = (1 + ci) / 2, shape1 = 1:n, shape2 = n:1))
  )
  log10Pe <- expression(paste("Expected -log"[10], "p-value"))
  log10Po <- expression(paste("Observed -log"[10], "p-value"))
  ggplot(df) +
    geom_point(aes(expected, observed), shape = 1, size = 1) +
    geom_abline(intercept = 0, slope = 1, alpha = 0.5) +
    geom_line(aes(expected, cupper), linetype = 2) +
    geom_line(aes(expected, clower), linetype = 2) +
    geom_hline(yintercept = 2.657577, size = 0.1) +
    theme_light(base_size = 26)+
    xlab(log10Pe) +
    ylab(log10Po)
}


# v <-  st_df %>% select(st, mutations) %>% mutate(mutations = map(mutations, ~ read_tsv(., col_types = cols(.default = col_character(), pos = col_integer())))) %>% unnest(cols = c(mutations))  %>% rename(ref = upstream_allele, alt = downstream_allele)

#when looking at intergenic regions only consider upstream gene variants
# wd <- "/Users/aweimann/Library/CloudStorage/OneDrive-UniversityofCambridge/cambridge_postdoc/current_results/burden_test_spectra_dataset/"
# phyloeffects_in <- "phyloeffects_out/S_aureus_CC398/"
# burden_test_out <- "burden_test_out/S_aureus_CC398/"
# dataset <- "S_aureus_CC398"
# setwd(wd)
# args <- c(dataset, phyloeffects_in, burden_test_out, wd)
options(show.error.locations = TRUE)


# make options for command line arguments
# include a flag to specify whether to remove recombination regions or not
parser <- ArgumentParser()
parser$add_argument("-d", "--dataset", type = "character", required = TRUE,
                    help = "Comma separated list of datasets to include in the analysis")
parser$add_argument("-w", "--working_dir", type = "character", required = TRUE,
                    help = "Working directory where the input files are located and output files will be saved")
parser$add_argument("-p", "--phyloeffects_in", type = "character", required = TRUE,
                    help = "Directory where the phyloeffects output files are located")
parser$add_argument("-b", "--burden_test_out", type = "character", required = TRUE,
                    help = "Directory where the burden test output files will be saved")
parser$add_argument("-a", "--gene_annotation", type = "character", required = TRUE,
                    help = "Gene annotation file in tsv format with columns: locus_tag, gene_name, feature, product, start, end")
parser$add_argument("-r", "--no_recombination_removal", action = "store_true", default = FALSE,
                    help = "Flag to specify whether to remove recombination regions")
parser$add_argument("-i", "--include_indels", action = "store_true", default = FALSE,
                    help = "Flag to specify whether to include indel events")
parser$add_argument("-t", "--test_intergenic", action = "store_true", default = FALSE,
                    help = "Flag to specify whether to test intergenic regions")

opt <- parser$parse_args()
dataset <- opt$dataset
# comma separated list of STs to include in the analysis
sts <- dataset
working_dir <- opt$working_dir
print(working_dir)
setwd(working_dir)
phyloeffects_in <- opt$phyloeffects_in
burden_test_out <- opt$burden_test_out
gene_annotation <- opt$gene_annotation
dataset <-  dataset %>% str_split(",") %>% unlist() %>% str_trim()
# flag to specify whether to remove recombination regions or not
no_recombination_removal <- opt$no_recombination_removal

annot <- read_tsv(gene_annotation)

# remove tRNA, rRNA and transposases
# annot <- annot %>% filter(!(feature %in% c("rRNA", "tRNA")) ) %>% 
# filter(!str_detect(product, "transposase"))

sts = c(dataset)
st_df <- tibble(st = sts)


st_df <-
    mutate(st_df, mutations = str_c(phyloeffects_in, st, "/", st, ".variant_effect_predictions.txt")) %>%
    mutate(st_df, recombination = str_c(phyloeffects_in, st, "/",st, ".recombination_prediction.txt")) %>%
    mutate(st_df, indel_events = str_c(phyloeffects_in, st, "/",st, ".indel_events.txt")) %>%
    mutate(st_df, indel_variant_effects = str_c(phyloeffects_in, st, "/",st, ".indel.variant_effect_predictions.txt")) %>%
    mutate(st_df, recombination_pos = str_c(phyloeffects_in, st, "/", st, ".recombination_pos.txt"))

v <-  st_df %>% select(st, mutations) %>% mutate(mutations = map(mutations, ~ read_tsv(., col_types = cols(.default = col_character(), pos = col_integer())))) %>% unnest(cols = c(mutations))  %>% rename(ref = upstream_aa, alt = downstream_aa)
v <- v %>% separate_rows(samples, sep = ",")


if (opt$include_indels) {
  # Read in indels from events and variant effect predictions
   events <- st_df %>% select(st, indel_events) %>%
     mutate(indel_events = map(indel_events, ~ read_tsv(., col_types = cols(.default = col_character())))) %>%
     unnest(cols = c(indel_events))
  
   effect_indels <- st_df %>% select(st, indel_variant_effects) %>%
     mutate(indel_variant_effects = map(indel_variant_effects, ~ read_tsv(., col_types = cols(.default = col_character(), pos = col_integer())))) %>%
     unnest(cols = c(indel_variant_effects)) %>%
     select(-samples)

    # # Parse variant_id from events (format: NC_009648:pos:ref:alt)
   events <- events %>%
     separate(variant_id, c("ref_genome", "pos", "ref", "alt"), sep = ":") %>%
     mutate(pos = as.integer(pos), ref = toupper(ref), alt = toupper(alt)) %>%
     rename(samples = node)
  
  # # Join events with effect predictions on st, pos, and alleles
   indels <- events %>%
     inner_join(effect_indels, by = c("st", "pos", "ref" = "upstream_allele", "alt" = "downstream_allele"), relationship = "many-to-many") %>%
     select(st, samples, pos, ref, alt, parent_node, node_state, parent_node_state, impact, mutation_type, locus_tag) %>%
     inner_join(annot) %>%
     group_by(st, pos, ref, alt) %>% add_tally() %>% filter(n == 1) %>% select(-n) %>% ungroup()
 
 v <- bind_rows(v, indels)
 }

# filter variants in intergenic regions
v <- v %>% filter(impact != "MODIFIER")

#remove recombination sites
#read in recombination position and interval 
if (!no_recombination_removal) {
  recombination <- st_df %>% select(st, recombination) %>%mutate(recombination = map(recombination, ~ read_tsv(., col_types = cols(.default = col_character(), start = col_integer(), stop = col_integer()))))  %>% unnest(cols = c(recombination)) %>% rename(end = stop)
  gubbins_embl <- st_df %>% select(st, recombination_pos) %>%mutate(recombination_pos = map(recombination_pos, ~ read_tsv(., col_types = cols(.default = col_character(), pos = col_integer()))))  %>% unnest(cols = c(recombination_pos)) 
  recomb_pos <- gubbins_embl %>% mutate(start = pos, end = pos + 1)  %>% genome_join(recombination,  by = c("st", "start", "end")) %>% filter(node.x == node.y) %>% select(-node.y, - start.y, -st.y) %>% rename(st = st.x) %>% select(pos, "st")
}

#remove recombination from multi codon substitutions and standard substitutions separately as MCS positions might not agree with Gubbins positions
if (!no_recombination_removal) {
v_mcs <- v  %>% filter(!is.na(multi_codon_substitution)) %>% mutate(start = pos, end =  pos + 2) %>%  genome_anti_join(recomb_pos %>% mutate(start = pos, end = pos + 1), by = c("st", "start", "end"))  %>% select(-start, -end)
v <- v %>% filter(is.na(multi_codon_substitution)) %>% anti_join(recomb_pos)
v <- bind_rows(v, v_mcs)
}


# Klebsiella
# v <- v %>% filter(samples != "NC_009648.1")
# v <- v %>% filter(st != "st25" & !(samples %in% c("NC_009648.1", "Node_37", "Node_474", "ERR775516", "Node_476", "Node_461", "Node_475")))
# v <- v %>% filter(st != "st14" & !(samples %in% c("Node_2421")))
# v <- v %>% filter(st != "st258" | !(samples %in% c("NC_009648.1", "Node_496", "Node_498", 
#                                                          "Node_497", "Node_495", "ERR1334502", "ERR1217461")))


#intergenic mutations
# intergenic_mutations <- v %>% group_by(st, node, pos, ref, alt) %>% filter(mutation_type == "upstream_gene_variant" | mutation_type == "downstream_gene_variant") %>% mutate(PAO1 = ifelse(length(PAO1) == 2, str_c(PAO1[1], PAO1[2], sep = ","), PAO1), gene_name = ifelse(length(gene_name) == 2,str_c(gene_name[1], gene_name[2], sep = ","), gene_name), n = length(PAO1)) %>% group_by(st, node, pos, ref, alt, PAO1, gene_name) %>% count() %>% group_by(PAO1, gene_name) %>% count() %>% arrange(-n)
# intergenic_regions <- read_tsv("intergenic_regions.txt")
# intergenic_mutations <- intergenic_mutations %>% inner_join(intergenic_regions) %>% mutate(gene_length = end - start)
# total_length_intergenic <- intergenic_regions %>% mutate(gene_length = end - start) %>% ungroup() %>% summarize(sum(gene_length))
# total_mutations <- intergenic_mutations %>% ungroup() %>% summarize(sum(n))
# intergenic_mutations <- intergenic_mutations   %>%  mutate(pval = poisson.test(n, r=total_mutations[[1]]*(gene_length/(total_length_intergenic[[1]])), alternative = "greater")['p.value'][[1]])
# intergenic_mutations$padj <- p.adjust(intergenic_mutations$pval, method = "BH")
# v <- v %>%  filter(mutation_type != "upstream_gene_variant" & mutation_type != "downstream_gene_variant")
v <- v %>%  filter(mutation_type != "upstream_gene_variant", mutation_type != "downstream_gene_variant")
# v <- v %>%  filter(mutation_type == "downstream_gene_variant")


#stratify by synonymous vs. non-synonymous variants 
mod <-  v  %>% group_by(locus_tag) %>% filter(impact == "LOW") %>% count() %>% arrange(-n) 
mod <- annot %>% select(locus_tag) %>% left_join(mod) %>% mutate(n = ifelse(is.na(n), 0, n))
high <-  v  %>% group_by(locus_tag) %>%   filter(impact != "LOW") %>% count() %>% arrange(-n) 

#count global number of mutations
comb <- mod %>% full_join(high, by = c("locus_tag")) %>% mutate(n.x = ifelse(is.na(n.x), 0, n.x), n.y = ifelse(is.na(n.y), 0, n.y))
comb <- comb %>% mutate(dn_ds = n.y/n.x)  %>% mutate(dn_ds = ifelse(is.infinite(dn_ds), n.y, dn_ds ))  %>% arrange(-dn_ds) 
#stratify between point mutations and structural variation
mutation_type <- v %>% mutate(mutation_type = ifelse(str_length(ref) != str_length(alt), "indel", "point_mutation")) %>%   group_by(locus_tag, mutation_type)  %>%  count()  %>% spread(mutation_type, n, fill = 0)
#stratify by impact 
impact <- v %>%   group_by(locus_tag, impact)  %>% count() %>% spread(impact, n, fill = 0)
#non-synonymous variants per ST
# per_patient <- v  %>%  filter(impact != "LOW") %>%  group_by(locus_tag, gene_name, st) %>% count() %>% spread(st, n, fill = 0)
#count number of STs mutations in a particular gene are found 
#per_st <- v  %>%  filter(impact != "LOW") %>%  group_by(locus_tag, st) %>% count()  %>%  group_by(PAO1, gene_name) %>% count() %>% rename(no_sts = n)
comb <- comb %>% left_join(impact) %>% left_join(mutation_type)
uq_variants <-  v %>% group_by(locus_tag, pos, upstream_allele, downstream_allele) %>% filter(impact != "LOW") %>% 
    summarize(unique_variants = 1) %>% group_by(locus_tag) %>% 
    summarize(unique_variants = sum(unique_variants))
# ref_genome_coverage <- read_tsv("reference_genes_coverage.txt")
# total length of genes in the genome 
total_length <- (mutate(annot, gene_length = end - start)  %>% ungroup() %>%  summarize(total_length = sum(gene_length)))$total_length
#subtract numebr of overlapping base pairs (0) 
total_length <- (mutate(annot, gene_length = end - start)  %>% 
                 ungroup() %>% 
                 summarize(total_length = sum(gene_length) - 0 ))$total_length
comb <- comb %>% left_join(uq_variants) %>%
    inner_join(annot) %>% 
    mutate(gene_length = end - start, position = start) 
if("HIGH" %in% colnames(comb)){
  comb <- comb %>% select(locus_tag, n.x, n.y, dn_ds, gene_length,  position, HIGH:point_mutation, unique_variants) 
}else if("MODERATE" %in% colnames(comb)){
  comb <- comb %>% select(locus_tag, n.x, n.y, dn_ds, gene_length,  position, MODERATE:point_mutation, unique_variants) 
}else{
  comb <- comb %>% select(locus_tag, n.x, n.y, dn_ds, gene_length,  position, unique_variants) 
}
comb <- comb %>% mutate(d.x_mod = n.y * 1000/gene_length) %>% arrange(-d.x_mod) 
no_mutations <- comb %>% ungroup() %>% summarize(sum(n.y)) %>% as_vector()

# comb <- comb   %>%  
#     rowwise() %>% 
#     mutate(r = no_mutations*(gene_length /(total_length)), pval = poisson.test(n.y,r=r, alternative = "greater")['p.value'][[1]])

poisson_midp <- function(x, lambda) {
  # P(X > x)
  p_greater <- ppois(x, lambda, lower.tail = FALSE)
  # P(X = x)
  p_equal <- dpois(x, lambda)
  
  # Mid-p formula
  mid_p <- p_greater + 0.5 * p_equal
  return(mid_p)
}

comb <- comb  %>% mutate(r = no_mutations*(gene_length /(total_length)), 
                         pval = poisson_midp(n.y, lambda = r))
comb$padj <- p.adjust(comb$pval, method = "BH")



comb <- comb %>% mutate(is_sig = ifelse(padj < 0.05, "sig", "non-sig"))
# comb <- comb %>% inner_join(per_st)
# ggplot(comb, aes(no_sts, padj, color = gene_length)) + geom_point(size = 1) + scale_y_log10() + scale_color_continuous(trans = 'log2') + geom_hline(yintercept = 0.05)
# ggsave("padj_vs_no_sts.png")


# comb %>% filter(padj < 0.05) %>% select(locus_tag)  %>% inner_join(v) %>% filter(impact != "LOW")%>% write_tsv("variants_in_burden_hits.txt")

comb <- comb  %>% arrange(padj) %>% inner_join(annot) %>% write_tsv(str_c(burden_test_out, "/poisson_test.txt"))
manhattan(comb)
ggsave(str_c(burden_test_out, "/manhattan_poisson_test.pdf"), width = 18)
gg_qqplot(comb$pval)

ggsave(str_c(burden_test_out, "/qqplot_poisson_test.pdf"))

burden_variants <- comb %>% filter(padj < 0.05) %>% select(locus_tag)  %>% inner_join(v) %>% filter(impact != "LOW")
v %>% write_tsv(str_c(burden_test_out, "/variants_in_burden_hits.txt"))

# test internal vs terminal nodes
is_transmitted <- v %>% filter(mutation_type != "LOW") %>% mutate(is_internal = str_detect(samples, "Node")) %>% 
  group_by(locus_tag, is_internal) %>% 
  count() %>% 
  pivot_wider(names_from = is_internal, values_from = n, values_fill = 0) 
is_transmitted <- rename(is_transmitted, transmitted = `TRUE`, untransmitted = `FALSE`) 
print("test")
  
total <- is_transmitted %>% ungroup () %>% summarise(untransmitted = sum(untransmitted), transmitted = sum(transmitted))
is_transmitted$untransmitted_total <- total$untransmitted
is_transmitted$transmitted_total <- total$transmitted
is_transmitted %>% inner_join(comb %>% filter(padj < 0.05))
f_test <- is_transmitted %>% inner_join(comb %>% filter(padj < 0.05))  %>% 
  mutate(ftest = list(fisher.test(matrix(c(transmitted, untransmitted, transmitted_total, untransmitted_total), nrow = 2))))
f_test <- f_test %>% mutate(transmission_pval = ftest[[1]]$p.value, transmission_odds_ratio = ftest[[1]]$estimate) 
f_test$transmission_padj <- p.adjust(f_test$transmission_pval, method = "BH")
comb <- comb %>% left_join(f_test %>% select(-ftest))
comb <- comb  %>% arrange(padj) %>% write_tsv(str_c(burden_test_out, "/poisson_test_w_transmission.txt"))


# test animal vs human
# is_human <- v %>% filter(mutation_type != "LOW") %>% 
#   filter(!str_detect(samples, "NODE")) %>% 
#   mutate(is_internal = str_detect(samples, "Human")) %>% 
#   group_by(locus_tag, is_internal) %>% 
#   count() %>% 
#   pivot_wider(names_from = is_internal, values_from = n, values_fill = 0) 
# is_human <- rename(is_human, human = `TRUE`, animal = `FALSE`) 
#   
# total <- is_human %>% ungroup () %>% summarise(animal = sum(animal), human = sum(human))
# is_human$animal_total <- total$animal
# is_human$human_total <- total$human
# is_human %>% inner_join(comb %>%  filter(padj < 0.05))
# f_test_human <- is_human %>% inner_join(comb %>% select(locus_tag, padj) %>% filter(padj < 0.05))  %>% 
#   mutate(ftest = list(fisher.test(matrix(c(human, animal, human_total, animal_total), nrow = 2))))
# f_test_human <- f_test_human %>% mutate(host_pval = ftest[[1]]$p.value, host_odds_ratio = ftest[[1]]$estimate) 
# f_test_human$host_padj <- p.adjust(f_test_human$host_pval, method = "BH")
# comb <- comb %>% left_join(f_test_human %>% select(-ftest))
# comb <- comb  %>% arrange(padj) %>% write_tsv(str_c(burden_test_out, "/", dataset, "_poisson_test_w_animal_vs_human.txt"))

