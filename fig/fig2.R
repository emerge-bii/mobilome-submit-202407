source(here::here('setup.R'))

#################
# data wrangling
##################

### dirctories
v2_mags_directory <- here::here("som-data", "fig-data", "emerge_mags_v2")
recombinase_directory <- here::here("som-data", "fig-data", "20230123_mge_recombinase", "recombinase")
v3_contig_tracking_directory <- here::here("som-data", "fig-data", "contig_tracking_v3")
recombinase_clustering_directory <- here::here("som-data", "fig-data", "recombinase_clustering_v1")

### other settings
contig_length_cutoff <- 3000
colour_brewer <- setNames(append(as.list(RColorBrewer::brewer.pal(12, "Paired")), c("#737373", "#FFFFFF", "#000000")), c("blue", "darkblue", "green", "darkgreen", "red", "darkred", "orange", "darkorange", "purple", "darkpurple", "yellow", "brown", "grey", "white", "black"))
habitat_levels <- c("Palsa", "Bog", "Fen")
colour_habitat <- c("#703C1B", "#058000", "#0001FF")
depth_levels <- c("0-9", "10-19", "20-29", "30-39")
depth_labels <- c("0", "10", "20", "30")
depth_fills   <- RColorBrewer::brewer.pal(5, "BuPu")[-1]
depth_colours <- depth_fills
year_levels <- c("2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017")
year_labels <- c("10", "11", "12", "13", "14", "15", "16", "17")
taxa_levels <- c("none", "phylum", "class", "order", "family", "genus", "species")
taxa_colours <- c("#777777", RColorBrewer::brewer.pal(6, "Dark2"))
phyla_levels <- c(
  "Acidobacteriota",
  "Actinomycetota",
  "Verrucomicrobiota",
  "Pseudomonadota",
  "Desulfobacterota",
  "Dormibacterota",
  "Bacteroidota",
  "Chloroflexota",
  "Nitrospirota",
  "Myxococcota",
  "Halobacteriota",
  "Bacillota_B",
    "Patescibacteria",
    "Planctomycetota",
    "Eremiobacterota",
    "Gemmatimonadota",
  "Other"
  )
phyla_colours <- ggsci::pal_d3(palette = 'category20', alpha=1)(20)[-4] %>% #remove red, not good w/ green
  head(16) %>% c("grey50")
phyla_colours[7] <- 'grey20' # change grey to black; was postion 8 but becomes 7 after removing red (4)
phyla_colours_lines <- phyla_colours
selected_phyla <- c(
  "Acidobacteriota",
  "Actinomycetota",
  "Verrucomicrobiota",
    "Pseudomonadota",
  "Desulfobacterota",
  "Dormibacterota",
  "Bacteroidota",
  "Chloroflexota",
  "Nitrospirota",
  "Myxococcota",
  "Halobacteriota",
  "Bacillota_B"
  )

mge_levels <- c("IS_Tn", "Phage", "CE", "Integron", "ambiguous")
mge_colours <- ggsci::pal_npg()(4) %>% c('grey50')


#################################
### MGE proportions per phyla ###
#################################

### load clustering files
read_recombinase_clustering <- function() {
    d <- tribble(
        ~cluster, ~filename,
        "90_AAI", "90_AAI_cluster.tsv",
        "100_AAI", "100_AAI_cluster.tsv"
        ) %>%
        mutate(
            data = map(filename, ~ here(recombinase_clustering_directory, .) %>% read_tsv(col_names = c("representative", "contig")), show_col_types = FALSE),
        ) %>%
        select(-filename)
    return(d)
}
recombinase_clustering <- read_recombinase_clustering()

### load recombinase all info file
read_recombinase_contig_info <- function() {
    d <- read_tsv(here(recombinase_directory, "recombinase.allinfo.tsv"), show_col_types = FALSE)
    return(d)
}
recombinase_contig_info <- read_recombinase_contig_info()

### filtered contigs with length minimal length cutoff
list_filtered_contigs <- function(recombinase_contig_info = read_recombinase_contig_info()) {
    d <- recombinase_contig_info %>%
        filter(contig_length >= contig_length_cutoff) %>%
        select(contig) %>%
        distinct() # there could be >1 rec in a contig
    return(d)
}
filtered_contigs <- list_filtered_contigs(recombinase_contig_info)

# load MGE contig tracking table
read_mge_to_mags_checkm2 <- function() {
    d <- read_tsv(here(v3_contig_tracking_directory, "mge_to_mags_checkm2.tsv"), show_col_types = FALSE)
    return(d)
}

mge_to_mags_checkm2 <- read_mge_to_mags_checkm2()
read_genome_info_checkm2 <- function() {
    d <- bind_rows(
            read_tsv(here(v2_mags_directory, "gtdbtk.bac120.summary.tsv"), show_col_types = FALSE),
            read_tsv(here(v2_mags_directory, "gtdbtk.ar53.summary.tsv"), show_col_types = FALSE)
            ) %>%
        select(genome = user_genome, taxonomy = classification, red_value)
    return(d)
}

genome_info_checkm2 <- read_genome_info_checkm2()
taxonomy_checkm2 <- genome_info_checkm2 %>%
  separate(taxonomy, sep = ";",
    into = c("domain", "phylum", "class", "order", "family", "genus", "species")
  ) %>%
  select(-red_value)

# load MAG cluster info
read_mag_derep_clusters_checkm2 <- function() {
    mag_path_pattern <- "(?<=/)[^/]*(?=.fna)"
    genome_set_pattern <- "^[^/]*(?=/)"
    d <- read_tsv(
        here(v2_mags_directory, "95_ANI_clusters.tsv"),
        show_col_types = FALSE,
        col_names = c("representative_path", "genome_path")) %>%
        mutate(
            representative = map_chr(representative_path, str_extract, pattern = mag_path_pattern),
            genome = map_chr(genome_path, str_extract, pattern = mag_path_pattern),
            genome_set = map_chr(genome_path, str_extract, pattern = genome_set_pattern)
            ) %>%
        select(representative, genome)
    
    return(d)
}

mag_derep_clusters_checkm2 <- read_mag_derep_clusters_checkm2()
clusters <- recombinase_clustering %>%
  filter(cluster %in% c("100_AAI", "90_AAI")) %>%
  unnest(data) %>%
  pivot_wider(names_from = cluster, values_from = representative) %>%
  rename(recombinase = contig, ANI_100 = `100_AAI`, AAI_90 = `90_AAI`)

### re-process contig tracking, by jiarong
df_contig_tracking_filt <- mge_to_mags_checkm2 %>%
  select(contig, genome_contig) %>%
  group_by(genome_contig) %>%
  filter(row_number()==1) %>% # make sure each MAG contig (genome_contig) has at most 1 match
  ungroup()

mge_to_mags_checkm2_filt <- mge_to_mags_checkm2 %>%
  inner_join(df_contig_tracking_filt)

# contigs with CheckM2 MAG taxonomy
mag_contigs <- mge_to_mags_checkm2_filt %>%
  separate(taxonomy, sep = ";",
    into = c("domain", "phylum", "class", "order", "family", "genus", "species")
    ) %>%
  filter(!is.na(recombinase)) %>%
  left_join(clusters) %>%
  rename(origin2 = type) %>%
  #distinct(contig, origin2, phylum, ANI_100, AAI_90) %>% ### NOTE: changed by jiarong - delete this line; counting total rec not unique here
  left_join(
    recombinase_contig_info %>%
      select(contig, contig_length) %>%
      distinct()
    ) %>%
  inner_join(filtered_contigs)

type_proportions <- mag_contigs %>%
  group_by(phylum, origin2) %>%
  summarise(n = n()) %>%
  mutate(prop = n / sum(n)) %>%
  mutate(n_total = sum(n), prop_total =  sum(prop))

phylum_species_counts <- mag_derep_clusters_checkm2 %>%
  left_join(taxonomy_checkm2, by = c("representative" = "genome")) %>%
  group_by(phylum) %>%
  summarise(
    n_species = unique(representative) %>% length(),
    n_genomes = n()
    )

mag_mge_plot_data <- type_proportions %>%
  left_join(phylum_species_counts) %>%
  mutate(
    phylum = str_remove(phylum, "p__"),
    phylum = ifelse(phylum %in% selected_phyla, phylum, "other"),
    phylum = factor(phylum, levels = phyla_levels),
    origin2 = factor(origin2, levels = mge_levels),
    per_genome = n / n_genomes,
    per_genome_total = n_total / n_genomes
    ) %>%
  filter(phylum != "other")


############
# Fig.2.A
############

mag_mge_plot_data2 <- mag_mge_plot_data %>%
  arrange(desc(per_genome_total)) %>%
  mutate(phylum = as.character(phylum))

phyla_levels2 <- unique(mag_mge_plot_data2$phylum)
mag_mge_plot_data2$phylum <- factor(mag_mge_plot_data2$phylum, levels = phyla_levels2)

mag_mge_plot <-
  mag_mge_plot_data2 %>%
  ggplot(aes(phylum, n)) +
  geom_col(aes(fill = origin2)) +
  geom_text(aes(y = n_total, label = per_genome_total %>% round(0)), size=9/.pt, hjust = 1.2, data = . %>% select(phylum, n_total, per_genome_total) %>% distinct()) +
  coord_flip() +
  scale_fill_manual("", breaks = mge_levels, values = mge_colours) +
  scale_y_reverse(labels = c('0', '100K', '200K', '300K', '400K'), breaks = c(0, 1e+5, 2e+5, 3e+5, 4e+5), limits = c(0, 3.1e+5) %>% rev) +
  scale_x_discrete(position = 'top') +
  guides(fill=guide_legend(nrow = 2, byrow = T)) +
  theme_classic() +
  theme(legend.position = 'bottom') +
  labs(x="", y="MGE recombinase number")

############
# Fig.2.B
###########

plot_alpha_diversity <- function(df, fileprefix, counting = contig, gcounting = n_species) {
  alpha_diversity <- df %>%
    group_by(origin2, phylum) %>%
    nest() %>%
    left_join(phylum_species_counts) %>%
    mutate(
      #richness = map_dbl(data, ~ .x %>% distinct({{counting}}) %>% nrow()), ### NOTE: changed by jiarong - delete
      richness = map_dbl(data, ~ .x %>% nrow()),
      phylum = str_remove(phylum, "p__"),
      phylum = ifelse(phylum %in% phyla_levels, phylum, "other"),
      #phylum = factor(phylum, levels = phyla_levels),
      phylum = factor(phylum, levels = phyla_levels2),
      origin2 = factor(origin2, levels = mge_levels)
      #origin2 = factor(origin2, levels = mge_levels %>% rev)
      ) %>%
    group_by(origin2, phylum) %>%
    summarise(
      n_species = sum({{gcounting}}),
      richness = sum(richness)
      ) %>%
    ungroup() %>%
    mutate(
      richness_per_species = richness / n_species,
      ) %>%
    complete(origin2, phylum, fill = list(richness = 0, richness_per_species = 0)) %>%
    filter(phylum != "other")
  tile_layers <- list(
    geom_tile(),
    #xlab("MGE type"),
    #ylab("Phylum"),
    labs(x=element_blank(), y=element_blank()),
    theme_classic()
    #cowplot::theme_cowplot()
  )

  base_plot <- alpha_diversity %>%
    ggplot(aes(origin2, phylum, fill = richness_per_species)) +
    tile_layers +
    scale_fill_distiller("MGE recombinase \nper genome", palette = "YlGn", direction = 1) +
    #scale_fill_distiller('', palette = "YlGn", direction = 1) +
    #guides(fill = guide_legend(title.position = 'bottom')) +
    guides(x=guide_axis(angle=45), fill=guide_colorbar(title.position = 'top'))

  max_richness = alpha_diversity %>% pull(richness_per_species) %>% max() %>% round(0)
  second_richness = alpha_diversity %>% arrange(desc(richness_per_species)) %>% slice(2) %>% pull(richness_per_species)
  trunc_plot <- alpha_diversity %>%
    mutate(richness_per_species = ifelse(richness_per_species == max(richness_per_species), second_richness, richness_per_species)) %>%
    ggplot(aes(origin2, phylum, fill = richness_per_species)) +
    tile_layers +
    geom_text(aes(x = "IS_Tn", y = "Nitrospirota", label = {{max_richness}}), size = 3) +
    scale_fill_distiller("recombinase clusters \nper species", palette = "YlGn", direction = 1)
  return(base_plot)
}

alpha_per_genomes_plot <- mag_contigs %>% 
    plot_alpha_diversity("mag", gcounting = n_genomes) + 
    theme(
      legend.position = "bottom",
      legend.justification = "centre",
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      )


############
# Fig.2.C
###########
average_mge_per_genome <- (mag_contigs %>% count() %>% pull(n)) / (phylum_species_counts %>% summarise(n_genomes = sum(n_genomes)) %>% pull(n_genomes))
genomes_comp_plot <- mag_contigs %>%
  count(phylum) %>%
  left_join(phylum_species_counts) %>%
  mutate(
    phylum = str_remove(phylum, "p__"),
    phylum = ifelse(phylum %in% phyla_levels, phylum, "Other")
    ) %>%
  mutate(phylum != 'Other') %>%
  ggplot(aes(n_genomes, n, colour = phylum, label = phylum)) +
  geom_abline(slope = average_mge_per_genome, intercept = 0, linetype = "dashed") +
  geom_point() +
  ggrepel::geom_text_repel(size = 9/.pt, data = . %>% filter(phylum != "Other"), force = 3, force_pull = 0.5, min.segment.length = 0.2) +
  scale_color_manual(values = phyla_colours_lines, breaks = phyla_levels, guide = "none") +
  scale_x_continuous(labels = scales::unit_format(unit='K', scale = 1e-3)) +
  scale_y_continuous(labels = scales::unit_format(unit='K', scale = 1e-3)) +
  xlab("Genome number") +
  ylab("MGE recombinase number") +
  theme_classic()


p <- mag_mge_plot + alpha_per_genomes_plot + genomes_comp_plot + plot_layout(widths = c(3.8, 1.6, 2)) + plot_annotation(tag_levels = 'A')
figdir <- here('fig.outdir')
dir.create(figdir)
figfile <- here(figdir, 'fig2.host_lineage.pdf')
ggsave(figfile, p, width = 7.2, height = 4, dpi = 300, device = 'pdf')
