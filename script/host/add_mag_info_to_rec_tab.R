
###
# MAG rec per genome
###

# set dirctories
v2_mags_directory <- here("som-data", "fig-data", "emerge_mags_v2")
recombinase_directory <- here("som-data")
v3_contig_tracking_directory <- here("som-data", "fig-data", "contig_tracking_v3")
recombinase_clustering_directory <- here("som-data", "fig-data", "recombinase_clustering_v1")

checkm2_report_f <- here("som-data", "fig-data", "emerge_mags_v2", "checkm2_v1.0.2_quality_report.tsv")

contig_length_cutoff <- 3000

### load clustering files
read_recombinase_clustering <- function() {
    d <- tribble(
        ~cluster, ~filename,
        "90_AAI", "90_AAI_cluster.tsv",
        "100_AAI", "100_AAI_cluster.tsv"
        ) %>%
        mutate(
            data = map(filename, ~ here(recombinase_clustering_directory, .) %>% read_tsv(col_names = c("representative", "contig"), show_col_types = F)),
        ) %>%
        select(-filename)
    return(d)
}
recombinase_clustering <- read_recombinase_clustering()


### load recombinase all info file
read_recombinase_contig_info <- function() {
    d <- read_tsv(here(recombinase_directory, "mge_recombinase.tsv"), show_col_types = FALSE)
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
  rename(recombinase = contig, AAI_100 = `100_AAI`, AAI_90 = `90_AAI`)


### select only MAGs binned from contigs in samples used in this study
### filter redudant MAGs wihtin the same sample
mag_metadata_f <- 'som-data/mag.tsv'
sample_metadata_f <- 'som-data/sample.metadata.tsv' # only for filed samples

sample_metadata_df <- read_tsv(sample_metadata_f, col_types = cols()) %>%
    mutate(folder2 = if_else(folder=='JGI', 'JGI', 'Cronin'))

mag_metadata_df_ori <- read_tsv(mag_metadata_f, col_types = cols()) %>% rename(mag_name=MAG) 

mag_metadata_df <- mag_metadata_df_ori %>%
    filter(stringr::str_detect(SampleID__, 'MainAutochamber')) %>%
    filter(folder %in% c('JGI', 'Cronin_v1', 'Cronin_v2')) %>%
    mutate(folder2 = if_else(folder %in% c('Cronin_v1', 'Cronin_v2'), 'Cronin', 'JGI')) %>%
    left_join(sample_metadata_df, by = c('SampleID__', 'folder2')) %>%
    left_join(
        mag_derep_clusters_checkm2 %>% select(mag_name=genome, mag_cluster=representative)
    ) %>%
    group_by(Sample, mag_cluster) %>%
    arrange(desc(Completeness)) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    filter(!is.na(Sample))

mag2sample_df <- mag_metadata_df %>% select(MAG=mag_name, mag_sample=Sample)


### re-process contig tracking, by jiarong
df_contig_tracking_filt <- mge_to_mags_checkm2 %>%
  select(contig, genome_contig, MAG) %>%
  distinct() %>%
  inner_join(
    recombinase_contig_info %>%
      select(contig, Sample) %>%
      distinct()
    ) %>%
  inner_join(mag2sample_df) %>%
  select(contig, genome_contig, Sample, mag_sample) %>%
  filter(Sample == mag_sample) %>%
  select(contig, genome_contig, Sample)

### add MAG contig match to MGE recombinase master table
df_contig_tracking_filt2 <- mge_to_mags_checkm2 %>%
  select(contig, genome_contig, MAG) %>%
  distinct() %>%
  inner_join(
    recombinase_contig_info %>%
      select(contig, Sample) %>%
      distinct()
    ) %>%
  inner_join(mag2sample_df) %>%
  select(contig, genome_contig, Sample, mag_sample, MAG) %>%
  filter(Sample == mag_sample) %>%
  select(contig, genome_contig, Sample, MAG)

df_rec_master <- recombinase_contig_info %>%
  left_join(
    df_contig_tracking_filt2 %>%
      left_join(
        mag_metadata_df %>% select(MAG=mag_name, Completeness, Classification) 
      ) %>%
      group_by(contig) %>%
      arrange(desc(Completeness)) %>%
      filter(row_number() == 1) %>%
      select(contig, mag=MAG, mag_contig=genome_contig, mag_lineage=Classification)
  )

df_rec_master %>% write_tsv('mge_recombinase.master.tsv')
