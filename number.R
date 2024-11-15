source(here::here('setup.R'))


######
# MGE recombinase stats
######
rec_f <- here::here('som-data/mge_recombinase.tsv')
rec_df <- read_tsv(rec_f, col_types = cols())

n_rec_total <- nrow(rec_df)
n_rec_uniq <- rec_df %>% pull(OTU) %>% discard(is.na) %>% unique %>% length
sprintf('[INFO] total MGE recombinase number: %d', n_rec_total)
sprintf('[INFO] unique MGE recombinase number: %d', n_rec_uniq)

rec_by_habitat_df <- rec_df %>% pull(Habitat) %>% table %>% tibble::enframe()
print(rec_by_habitat_df)

rec_by_type_df <- rec_df %>% pull(origin2) %>% table %>% tibble::enframe() %>%
    rename(count = value) %>%
    mutate(percentage = count / sum(count))
print(rec_by_type_df)


### MAG related stats
n_contig_ge3kb <- rec_df %>% filter(contig_length>=3000) %>% pull(contig) %>% unique %>% length
n_contig_binned <- rec_df %>% filter(contig_length>=3000 & (!is.na(mag))) %>% pull(contig) %>% unique %>% length
sprintf("[INFO] %% of MGE recombinase encoding contigs binned: %d / %d (%.1f%%)", n_contig_binned, n_contig_ge3kb, n_contig_binned / n_contig_ge3kb * 100)

tmp_rec_df <- rec_df %>% filter(contig_length >= 3000 & (!is.na(mag))) %>% select(contig, phylum, mag_lineage) %>% distinct() %>% separate(mag_lineage, sep=';', into=c('mag_domain', 'mag_phylum', 'mag_class', 'mag_order', 'mag_family', 'mag_genus', 'mag_species')) %>% filter(phylum!='p__Other')

n_total <- tmp_rec_df %>% nrow()
n_matched <- tmp_rec_df %>% filter(phylum == mag_phylum) %>% nrow()
sprintf("[INFO] %% of contig phylum (CAT) matched mag phylum: %d / %d (%.1f%%)", n_matched, n_total, n_matched/n_total*100)

### MGE recombinase per genome
rec_per_genome_mag <- 23 # from Fig.2
rec_recovery_low <- 0.35
rec_recovery_high <- 0.6
sprintf('[INFO] adjusted recombinase per genome from MAG: %0.0f - %0.0f', 23/0.6, 23 / 0.35)

### host lineage stats
n_phylum <- rec_df %>% filter(contig_length>=3000) %>% pull(phylum) %>% unique %>% discard(function(x) x=='p__Other') %>% length
sprintf('[INFO] total host phylum number: %d', n_phylum)


### MGE metabolism stats
integron_f <- here::here('som-data/impact_function.integron.all_gene.tsv')
ce_f <- here::here('som-data/impact_function.ce.all_gene.tsv')
phl_f <- here::here('som-data/impact_function.phl.all_gene.tsv')
is_f <- here::here('som-data/impact_function.is.interrupted_gene.tsv')

integron_kegg_cnt <- read_tsv(integron_f, col_types = cols()) %>%
    filter(!is.na(kegg_hit)) %>%
    filter(stringr::str_detect(kegg_hit, fixed('hypothetical', ignore_case=T))) %>%
    pull(kegg_hit) %>% unique() %>% length
ce_kegg_cnt <- read_tsv(ce_f, col_types = cols()) %>%
    filter(!is.na(kegg_hit)) %>%
    filter(stringr::str_detect(kegg_hit, fixed('hypothetical', ignore_case=T))) %>%
    pull(kegg_hit) %>% unique() %>% length
phl_kegg_cnt <- read_tsv(phl_f, col_types = cols()) %>%
    filter(!is.na(kegg_hit)) %>%
    filter(stringr::str_detect(kegg_hit, fixed('hypothetical', ignore_case=T))) %>%
    pull(kegg_hit) %>% unique() %>% length
is_kegg_cnt <- read_tsv(is_f, col_types = cols()) %>%
    filter(!is.na(kegg_hit)) %>%
    filter(stringr::str_detect(kegg_hit, fixed('hypothetical', ignore_case=T))) %>%
    pull(kegg_hit) %>% unique() %>% length

total_kegg_cnt <- integron_kegg_cnt + ce_kegg_cnt + phl_kegg_cnt + is_kegg_cnt
sprintf('[INFO] total cargo genes and IS interrupted genes with KEGG annotation: %d', total_kegg_cnt)

curated_kegg_cnt <- read_tsv(here::here("som-data/impact_function.curated.tsv"), col_types = cols()) %>% nrow()
sprintf('[INFO] %% of curated KEGG annotations: %d / %d (%.1f%%)', curated_kegg_cnt, total_kegg_cnt, curated_kegg_cnt/total_kegg_cnt*100)

### activity and mobility
nh_pw_f <- here::here('som-data/genomic_neighborhood.pairwise_ani.tsv')
nh_pw_df <- read_tsv(nh_pw_f, col_types=cols())
colnames(nh_pw_df) <- c('seqname1', 'seqname2', 'low_ani', 'high_ani', 'OTU')
rec_vec <- rec_df %>% filter(origin2 != 'ambiguous') %>% pull(recombinase)
df_pw <- nh_pw_df %>% filter(seqname1 %in% rec_vec & seqname2 %in% rec_vec)
n_otu_ge2 <- df_pw %>% pull(OTU) %>% unique() %>% length()

sprintf("[INFO] # of unique MGE recombinases (in protein seq) occurred >=2 times: %d", n_otu_ge2)

# partial carriage
par_car_df <- read_tsv(here::here("som-data/partial_carriage.info.tsv"), col_types = cols()) %>% 
    dplyr::filter(Habitat %in% c('Palsa', 'Bog', 'Fen')) %>%
    dplyr::mutate(origin=origin2, Year=as.factor(Year), Habitat=factor(Habitat, levels=c('Palsa', 'Bog', 'Fen'))) %>%
    dplyr::mutate(cohen_d_genome = if_else(cohen_d_upstream >= cohen_d_downstream, cohen_d_upstream, cohen_d_downstream)) %>%
    dplyr::mutate(carriage_ratio_genome = if_else(cohen_d_upstream >= cohen_d_downstream, carriage_ratio_upstream, carriage_ratio_downstream)) %>%
    dplyr::mutate(partial = if_else((cohen_d_genome >= 1 & carriage_ratio_genome<1), 1, 0)) %>%
    dplyr::mutate(genome_mean = if_else(cohen_d_upstream >= cohen_d_downstream, upstream_mean, downstream_mean)) %>%
    dplyr::filter(upstream_mean>=1 | downstream_mean>=1) %>%
    dplyr::filter(origin != 'ambiguous')

n_otu <- par_car_df %>% pull(OTU) %>% unique() %>% length()
sprintf("[INFO] # of unique MGE recombinase (in protein seq) with partial carriage detected: %d",  n_otu)

is_interrupt_df <- read_tsv(here::here("som-data/impact_function.is.interrupted_gene.tsv"), col_types = cols()) %>% rename(recombinase=rec_gene_id)
n_is_interrupt <- is_interrupt_df %>% inner_join(rec_df %>% select(recombinase, OTU)) %>% pull(OTU) %>% unique() %>% length()
n_is_interrupt_partial_analysis <- is_interrupt_df %>% inner_join(par_car_df %>% select(recombinase, OTU)) %>% pull(OTU) %>% unique() %>% length()
n_is_interrupt_partial_carry <- is_interrupt_df %>% inner_join(par_car_df %>% filter(partial == 1) %>% select(recombinase, OTU)) %>% pull(OTU) %>% unique() %>% length()
sprintf("[INFO] # of unique IS recombinase interrupting genes: %d", n_is_interrupt)
sprintf("[INFO] %% of unique IS recombinase interrupting genes that are partial carried: %.1f%% (%d / %d)", n_is_interrupt_partial_carry / n_is_interrupt_partial_analysis * 100, n_is_interrupt_partial_carry, n_is_interrupt_partial_analysis)

active_rec_vec <- readLines("som-data/fig-data/metat/add-hiseq-shared/all.sample.cov_0d9.active.list")
metat_samples <- read_tsv(here::here("som-data/sample.metat.qc.tsv"), col_types=cols()) %>% filter(read_reverse_mapped >= 0.8) %>% pull(Sample)
metat_otu_df <- rec_df %>% filter(Sample %in% metat_samples) %>% mutate(active=if_else(recombinase %in% active_rec_vec, 1, 0)) %>% group_by(OTU) %>% summarise(active = sum(active))
par_car_otu_df <- par_car_df %>% group_by(OTU) %>% summarise(partial=sum(partial))
nh_otu_df <- df_pw %>% mutate(unstable = if_else(low_ani < 1, 1, 0)) %>% group_by(OTU) %>% summarise(unstable = sum(unstable))

overlap_otu_df <- metat_otu_df %>% inner_join(par_car_otu_df) %>% inner_join(nh_otu_df)
sprintf("[INFO] recombinase OTUs shared across three activity/mobility analyses: %d", overlap_otu_df %>% nrow())
n_active <- overlap_otu_df %>% filter(active != 0) %>% nrow()
n_active_unstable <- overlap_otu_df %>% filter(active != 0) %>% filter(unstable != 0) %>% nrow()
n_active_partial <- overlap_otu_df %>% filter(active != 0) %>% filter(partial != 0) %>% nrow()
sprintf("[INFO] Among %d recombinase OTUs expressed in metaT, %d (%.1f%%) were detected in >1 genomic neighborhoods and %d (%.1f%%) were associated with partial carriage", n_active, n_active_unstable, n_active_unstable/n_active*100, n_active_partial, n_active_partial/n_active*100)

### conclusion

n_is_interrupt <- read_tsv(is_f, col_types=cols()) %>% nrow()
n_is <- rec_by_type_df %>% filter(name == "IS_Tn") %>% pull(count) %>% .[[1]]
sprintf("[INFO] %% of IS_Tn recombinase interrupting genes: %d / %d (%.1f%%)", n_is_interrupt, n_is, n_is_interrupt/n_is*100)
