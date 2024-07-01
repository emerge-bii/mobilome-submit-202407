source(here::here('setup.R'))

### 
# Fig.1.A is from external images
###



###
# Fig.1.B
###

tab_f <- 'som-data/fig-data/20230123_mge_recombinase/recombinase/recombinase.allinfo.tsv'
df <- read_tsv(tab_f, col_types = cols())

##### add distribution of all identified recombinase w/o any length or habitat filtering

### MGE type distribution across habitats

pal <- ggpubr::get_palette(palette='npg', 4)
pal <- c(pal, 'grey50')

# sort MGE recombinase by counts
df_tmp <- df %>% 
    dplyr::group_by(origin2) %>%
    dplyr::summarise(count=n()) %>%
    dplyr::mutate(tmp_sum=sum(count)) %>%
    dplyr::mutate(perc=count/tmp_sum) %>%
    dplyr::arrange(desc(perc))


select_vec <- df_tmp %>% dplyr::pull(origin2)

select_vec <- c(select_vec[ !select_vec=='ambiguous' ], 'ambiguous')

names(pal) <- select_vec

df1 <- df %>% 
    dplyr::filter(!is.na(Habitat)) %>%
    group_by(Habitat, origin2) %>% summarise(count=n()) %>%
    mutate(type_sum=sum(count)) %>% mutate(perc=count/type_sum) %>%
    mutate(label_pos=1.02) %>%
    mutate(origin=origin2) %>%
    dplyr::mutate(Habitat=factor(Habitat, levels=c('Palsa', 'Bog', 'Fen'))) %>%
    dplyr::mutate(origin=factor(origin, levels=select_vec)) %>%
    dplyr::filter(!is.na(Habitat)) %>%
    dplyr::filter(Habitat!='Collapsed Palsa')


gg <- (ggplot(data=df1, aes(x=Habitat, y=perc, fill=origin))
       + geom_col() + scale_fill_manual(values=pal, name='')
       + geom_text(aes(y=label_pos, label=type_sum), angle=0)
       + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
       + theme_classic()
       + ylim(0, 1.04)
       + scale_y_continuous(labels = scales::percent_format(accuracy = 1L))
       + guides(fill=guide_legend(nrow = 2, byrow = T))
       + theme(legend.position = 'top', legend.box.spacing = unit(0, 'pt'),  legend.box.margin = margin(0, 0, 0, 0), legend.spacing = unit(0, 'pt'), legend.margin = margin(0, 0, 0, 0))
       + labs(x='', y='Percentage of MGE recombinases')
)
fig1b <- gg



###
# Fig.1.C
###

# sep rec origin

allinfo_tab <- 'som-data/fig-data/metat/add-hiseq-shared/allinfo.sep_rec_ori.metat_cov_ge0d9.metag_cov_ge0.metat_vs_metag.tsv'
sanity_check_tab <- 'som-data/fig-data/metat/add-hiseq-shared/all.sample.read_mapped_dirseq.contig_ge10k.stats'

df <- read_tsv(file=allinfo_tab, col_types = cols()) %>%
    dplyr::mutate(Habitat=factor(Habitat, levels = c('Palsa', 'Bog', 'Fen'))) %>%
    dplyr::mutate(Year=factor(Year)) %>%
    dplyr::mutate(DepthLumping =  ifelse(DepthAvg >= 0 & DepthAvg < 10, "0-9", 
                           ifelse(DepthAvg >= 10 & DepthAvg < 20, "10-19", 
                           ifelse(DepthAvg >= 20 & DepthAvg < 30, "20-29",
                           ifelse(DepthAvg >= 30 & DepthAvg < 40, "30-39",
                           ifelse(DepthAvg >= 40 & DepthAvg < 50, "40-49",
                           ifelse(DepthAvg >= 50 & DepthAvg < 60, "50-59",
                           ifelse(DepthAvg >= 60 & DepthAvg < 70, "60-69",
                           ifelse(DepthAvg >= 70 & DepthAvg < 80, "70-79",
                           NA)))))))))


df_sanity_check <- read_tsv(file=sanity_check_tab, col_types = cols()) %>%
    #dplyr::filter(read_forward_mapped >=0.8 | read_reverse_mapped >=0.8)
    dplyr::filter(read_reverse_mapped >=0.8)

sample_vec <- df_sanity_check$Sample

df <- df %>% dplyr::filter(sample %in% sample_vec) %>%
    dplyr::mutate(ko=if_else(stringr::str_starts(ko, 'rec__'), 'recombinase', ko)) %>% # this line for combining sep rec
    dplyr::mutate(anno=if_else(ko=='recombinase', 'recombinase', anno)) %>% 
    dplyr::mutate(anno=if_else(is.na(anno), ko, anno)) %>%
    dplyr::filter(cnt_metag>=20) %>%
    dplyr::mutate(anno=if_else(anno=='recombinase', 'MGE recombinase',
                              if_else(stringr::str_detect(anno, 'glnA'), 'glnA\nglutamine synthetase',
                                     if_else(stringr::str_detect(anno, 'hupB'), 'hupB\nDNA binding protein HU-beta',
                                            if_else(stringr::str_detect(anno, 'rpoE'), 'rpoE\nRNA polymerase sigma-70 factor',
                                                   if_else(stringr::str_detect(anno, 'cspA'), 'cspA\ncold shock protein', anno))))))


# convert gene names for better plots
df_tmp3 <- df %>% dplyr::group_by(ko, anno) %>% 
    dplyr::summarise(cnt_metag=sum(cnt_metag)) %>% ungroup() %>% 
    dplyr::arrange(desc(cnt_metag)) %>%
    #dplyr::filter(!is.na(anno)) %>%
    dplyr::slice_head(n=5)


df_tmp3 <- df_tmp3 %>% 
    mutate(fill=if_else(anno=='MGE recombinase', 'grey50', 'NA'))
color_pal <- df_tmp3$fill
names(color_pal) <- df_tmp3$anno

gg <- (ggplot(data=df_tmp3, aes(x=factor(anno, levels=(anno)), y=cnt_metag, fill=anno)) + geom_col(color='grey20') 
       + theme_classic() 
       + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1), legend.position = 'none')
       + scale_fill_manual(values = color_pal)
       + scale_y_continuous(labels = scales::comma)
       #+ scale_x_discrete(labels = c('MGE\nrecombinase', 'K03088', 'K09687', 'K08884', 'K02004'))
       + scale_x_discrete(labels = c('MGE\nrecombinase', 'rpoE', 'drrA', 'pknB', 'ABC.CD.P'))
#        + scale_x_discrete(labels = c('MGE recombinase', 'K03088\nRNA polymerase sigma-70 factor', 
#                                      'K09687\nantibiotic transport system ATP-binding protein', 'K08884\nserine/threonine protein kinase',
#                                     'K02004\nputative ABC transport system permease'))
       + labs(x="", y="Gene count")
       #+ coord_flip()
      )
fig1c <- gg


###
# Fig.1.D
###

tab_f <- here::here('som-data', 'fig-data', '20230123_mge_recombinase', 'recombinase', 'recombinase.allinfo.tsv')
scg_tab_f <- here::here('som-data', 'fig-data', 'contig_taxa', 'rp', 'all_info.all_gene.fix_domain.tsv')
rp_gene_lst_f <- here::here('som-data', 'fig-data', 'contig_taxa', 'rp', 'rp.ko.bac_arc_shared.list')

rp_gene_vec <- readLines(rp_gene_lst_f)
df <- read_tsv(tab_f, col_types = cols())

df_scg <- read_tsv(scg_tab_f, col_types=cols()) %>%
    dplyr::filter(ko %in% rp_gene_vec) %>%
    group_by(Sample, Habitat, Year, DepthAvg, ko) %>% 
    summarise(rp_cnt=n()) %>%
    #summarise(rp_cnt=mean(rp_cnt)) %>%
    summarise(rp_cnt = median(rp_cnt)) %>%
    ungroup()


df_add_scg <- df %>% dplyr::filter(contig_length>=3000 & (!is.na(Habitat))) %>%
  dplyr::mutate(origin=if_else(origin %in% c('Cell', 'IS_Tn;CE', 'IS_Tn;Phage', 'Phage;CE'), 'Other', origin)) %>%
  dplyr::mutate(origin2=if_else(origin2 %in% c('Phage', 'PhageLike', 'PhageOther'), 'Phage', origin2)) %>%
  group_by(origin2, Sample, Habitat, Year, DepthAvg) %>% summarise(count=n()) %>%
  left_join(df_scg %>% select(Sample, rp_cnt), by='Sample') %>% filter(!is.na(rp_cnt)) %>%
  dplyr::mutate(Habitat=factor(Habitat, levels=c('Palsa', 'Bog', 'Fen')), Year=as.factor(Year)) %>%
  dplyr::filter(!is.na(Habitat)) %>%
  dplyr::mutate(Depth =  ifelse(DepthAvg >= 1 & DepthAvg < 9, "0-9", 
       ifelse(DepthAvg >= 10 & DepthAvg < 20, "10-19", 
       ifelse(DepthAvg >= 20 & DepthAvg < 30, "20-29",
       ifelse(DepthAvg >= 30 & DepthAvg < 40, "30-39",
       ifelse(DepthAvg >= 40 & DepthAvg < 50, "40-49",
       ifelse(DepthAvg >= 50 & DepthAvg < 60, "50-59",
       ifelse(DepthAvg >= 60 & DepthAvg < 70, "60-69",
       ifelse(DepthAvg >= 70 & DepthAvg < 80, "70-79",
       NA))))))))) %>% 
    dplyr::filter(!is.na(Depth)) %>%
    dplyr::mutate(mge_cnt=count, scg_cnt=rp_cnt)


####
## MAG rec per genome
####

# set dirctories
v2_mags_directory <- here("som-data", "fig-data", "emerge_mags_v2")
recombinase_directory <- here("som-data", "fig-data", "20230123_mge_recombinase", "recombinase")
v3_contig_tracking_directory <- here("som-data", "fig-data", "contig_tracking_v3")
recombinase_clustering_directory <- here("som-data", "fig-data", "recombinase_clustering_v1")

checkm2_report_f <- here("som-data", "fig-data", "emerge_mags_v2", "checkm2_v1.0.2_quality_report.tsv")

contig_length_cutoff <- 3000

### load clustering files
read_recombinase_clustering <- function() {
    d <- tribble(
        ~cluster, ~filename,
        "100_ANI", "100_ANI_cluster.tsv",
        "70_AAI", "70_AAI_cluster.tsv",
        "75_AAI", "75_AAI_cluster.tsv",
        "80_AAI", "80_AAI_cluster.tsv",
        "85_AAI", "85_AAI_cluster.tsv",
        "90_AAI", "90_AAI_cluster.tsv",
        "95_AAI", "95_AAI_cluster.tsv",
        "99_AAI", "99_AAI_cluster.tsv",
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
  filter(cluster %in% c("100_ANI", "90_AAI")) %>%
  unnest(data) %>%
  pivot_wider(names_from = cluster, values_from = representative) %>%
  rename(recombinase = contig, ANI_100 = `100_ANI`, AAI_90 = `90_AAI`)


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
  #distinct(contig, origin2, phylum, ANI_100, AAI_90) %>% ### NOTE: changed by jiarong - delete this line
  left_join(
    recombinase_contig_info %>%
      select(contig, contig_length) %>%
      distinct()
    )


df_mag_per_cluster <- mag_contigs %>% group_by(MAG) %>%
    summarise(mge_per_cell = n()) %>% ungroup() %>%
    right_join(mag_derep_clusters_checkm2, by=c('MAG'='genome')) %>%
    right_join(
        read_tsv(checkm2_report_f, col_types = cols()) %>% select(MAG=Name, completeness=Completeness, contamination=Contamination),
        by = c('MAG')
    ) %>%
    mutate(mge_per_cell = mge_per_cell * 100/completeness) %>%
    filter(stringr::str_detect(MAG, regex('^201|^33'))) %>% # select data used for rec identification
    mutate(mge_per_cell=if_else(is.na(mge_per_cell), 0, mge_per_cell)) %>%
    mutate(source='mag') %>%
    select(source, mge_per_cell)


###
#  khedkar dataset
###

rec_per_genome_khedkar <- c(
    rep(x = 29, times = 95), rep(53, 89), rep(19, 19), rep(14, 375), rep(32, 570), rep(13, 271),
    rep(33, 501), rep(11, 53), rep(26, 9859), rep(30, 36), rep(8, 240), rep(33, 49), rep(16, 128),
    rep(21, 108), rep(25, 25644), rep(24, 41), rep(35, 1621), rep(16, 47), rep(29, 65), rep(34, 2001),
    rep(34, 3998), rep(27, 27249), rep(8, 2693)
)

df_khedkar <- tibble::tibble(mge_per_cell=rec_per_genome_khedkar, source='khedkar et al.')

df_add_scg8 <- df_add_scg %>% 
    group_by(Sample) %>% summarise(mge_cnt=sum(mge_cnt), scg_cnt=mean(scg_cnt)) %>% # each MGE orign per sample are in different rows so need sum(mge_cnt), scg has already been summed
    filter(scg_cnt >=5) %>%
    mutate(mge_per_cell=mge_cnt/scg_cnt) %>%
    mutate(source='contig') %>%
    select(source, mge_per_cell) %>%
    bind_rows(df_mag_per_cluster, df_khedkar) %>%
    mutate(source=factor(source, levels=c('contig', 'mag', 'khedkar et al.')))


# show the stats needed for main text
#tapply(df_add_scg8$mge_per_cell, df_add_scg8$source, summary)


my_comparisons <- list(c("contig", "mag"), c("mag", "khedkar et al."), c("contig", "khedkar et al."))

options(repr.plot.width=3, repr.plot.height=2.3, repr.plot.res=300)
gg <- (ggplot(data=df_add_scg8, aes(x=source, y=mge_per_cell)) 
       + geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.6) 
       + geom_hline(yintercept = 26, color='red', linetype='dashed')
       + ggpubr::stat_compare_means(aes(label = after_stat(p.signif)), comparisons = my_comparisons, vjust = 0.5)
       + theme_classic()
       + scale_x_discrete(labels=c('Stordalen\nMire\ncontig', 'Stordalen\nMire\nMAG', 'Isolate\ngenome\nKhedkar et al.'))
       + scale_y_continuous(breaks = c(0, 10, 100, 1000), trans = scales::pseudo_log_trans(base = 10))
       + annotation_logticks(sides = 'l')
       + labs(x='', y='MGE recombinase\nnumber per genome')
       )

fig1d <- gg



layout <- "
AAAA
AAAA
BBCC
BBDD
"

gg <- ggplot() + theme_void()
p <- gg + fig1b + fig1c + fig1d + plot_layout(design=layout) + plot_annotation(tag_levels = 'A')

figdir <- here::here('fig.outdir')
dir.create(figdir)

figfile <- here::here(figdir, 'fig1.rec_overview.pdf')
ggsave(figfile, p, width = 7.2, height = 8, dpi = 300, device = 'pdf')

#figfile <- here::here(figdir, 'fig1.rec_overview.png')
#ggsave(figfile, p, width = 7.2, height = 8, dpi = 300, device = 'png')
