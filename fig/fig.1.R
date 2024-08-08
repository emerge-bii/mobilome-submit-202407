source(here::here('setup.R'))

###
# Fig.1.A is from external images
###

###
# Fig.1.B
###

tab_f <- 'som-data/fig-data/20230123_mge_recombinase/recombinase/recombinase.allinfo.tsv'

df <- read_tsv(tab_f, col_types = cols())

# Add distribution of all identified recombinase w/o any length or habitat filtering

# MGE type distribution across habitats

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

options(repr.plot.width=3, repr.plot.height=4, repr.plot.res=300)
gg <- (ggplot(data=df1, aes(x=Habitat, y=perc, fill=origin))
       + geom_col() + scale_fill_manual(values=pal, name='')
       + geom_text(data = df1 %>% group_by(Habitat) %>% filter(row_number() == 1), aes(y=label_pos, label=type_sum), angle=25, size=3)
       + theme_classic()
       + scale_y_continuous(labels = scales::percent_format(accuracy = 1L), limits = c(0, 1.08))
       + guides(fill=guide_legend(ncol = 2))
       + theme(legend.position = 'bottom', 
               legend.box.spacing = unit(0, 'pt'), legend.box.margin = margin(0, 0, 0, 0), 
               legend.spacing = unit(0, 'pt'), legend.margin = margin(0, 0, 0, 0), legend.key.size = unit(8, 'pt'))
       + labs(x='', y='MGE recombinases (%)')
)

fig1b <- gg



#########
# Fig.1.C
#########

# donut plot of rec vs. all genes

tab_f <- 'som-data/fig-data/20230123_mge_recombinase/recombinase/all_gene.ko_cnt_w_clu.tsv'

df_donut <- read_tsv(tab_f, col_types = cols(), col_names = c('gene', 'count')) %>%
    mutate(group = if_else(stringr::str_starts(gene, 'K'), 'Annotated',
                          if_else(stringr::str_starts(gene, 'OTU'), 'Unannotated',
                                 if_else(gene == 'recombinase', 'MGE recombinase', gene)))) %>%
    group_by(group) %>%
    summarise(count = sum(count)) %>%
    ungroup()

df_donut2 <- df_donut

df_donut2 <- df_donut2 %>%
    mutate(perc = round(100 * count / sum(count), digits = 1)) %>%
    mutate(group = factor(group, levels = c('MGE recombinase', 'Annotated', 'Unannotated'))) %>%
    mutate(group2 = stringr::str_c(group, '\n(', perc, '%)'))


total_cnt <- sum(df_donut2$count)
total_cnt <- sprintf('%.1fM', total_cnt / 1e+6)

pal <- c(RColorBrewer::brewer.pal(n = 3, name='BuPu') %>% rev %>% head(n=2), 'grey50')[1:3]

names(pal) <- c('MGE recombinase', 'Annotated', 'Unannotated')

options(repr.plot.width = 2.5, repr.plot.height = 2.5, repr.plot.res = 300)
gg <- ggpubr::ggdonutchart(df_donut2, 'count', label = 'group2', lab.pos = 'out', fill = 'group', color= 'group', lab.font = c(3, 'plain', 'black'), palette = pal) +
    theme(legend.position = 'none') +
    annotate(geom = 'text', x = 0.5, y = 0, label = total_cnt, size = 3)

fig1c <- gg


###
# Fig.1.D
###

tab_f <- 'som-data/fig-data/20230123_mge_recombinase/recombinase/all_gene.ko_cnt_w_clu.tsv'
df <- read_tsv(tab_f, col_types = cols(), col_names = c('gene', 'count'))
df_sub <- df %>% head(n=5) 

col_purple <- RColorBrewer::brewer.pal(n = 3, name='BuPu') %>% rev %>% head(n=1)

df_tmp3 <- df_sub %>% mutate(anno=gene)
df_tmp3 <- df_tmp3 %>% 
    mutate(fill=if_else(anno=='recombinase', col_purple, 'NA'))

color_pal <- df_tmp3$fill
names(color_pal) <- df_tmp3$anno
options(repr.plot.width=2.8, repr.plot.height=2, repr.plot.res=300)
gg <- (ggplot(data=df_tmp3, aes(x=factor(anno, levels=(anno)), y=count, fill=anno)) + geom_col(color='grey20') 
       + theme_classic() 
       + theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1), legend.position = 'none')
       + scale_fill_manual(values = color_pal)
       + scale_y_continuous(breaks = c(500000, 1000000, 1500000, 2000000), labels = c('0.5M', '1M', '1.5M', '2M'))
       + scale_x_discrete(labels = c('MGE\nrecombinase', 'rpoE', 'drrA', 'ABC.CD.P', 'pknB'))
#        + scale_x_discrete(labels = c('MGE recombinase', 'K03088\nRNA polymerase sigma-70 factor', 
#                                      'K09687\nantibiotic transport system ATP-binding protein',
#                                     'K02004\nputative ABC transport system permease', 'K08884\nserine/threonine protein kinase'))
       + labs(x="", y="Gene count")
      )

fig1d <- gg


###
# Fig.1.E
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
    summarise(rp_cnt = median(rp_cnt)) %>%
    ungroup()

    
df_add_scg <- df %>% dplyr::filter(contig_length>=3000 & (!is.na(Habitat))) %>%
  dplyr::mutate(origin=if_else(origin %in% c('Cell', 'IS_Tn;CE', 'IS_Tn;Phage', 'Phage;CE'), 'Other', origin)) %>%
  dplyr::mutate(origin2=if_else(origin2 %in% c('Phage', 'PhageLike', 'PhageOther'), 'Phage', origin2)) %>%
  group_by(origin2, Sample, Habitat, Year, DepthAvg) %>% summarise(count=n()) %>%
  left_join(df_scg %>% select(Sample, rp_cnt), by='Sample') %>% filter(!is.na(rp_cnt)) %>%
  dplyr::mutate(Habitat=factor(Habitat, levels=c('Palsa', 'Bog', 'Fen')), Year=as.factor(Year)) %>%
  dplyr::filter(!is.na(Habitat)) %>%
  dplyr::filter(!is.na(DepthAvg)) %>%
  dplyr::mutate(Depth =  ifelse(DepthAvg >= 1 & DepthAvg < 9, "0-9", 
       ifelse(DepthAvg >= 10 & DepthAvg < 20, "10-19", 
       ifelse(DepthAvg >= 20 & DepthAvg < 30, "20-29",
       ifelse(DepthAvg >= 30 & DepthAvg < 40, "30-39",
       ifelse(DepthAvg >= 40 & DepthAvg < 50, "40-49",
       ifelse(DepthAvg >= 50 & DepthAvg < 60, "50-59",
       ifelse(DepthAvg >= 60 & DepthAvg < 70, "60-69",
       ifelse(DepthAvg >= 70 & DepthAvg < 80, "70-79",
       "80+"))))))))) %>% 
  dplyr::mutate(mge_cnt=count, scg_cnt=rp_cnt)


###
# MAG rec per genome
###

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
  rename(recombinase = contig, AAI_100 = `100_AAI`, AAI_90 = `90_AAI`)


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

khedkar_rec_per_genome_f <- here("som-data", "fig-data", "khedkar_data", "mge_bins_per_genome_final.txt.gz")

df_khedkar <- read_tsv(khedkar_rec_per_genome_f, col_types = cols()) %>% 
    select(-Hotspot, -Cellular) %>% 
    mutate(mge_per_cell = rowSums(across(where(is.numeric)))) %>%
    mutate(source = 'khedkar et al.') %>%
    select(source, mge_per_cell)


df_add_scg8 <- df_add_scg %>% 
    group_by(Sample) %>% summarise(mge_cnt=sum(mge_cnt), scg_cnt=mean(scg_cnt)) %>% # each MGE orign per sample are in different rows so need sum(mge_cnt), scg has already been summed
    filter(scg_cnt >=5) %>%  # NOTE: require at least 5 single copy core gene since low sequencing depth samples are not reliable
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
       #+ ggpubr::stat_compare_means(aes(label = after_stat(p.signif)), comparisons = my_comparisons, vjust = 0.5)
       + theme_classic()
       + scale_x_discrete(labels=c('Stordalen\nMire\ncontig', 'Stordalen\nMire\nMAG', 'Isolate\ngenome\nKhedkar et al.'))
       + scale_y_continuous(breaks = c(0, 10, 100, 1000), trans = scales::pseudo_log_trans(base = 10))
       + annotation_logticks(sides = 'l')
       + labs(x='', y='MGE recombinases\nper genome')
       )

fig1e <- gg


########
# Fig.1.F
########

# short vs. hybrid MAG to get short read MAG rec recovery rate

rec_allinfo_f <- 'som-data/fig-data/short_vs_long/recombinase.allinfo.add_mag.tsv'
long_vs_short_f <- 'som-data/fig-data/short_vs_long/long_read_comparison_genomes.tsv'
mapping_info_f <- 'som-data/fig-data/short_vs_long/gene_mapping.allinfo.final.add_longonly.tsv'


old_taxon_vec <- c('Actinobacteriota')
old2new_taxon <- c('Actinomycetota')
names(old2new_taxon) <- old_taxon_vec

df_taxa <- read_tsv(mapping_info_f, col_types = cols()) %>%
    dplyr::select(name=genome, classification) %>%
    dplyr::distinct() %>%
    tidyr::separate(classification, sep = ';', into = c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species')) %>%
    dplyr::mutate(phylum = stringr::str_remove(phylum, pattern = 'p__')) %>%
    dplyr::mutate(phylum = if_else(phylum %in% old_taxon_vec, old2new_taxon[phylum], phylum))


df_taxa$color <- color_pal[df_taxa$phylum]

df_mag <- read_tsv(long_vs_short_f, col_types = cols()) %>%
    dplyr::mutate(mag_idx=sprintf('MAG%02d', row_number())) %>%
    dplyr::left_join(df_taxa, by = c('name')) %>%
    dplyr::group_by(phylum) %>%
    dplyr::arrange(mag_idx, .by_group = T) %>%
    dplyr::mutate(phylum_idx = row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(mag_idx2=stringr::str_c(mag_idx, '|', phylum, phylum_idx)) %>%
    dplyr::mutate(mag_idx2=stringr::str_c('<span style = \"color:', color, '\">', mag_idx2, '</span>', sep=' '))

#df_mag %>% write_tsv(file = 'som-data/fig-data/short_vs_long/17mag/long_read_comparison_genomes.add_new_name.tsv')

df_rec <- read_tsv(rec_allinfo_f, col_types = cols()) %>%
    dplyr::select(recombinase, origin, OTU_100, OTU_90, mag) %>%
    dplyr::mutate(origin=if_else(stringr::str_detect(origin, ";"), "Ambiguous", origin)) %>%
    dplyr::filter(!is.na(mag)) %>%
    dplyr::filter(mag %in% df_mag$name | mag %in% df_mag$name_short)


color_pal2 <- ggsci::pal_uchicago(palette = 'default')(2)
names(color_pal2) <- c('Short read\nrecovered', 'Short read\nmissed')

gg <- df_mag %>% dplyr::select(name, name_short) %>%
    dplyr::left_join(df_rec %>% select(name=mag, recombinase, origin, OTU_100) %>% tidyr::nest(.by=name, .key="data")) %>%
    dplyr::left_join(df_rec %>% select(name_short=mag, recombinase, origin, OTU_100) %>% tidyr::nest(.by=name_short, .key="data_short")) %>%
    tidyr::unnest(data) %>%
    tidyr::nest(data=c(recombinase, OTU_100)) %>%
    dplyr::mutate(data_short=map2(data_short, origin, \(x, y) x %>% dplyr::filter(origin==y))) %>%
    dplyr::mutate(
        shared=map2(data, data_short, \(x, y) length(intersect(x$OTU_100, y$OTU_100))),
        long_only=map2(data, data_short, \(x, y) length(setdiff(x$OTU_100, y$OTU_100))),
        short_only=map2(data, data_short, \(x, y) length(setdiff(y$OTU_100, x$OTU_100)))
    ) %>%
    dplyr::select(name, name_short, origin, shared, long_only, short_only) %>%
    tidyr::unnest(c(shared, long_only, short_only)) %>%
    dplyr::group_by(origin) %>% 
    dplyr::summarize(shared=sum(shared), long_only=sum(long_only), short_only=sum(short_only)) %>%
    #dplyr::mutate(recovery_rate=format(round(shared/(long_only+shared), 2), nsmall=2)) %>%
    dplyr::mutate(recovery_rate=(short_only+shared)/(short_only+long_only+shared)) %>%
    dplyr::mutate(recovery_rate=scales::percent(recovery_rate, accuracy = 1L)) %>%
    dplyr::mutate(rec_total = short_only+long_only+shared) %>%
    dplyr::arrange(desc(rec_total)) %>%
    dplyr::mutate(origin = factor(origin, levels = c('IS_Tn', 'Phage', 'CE', 'Integron'))) %>%
    #dplyr::mutate(origin = stringr::str_c(recovery_rate, origin, sep = ' | ')) %>%
    tidyr::pivot_longer(cols = c('shared', 'long_only', 'short_only'), names_to = 'groups', values_to = 'values') %>%
    dplyr::mutate(groups=if_else(groups=='long_only', 'Short read\nmissed', groups)) %>%
    dplyr::mutate(groups=if_else(groups=='shared' | groups=='short_only', 'Short read\nrecovered', groups)) %>%
    #dplyr::mutate(recovery_rate = if_else(groups=='long_only', recovery_rate, NA)) %>%
    dplyr::filter(origin!='Ambiguous') %>%
    dplyr::filter(origin!='Integron') %>%
    ggplot(aes(y=values, x=origin, fill=groups)) + geom_col() + scale_fill_manual(name='', values = color_pal2) + theme_classic() +
        geom_text(aes(label = recovery_rate, y=rec_total), vjust=-0.2, size=3) +
        ylim(c(0, 500)) +
        guides(fill=guide_legend(nrow = 2, reverse = F)) +
        theme(legend.title=element_blank(), legend.position='bottom', 
              legend.box.spacing = unit(0, 'pt'), legend.margin = margin(0,0,0,0), legend.key.size = unit(8, 'pt')) + 
        labs(y='Unique MGE\nrecombinases', x='')


leg <- cowplot::get_legend(gg)
gg2 <- gg + theme(legend.position='none')
gg3 <- gg2 + inset_element(leg, 0.6, 0.8, 0.75, 0.99, align_to = 'panel', clip=T, ignore_tag = T)

fig1f <- gg3

layout <- "
AAA
CBE
DBF
"

gg <- ggplot() + theme_void()


fig1b.leg <- ggpubr::as_ggplot(ggpubr::get_legend(fig1b))
fig1b.wo_leg <- fig1b & theme(legend.position = 'none') & labs(tag = '')
fig1b.new <- wrap_elements(panel=fig1b.leg) / fig1b.wo_leg + plot_layout(heights = c(1,4))

p <- gg + 
    fig1b.new +
    wrap_elements(full=fig1c) + 
    fig1d +  fig1e + (fig1f & labs(tag = 'F', y='Unique MGE\nrecombinases')) + 
    plot_layout(design=layout, heights = c(1.5,0.6,0.4), widths = c(1,1,1.4)) + 
    plot_annotation(tag_levels = list(c('A', 'D', '', 'B', 'C', 'E', 'F')))

figdir <- here::here('fig.outdir')
dir.create(figdir)
figfile <- here::here(figdir, 'fig.1.rec_overview.pdf')
ggsave(figfile, p, width = 7.2, height = 7, dpi = 300, device = 'pdf')
#figfile <- here::here(figdir, 'fig.1.rec_overview.png')
#ggsave(figfile, p, width = 7.2, height = 7, dpi = 300, device = 'png')
