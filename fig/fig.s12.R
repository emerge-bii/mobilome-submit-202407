source(here::here('setup.R'))


### partial carriage


tab_f <- 'som-data/partial_carriage.info.tsv'
rec_f <- 'som-data/mge_recombinase.tsv'

stable <- 'N'
unstable <- 'Y'
cohen_d_cutoff <- 1

df_rec <- read_tsv(rec_f, col_types = cols())
df <- read_tsv(tab_f, col_types = cols()) %>%
    select(-(Sample:origin2)) %>%
    left_join(df_rec, by = c("recombinase")) %>%
    dplyr::mutate(origin=origin2, Year=as.factor(Year), Habitat=factor(Habitat, levels=c('Palsa', 'Bog', 'Fen'))) %>%
    dplyr::filter(!is.na(Year) & Year != 2010 & Habitat != 'Collapsed Palsa') %>%
    dplyr::mutate(cohen_d_genome = if_else(cohen_d_upstream >= cohen_d_downstream, cohen_d_upstream, cohen_d_downstream)) %>%
    dplyr::mutate(carriage_ratio_genome = if_else(cohen_d_upstream >= cohen_d_downstream, carriage_ratio_upstream, carriage_ratio_downstream)) %>%
    dplyr::mutate(partial = if_else((cohen_d_genome >= cohen_d_cutoff & carriage_ratio_genome<1), 1, 0)) %>%
    dplyr::mutate(genome_mean = if_else(cohen_d_upstream >= cohen_d_downstream, upstream_mean, downstream_mean)) %>%
    dplyr::filter(upstream_mean>=1 | downstream_mean>=1) %>%
    dplyr::filter(origin != 'ambiguous') %>%
    dplyr::filter(!is.na(DepthAvg)) %>%
    mutate(DepthLumping =  ifelse(DepthAvg >= 0 & DepthAvg < 10, "0-9", 
                           ifelse(DepthAvg >= 10 & DepthAvg < 20, "10-19", 
                           ifelse(DepthAvg >= 20 & DepthAvg < 30, "20-29",
                           ifelse(DepthAvg >= 30 & DepthAvg < 40, "30-39",
                           ifelse(DepthAvg >= 40 & DepthAvg < 50, "40-49",
                           ifelse(DepthAvg >= 50 & DepthAvg < 60, "50-59",
                           ifelse(DepthAvg >= 60 & DepthAvg < 70, "60-69",
                           ifelse(DepthAvg >= 70 & DepthAvg < 80, "70-79",
                           "80+")))))))))


###
# all data per MGE type
###

df_all <- df %>% 
    dplyr::group_by(OTU, origin) %>%
    dplyr::summarise(partial=sum(partial), total=n()) %>%
    dplyr::mutate(partial_ratio = partial/total) %>%
    dplyr::mutate(partial_carriage = if_else(partial>0, unstable, stable))


df_all_summary <- df_all %>% 
    group_by(origin, partial_carriage) %>% 
    summarise(carriage_type_cnt=n()) %>% 
    mutate(carriage_type_prop=carriage_type_cnt/sum(carriage_type_cnt)) %>%
    filter(partial_carriage==unstable)


options(repr.plot.width=3, repr.plot.height=4, repr.plot.res=300)
df_cnt <- df_all %>% dplyr::group_by(origin) %>% dplyr::summarise(n=n())
gg <- (ggplot(data=df_all_summary, aes(x=origin, y=carriage_type_prop, fill=partial_carriage))
       + geom_col(color='grey20')
       #+ geom_text(data=df_cnt, mapping=aes(x=origin, y=max(df_all_summary %>% pull(carriage_type_prop)) * 1.05, label=n, fill=unstable), size=1.5)
       + theme_classic()
       + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
       #+ scale_fill_brewer(palette = 'Dark2', name=element_blank(), guide='none')
       + scale_fill_manual(values = c('white'), name=element_blank(), guide='none')
       + scale_y_continuous(labels=scales::percent)
       + labs(x='', y="OTUs with partial carriage by host population (%)")
       )

#fig8a <- gg

df_20k <- df_all_summary %>% mutate(cutoff = 20000)



###
# by Year; fig8b
###

df_dummy <- df %>% 
    dplyr::group_by(OTU, Year, origin) %>%
    dplyr::summarise(partial=sum(partial), total=n()) %>%
    dplyr::mutate(partial_ratio = partial/total) %>%
    dplyr::mutate(partial_carriage = if_else(partial>0, unstable, stable))

df_dummy_summary <- df_dummy %>% 
    group_by(origin, Year, partial_carriage) %>% 
    summarise(carriage_type_cnt=n()) %>% 
    mutate(carriage_type_prop=carriage_type_cnt/sum(carriage_type_cnt)) %>%
    filter(partial_carriage==unstable)

options(repr.plot.width=3, repr.plot.height=6, repr.plot.res=300)
df_cnt <- df_dummy %>% dplyr::group_by(Year, origin) %>% dplyr::summarise(n=n())
gg <- (ggplot(data=df_dummy_summary, aes(x=Year, y=carriage_type_prop, fill=partial_carriage))
       + geom_col(color='grey20')
       #+ geom_text(data=df_cnt, mapping=aes(x=Year, y=max(df_dummy_summary %>% pull(carriage_type_prop)) * 1.05, label=n, fill=unstable), size=1.5)
       + facet_wrap(~origin, ncol=1, scale='free_y')
       + theme_classic()
       + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
       #+ scale_fill_brewer(palette = 'Dark2', name=element_blank(), guide='none')
       + scale_fill_manual(values = c('white'), name=element_blank(), guide='none')
       + scale_y_continuous(labels=scales::percent)
       #+ labs(x='', y="OTUs with partial carriage by host population (%)")
       + labs(x='', y="")
       )

fig8b <- gg


### by Habitat

df_dummy <- df %>% 
    dplyr::group_by(OTU, Habitat, origin) %>%
    dplyr::summarise(partial=sum(partial), total=n()) %>%
    dplyr::mutate(partial_ratio = partial/total) %>%
    dplyr::mutate(partial_carriage = if_else(partial>0, unstable, stable))


df_dummy_summary <- df_dummy %>% 
    group_by(origin, Habitat, partial_carriage) %>% 
    summarise(carriage_type_cnt=n()) %>% 
    mutate(carriage_type_prop=carriage_type_cnt/sum(carriage_type_cnt)) %>%
    filter(partial_carriage==unstable)


options(repr.plot.width=1.5, repr.plot.height=6, repr.plot.res=300)
df_cnt <- df_dummy %>% dplyr::group_by(Habitat, origin) %>% dplyr::summarise(n=n())
gg <- (ggplot(data=df_dummy_summary, aes(x=Habitat, y=carriage_type_prop, fill=partial_carriage))
       + geom_col(color='grey20')
       #+ geom_text(data=df_cnt, mapping=aes(x=Habitat, y=max(df_dummy_summary %>% pull(carriage_type_prop)) * 1.05, label=n, fill=unstable), size=1.5)
       + facet_wrap(~origin, ncol=1, scale='free_y')
       + theme_classic()
       + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
       #+ scale_fill_brewer(palette = 'Dark2', name=element_blank(), guide='none')
       + scale_fill_manual(values = c('white'), name=element_blank(), guide='none')
       + scale_y_continuous(labels=scales::percent)
       #+ labs(x='', y="Percentage of OTUs w/ partial carriage by host")
       + labs(x='', y="")
       )

fig8c <- gg



### by host phylum
df_dummy <- df %>% 
    dplyr::mutate(phylum = stringr::str_remove_all(phylum, 'p__')) %>%
    dplyr::mutate(domain = stringr::str_remove_all(domain, 'd__')) %>%
    dplyr::filter(phylum != 'Other') %>%
    dplyr::group_by(OTU, domain, phylum, origin) %>%
    dplyr::summarise(partial=sum(partial), total=n()) %>%
    dplyr::mutate(partial_ratio = partial/total) %>%
    dplyr::mutate(partial_carriage = if_else(partial>0, unstable, stable))


df_dummy_summary <- df_dummy %>% 
    group_by(origin, domain, phylum, partial_carriage) %>% 
    summarise(carriage_type_cnt=n()) %>% 
    mutate(carriage_type_prop=carriage_type_cnt/sum(carriage_type_cnt), carriage_type_sum=sum(carriage_type_cnt)) %>%
    dplyr::filter(carriage_type_sum>=10) %>%
    filter(partial_carriage==unstable)


col_pal <- ggsci::pal_nejm()(2)
names(col_pal) <- c('Archaea', 'Bacteria')
df_dummy_summary$color <- col_pal[df_dummy_summary$domain]
df_dummy_summary <- df_dummy_summary %>%
    mutate(phylum = stringr::str_c('<span style = \"color:', color, '\">', phylum, '<span>', sep=' '))

phylum_level <- df_dummy_summary %>% 
    group_by(domain, phylum, partial_carriage) %>%
    summarise(carriage_type_prop=sum(carriage_type_prop)) %>% 
    arrange(desc(carriage_type_prop)) %>% 
    pull(phylum)

options(repr.plot.width=7, repr.plot.height=6, repr.plot.res=300)
gg <- (ggplot(data=df_dummy_summary, aes(x=factor(phylum, levels=phylum_level), y=carriage_type_prop, fill=partial_carriage))
       + geom_col(color='grey20')
       #+ geom_text(aes(x=phylum, y=max(df_dummy_summary %>% pull(carriage_type_prop)) * 1.05, label=carriage_type_sum, fill=unstable), size=1.5)
       + facet_wrap(~origin, ncol=1, scale='free_y')
       + theme_classic()
       + theme(axis.text.x = ggtext::element_markdown(angle=45, hjust=1, vjust=1))
       #+ scale_fill_brewer(palette = 'Dark2', name=element_blank(), guide='none')
       + scale_fill_manual(values = c('white'), name=element_blank(), guide='none')
       + scale_y_continuous(labels=scales::percent)
       + labs(x='', y="OTUs with partial carriage by host population (%)")
       )

fig8d <- gg


###
# 10000 bp up/down-stream neighborhood required
###



tab_f <- 'som-data/fig-data/partial_carriage/all_recombinase.coord.min_up_down_10000.rec_vs_contig_depth_cov.add_allinfo.rm_other_rec.dedup.add_cluster.tsv'
rec_f <- 'som-data/mge_recombinase.tsv'

stable <- 'N'
unstable <- 'Y'
cohen_d_cutoff <- 1

df_rec <- read_tsv(rec_f, col_types = cols())
df <- read_tsv(tab_f, col_types = cols()) %>%
    select(-(Sample:origin2)) %>%
    left_join(df_rec, by = c("recombinase")) %>%
    dplyr::mutate(origin=origin2, Year=as.factor(Year), Habitat=factor(Habitat, levels=c('Palsa', 'Bog', 'Fen'))) %>%
    dplyr::filter(!is.na(Year) & Year != 2010 & Habitat != 'Collapsed Palsa') %>%
    dplyr::mutate(cohen_d_genome = if_else(cohen_d_upstream >= cohen_d_downstream, cohen_d_upstream, cohen_d_downstream)) %>%
    dplyr::mutate(carriage_ratio_genome = if_else(cohen_d_upstream >= cohen_d_downstream, carriage_ratio_upstream, carriage_ratio_downstream)) %>%
    dplyr::mutate(partial = if_else((cohen_d_genome >= cohen_d_cutoff & carriage_ratio_genome<1), 1, 0)) %>%
    dplyr::mutate(genome_mean = if_else(cohen_d_upstream >= cohen_d_downstream, upstream_mean, downstream_mean)) %>%
    dplyr::filter(upstream_mean>=1 | downstream_mean>=1) %>%
    dplyr::filter(origin != 'ambiguous') %>%
    dplyr::filter(!is.na(DepthAvg)) %>%
    mutate(DepthLumping =  ifelse(DepthAvg >= 0 & DepthAvg < 10, "0-9", 
                           ifelse(DepthAvg >= 10 & DepthAvg < 20, "10-19", 
                           ifelse(DepthAvg >= 20 & DepthAvg < 30, "20-29",
                           ifelse(DepthAvg >= 30 & DepthAvg < 40, "30-39",
                           ifelse(DepthAvg >= 40 & DepthAvg < 50, "40-49",
                           ifelse(DepthAvg >= 50 & DepthAvg < 60, "50-59",
                           ifelse(DepthAvg >= 60 & DepthAvg < 70, "60-69",
                           ifelse(DepthAvg >= 70 & DepthAvg < 80, "70-79",
                           "80+")))))))))


###
# all data per MGE type
###

df_all <- df %>% 
    dplyr::group_by(OTU, origin) %>%
    dplyr::summarise(partial=sum(partial), total=n()) %>%
    dplyr::mutate(partial_ratio = partial/total) %>%
    dplyr::mutate(partial_carriage = if_else(partial>0, unstable, stable))


df_all_summary <- df_all %>% 
    group_by(origin, partial_carriage) %>% 
    summarise(carriage_type_cnt=n()) %>% 
    mutate(carriage_type_prop=carriage_type_cnt/sum(carriage_type_cnt)) %>%
    filter(partial_carriage==unstable)


options(repr.plot.width=3, repr.plot.height=4, repr.plot.res=300)
df_cnt <- df_all %>% dplyr::group_by(origin) %>% dplyr::summarise(n=n())
gg <- (ggplot(data=df_all_summary, aes(x=origin, y=carriage_type_prop, fill=partial_carriage))
       + geom_col(color='grey20')
       #+ geom_text(data=df_cnt, mapping=aes(x=origin, y=max(df_all_summary %>% pull(carriage_type_prop)) * 1.05, label=n, fill=unstable), size=1.5)
       + theme_classic()
       + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
       #+ scale_fill_brewer(palette = 'Dark2', name=element_blank(), guide='none')
       + scale_fill_manual(values = c('white'), name=element_blank(), guide='none')
       + scale_y_continuous(labels=scales::percent)
       + labs(x='', y="OTUs with partial carriage by host population (%)")
       )

#fig8a <- gg

df_10k <- df_all_summary %>% mutate(cutoff = 10000)



###
# 5000 bp up/down-stream neighborhood required
###

tab_f <- 'som-data/fig-data/partial_carriage/all_recombinase.coord.min_up_down_5000.rec_vs_contig_depth_cov.add_allinfo.rm_other_rec.dedup.add_cluster.tsv'
rec_f <- 'som-data/mge_recombinase.tsv'

stable <- 'N'
unstable <- 'Y'
cohen_d_cutoff <- 1

df_rec <- read_tsv(rec_f, col_types = cols())
df <- read_tsv(tab_f, col_types = cols()) %>%
    select(-(Sample:origin2)) %>%
    left_join(df_rec, by = c("recombinase")) %>%
    dplyr::mutate(origin=origin2, Year=as.factor(Year), Habitat=factor(Habitat, levels=c('Palsa', 'Bog', 'Fen'))) %>%
    dplyr::filter(!is.na(Year) & Year != 2010 & Habitat != 'Collapsed Palsa') %>%
    dplyr::mutate(cohen_d_genome = if_else(cohen_d_upstream >= cohen_d_downstream, cohen_d_upstream, cohen_d_downstream)) %>%
    dplyr::mutate(carriage_ratio_genome = if_else(cohen_d_upstream >= cohen_d_downstream, carriage_ratio_upstream, carriage_ratio_downstream)) %>%
    dplyr::mutate(partial = if_else((cohen_d_genome >= cohen_d_cutoff & carriage_ratio_genome<1), 1, 0)) %>%
    dplyr::mutate(genome_mean = if_else(cohen_d_upstream >= cohen_d_downstream, upstream_mean, downstream_mean)) %>%
    dplyr::filter(upstream_mean>=1 | downstream_mean>=1) %>%
    dplyr::filter(origin != 'ambiguous') %>%
    dplyr::filter(!is.na(DepthAvg)) %>%
    mutate(DepthLumping =  ifelse(DepthAvg >= 0 & DepthAvg < 10, "0-9", 
                           ifelse(DepthAvg >= 10 & DepthAvg < 20, "10-19", 
                           ifelse(DepthAvg >= 20 & DepthAvg < 30, "20-29",
                           ifelse(DepthAvg >= 30 & DepthAvg < 40, "30-39",
                           ifelse(DepthAvg >= 40 & DepthAvg < 50, "40-49",
                           ifelse(DepthAvg >= 50 & DepthAvg < 60, "50-59",
                           ifelse(DepthAvg >= 60 & DepthAvg < 70, "60-69",
                           ifelse(DepthAvg >= 70 & DepthAvg < 80, "70-79",
                           "80+")))))))))


###
# all data per MGE type
###

df_all <- df %>% 
    dplyr::group_by(OTU, origin) %>%
    dplyr::summarise(partial=sum(partial), total=n()) %>%
    dplyr::mutate(partial_ratio = partial/total) %>%
    dplyr::mutate(partial_carriage = if_else(partial>0, unstable, stable))

df_all_summary <- df_all %>% 
    group_by(origin, partial_carriage) %>% 
    summarise(carriage_type_cnt=n()) %>% 
    mutate(carriage_type_prop=carriage_type_cnt/sum(carriage_type_cnt)) %>%
    filter(partial_carriage==unstable)


options(repr.plot.width=3, repr.plot.height=4, repr.plot.res=300)
df_cnt <- df_all %>% dplyr::group_by(origin) %>% dplyr::summarise(n=n())
gg <- (ggplot(data=df_all_summary, aes(x=origin, y=carriage_type_prop, fill=partial_carriage))
       + geom_col(color='grey20')
       #+ geom_text(data=df_cnt, mapping=aes(x=origin, y=max(df_all_summary %>% pull(carriage_type_prop)) * 1.05, label=n, fill=unstable), size=1.5)
       + theme_classic()
       + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
       #+ scale_fill_brewer(palette = 'Dark2', name=element_blank(), guide='none')
       + scale_fill_manual(values = c('white'), name=element_blank(), guide='none')
       + scale_y_continuous(labels=scales::percent)
       + labs(x='', y="OTUs with partial carriage by host population (%)")
       )

#fig8a <- gg

df_5k <- df_all_summary %>% mutate(cutoff = 5000)


###
# 2500 bp up/down-stream neighborhood required
###

tab_f <- 'som-data/fig-data/partial_carriage/all_recombinase.coord.min_up_down_2500.rec_vs_contig_depth_cov.add_allinfo.rm_other_rec.dedup.add_cluster.tsv'
rec_f <- 'som-data/mge_recombinase.tsv'

stable <- 'N'
unstable <- 'Y'
cohen_d_cutoff <- 1

df_rec <- read_tsv(rec_f, col_types = cols())
df <- read_tsv(tab_f, col_types = cols()) %>%
    select(-(Sample:origin2)) %>%
    left_join(df_rec, by = c("recombinase")) %>%
    dplyr::mutate(origin=origin2, Year=as.factor(Year), Habitat=factor(Habitat, levels=c('Palsa', 'Bog', 'Fen'))) %>%
    dplyr::filter(!is.na(Year) & Year != 2010 & Habitat != 'Collapsed Palsa') %>%
    dplyr::mutate(cohen_d_genome = if_else(cohen_d_upstream >= cohen_d_downstream, cohen_d_upstream, cohen_d_downstream)) %>%
    dplyr::mutate(carriage_ratio_genome = if_else(cohen_d_upstream >= cohen_d_downstream, carriage_ratio_upstream, carriage_ratio_downstream)) %>%
    dplyr::mutate(partial = if_else((cohen_d_genome >= cohen_d_cutoff & carriage_ratio_genome<1), 1, 0)) %>%
    dplyr::mutate(genome_mean = if_else(cohen_d_upstream >= cohen_d_downstream, upstream_mean, downstream_mean)) %>%
    dplyr::filter(upstream_mean>=1 | downstream_mean>=1) %>%
    dplyr::filter(origin != 'ambiguous') %>%
    mutate(DepthLumping =  ifelse(DepthAvg >= 0 & DepthAvg < 10, "0-9", 
                           ifelse(DepthAvg >= 10 & DepthAvg < 20, "10-19", 
                           ifelse(DepthAvg >= 20 & DepthAvg < 30, "20-29",
                           ifelse(DepthAvg >= 30 & DepthAvg < 40, "30-39",
                           ifelse(DepthAvg >= 40 & DepthAvg < 50, "40-49",
                           ifelse(DepthAvg >= 50 & DepthAvg < 60, "50-59",
                           ifelse(DepthAvg >= 60 & DepthAvg < 70, "60-69",
                           ifelse(DepthAvg >= 70 & DepthAvg < 80, "70-79",
                           "80+")))))))))



###
# all data per MGE type
###

df_all <- df %>% 
    dplyr::group_by(OTU, origin) %>%
    dplyr::summarise(partial=sum(partial), total=n()) %>%
    dplyr::mutate(partial_ratio = partial/total) %>%
    dplyr::mutate(partial_carriage = if_else(partial>0, unstable, stable))


df_all_summary <- df_all %>% 
    group_by(origin, partial_carriage) %>% 
    summarise(carriage_type_cnt=n()) %>% 
    mutate(carriage_type_prop=carriage_type_cnt/sum(carriage_type_cnt)) %>%
    filter(partial_carriage==unstable)


options(repr.plot.width=3, repr.plot.height=4, repr.plot.res=300)
df_cnt <- df_all %>% dplyr::group_by(origin) %>% dplyr::summarise(n=n())
gg <- (ggplot(data=df_all_summary, aes(x=origin, y=carriage_type_prop, fill=partial_carriage))
       + geom_col(color='grey20')
       #+ geom_text(data=df_cnt, mapping=aes(x=origin, y=max(df_all_summary %>% pull(carriage_type_prop)) * 1.05, label=n, fill=unstable), size=1.5)
       + theme_classic()
       + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
       #+ scale_fill_brewer(palette = 'Dark2', name=element_blank(), guide='none')
       + scale_fill_manual(values = c('white'), name=element_blank(), guide='none')
       + scale_y_continuous(labels=scales::percent)
       + labs(x='', y="OTUs with partial carriage by host population (%)")
       )

#fig8a <- gg
df_2500 <- df_all_summary %>% mutate(cutoff = 2500)

df_merged <- rbind(df_20k, df_10k, df_5k, df_2500)

options(repr.plot.width = 3, repr.plot.height = 4, repr.plot.res=300)
gg_sensitivity <- ggplot(data=df_merged, aes(x=origin, y=carriage_type_prop)) +
    geom_col(fill='white', color='grey20') + facet_wrap(~cutoff, ncol = 1) +
    scale_y_continuous(labels=scales::percent) +
    labs(x='', y="OTUs with partial carriage by host population (%)") +
    theme_classic() +
    guides(x = guide_axis(angle = 45))


options(repr.plot.width=7.2, repr.plot.height=8.6, repr.plot.res=300)
layout <- '
ABC
DDD
'
p <- gg_sensitivity + fig8b + fig8c + fig8d + plot_layout(design = layout) + plot_annotation(tag_levels = 'A')

figdir <- here::here('fig.outdir')
dir.create(figdir)

figfile <- here::here(figdir, 'fig.s12.partial_carriage.pdf')
ggsave(figfile, p, width = 7.2, height = 8.6, dpi = 300, device = 'pdf')

#figfile <- here::here(figdir, 'fig.s12.partial_carriage.png')
#ggsave(figfile, p, width = 7.2, height = 8.6, dpi = 300, device = 'png')
