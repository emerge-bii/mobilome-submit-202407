source(here::here('setup.R'))

### metaT on recomibnase OTU stability or activity
### cov_0d9

allinfo_tab <- '../metat.gene_cnt.tsv'
rp_gene_lst_f <- 'som-data/fig-data/contig_taxa/rp/rp.ko.bac_arc_shared.list'
rec_f <- 'som-data/mge_recombinase.tsv'
clust_f <- 'som-data/mge_recombinase.clustering.100aai.tsv'
active_gene_list_f <- 'som-data/fig-data/metat/add-hiseq-shared/all.sample.cov_0d9.active.list'
sanity_check_tab <- 'som-data/sample.metat.qc.tsv'

stable <- 'N'
unstable <- 'Y'

df_sanity_check <- read_tsv(file=sanity_check_tab, col_types = cols()) %>%
    #dplyr::filter(read_forward_mapped >=0.8 | read_reverse_mapped >=0.8)
    dplyr::filter(read_reverse_mapped >=0.8)

sample_vec <- df_sanity_check$Sample
df_rec <- read_tsv(rec_f, col_types = cols())

df <- read_tsv(clust_f, col_types=cols()) %>%
    dplyr::inner_join(df_rec, c('mem'='recombinase')) %>%
    dplyr::filter(Sample %in% sample_vec)

active_gene_vec <- readLines(active_gene_list_f)
active_rec_vec <- intersect(active_gene_vec, df$mem)

# partial is meant to be "active"; used here for re-using the code

df <- df %>% 
    dplyr::filter(origin != 'Cell') %>%
    dplyr::mutate(origin = if_else(stringr::str_detect(origin, ';'), 'ambiguous', origin)) %>%
    dplyr::filter(origin != 'ambiguous') %>%
    dplyr::mutate(Habitat = factor(Habitat, levels = c('Palsa', 'Bog', 'Fen')), Year=as.factor(Year)) %>%
    dplyr::mutate(partial=if_else(mem %in% active_rec_vec, 1, 0)) %>%
    dplyr::filter(!is.na(DepthAvg)) %>%
    dplyr::mutate(DepthLumping =  ifelse(DepthAvg >= 0 & DepthAvg < 10, "0-9", 
                           ifelse(DepthAvg >= 10 & DepthAvg < 20, "10-19", 
                           ifelse(DepthAvg >= 20 & DepthAvg < 30, "20-29",
                           ifelse(DepthAvg >= 30 & DepthAvg < 40, "30-39",
                           ifelse(DepthAvg >= 40 & DepthAvg < 50, "40-49",
                           ifelse(DepthAvg >= 50 & DepthAvg < 60, "50-59",
                           ifelse(DepthAvg >= 60 & DepthAvg < 70, "60-69",
                           ifelse(DepthAvg >= 70 & DepthAvg < 80, "70-79",
                           "80+")))))))))


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
       + geom_col()
       #+ geom_text(data=df_cnt, mapping=aes(x=origin, y=max(df_all_summary %>% pull(carriage_type_prop)) * 1.05, label=n, fill=unstable), size=1.5)
       + theme_classic()
       + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
       + scale_fill_brewer(palette = 'Dark2', name=element_blank(), guide='none')
       + scale_y_percent(accuracy = 0.1L)
       + labs(x='', y="active OTUs (%)")
       )


### by Year

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
gg <- (ggplot(data=df_dummy_summary, aes(x=Year, y=carriage_type_prop))
       + geom_col(fill='white', color='grey20')
       #+ geom_text(data=df_cnt, mapping=aes(x=Year, y=max(df_dummy_summary %>% pull(carriage_type_prop)) * 1.05, label=n, fill=unstable), size=1.5)
       + facet_wrap(~origin, ncol=1, scales='free_y')
       + theme_classic()
       + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
       + scale_y_percent(accuracy = 0.1)
       + labs(x='', y="Active OTUs (%)")
       )

fig.a <- gg


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
gg <- (ggplot(data=df_dummy_summary, aes(x=Habitat, y=carriage_type_prop))
       + geom_col(fill='white', color='grey20')
       #+ geom_text(data=df_cnt, mapping=aes(x=Habitat, y=max(df_dummy_summary %>% pull(carriage_type_prop)) * 1.05, label=n, fill=unstable), size=1.5)
       + facet_wrap(~origin, ncol=1, scales='free_y')
       + theme_classic()
       + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
       + scale_y_percent(accuracy = 0.1)
       #+ labs(x='', y="Active OTUs (%)")
       + labs(x='', y="")
       )

fig.b <- gg


### by Depth

df_dummy <- df %>% 
    dplyr::group_by(OTU, DepthLumping, origin) %>%
    dplyr::summarise(partial=sum(partial), total=n()) %>%
    dplyr::mutate(partial_ratio = partial/total) %>%
    dplyr::mutate(partial_carriage = if_else(partial>0, unstable, stable))


df_dummy_summary <- df_dummy %>% 
    group_by(origin, DepthLumping, partial_carriage) %>% 
    summarise(carriage_type_cnt=n()) %>% 
    mutate(carriage_type_prop=carriage_type_cnt/sum(carriage_type_cnt)) %>%
    filter(partial_carriage==unstable)


options(repr.plot.width=2, repr.plot.height=6, repr.plot.res=300)
df_cnt <- df_dummy %>% dplyr::group_by(DepthLumping, origin) %>% dplyr::summarise(n=n())
gg <- (ggplot(data=df_dummy_summary, aes(x=DepthLumping, y=carriage_type_prop))
       + geom_col(fill='white', color='grey20')
       #+ geom_text(data=df_cnt, mapping=aes(x=DepthLumping, y=max(df_dummy_summary %>% pull(carriage_type_prop)) * 1.05, label=n, fill=unstable), size=1.5)
       + facet_wrap(~origin, ncol=1, scales='free_y')
       + theme_classic()
       + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
       + scale_y_percent(accuracy = 0.1L)
       + labs(x='', y="Active OTUs (%)")
       )


### by host domain
df_dummy <- df %>%
    dplyr::mutate(domain = stringr::str_remove_all(domain, 'd__')) %>%
    dplyr::filter(domain != 'Other') %>%
    dplyr::group_by(OTU, domain, origin) %>%
    dplyr::summarise(partial=sum(partial), total=n()) %>%
    dplyr::mutate(partial_ratio = partial/total) %>%
    dplyr::mutate(partial_carriage = if_else(partial>0, unstable, stable))


df_dummy_summary <- df_dummy %>% 
    group_by(origin, domain, partial_carriage) %>% 
    summarise(carriage_type_cnt=n()) %>% 
    mutate(carriage_type_prop=carriage_type_cnt/sum(carriage_type_cnt), , carriage_type_sum=sum(carriage_type_cnt)) %>%
    filter(partial_carriage==unstable)


col_pal <- ggsci::pal_nejm()(2)
names(col_pal) <- c('Archaea', 'Bacteria')
df_dummy_summary$color <- col_pal[df_dummy_summary$domain]
df_dummy_summary <- df_dummy_summary %>%
    mutate(domain = stringr::str_c('<span style = \"color:', color, '\">', domain, '<span>', sep=' '))


options(repr.plot.width=1.5, repr.plot.height=6, repr.plot.res=300)
gg <- (ggplot(data=df_dummy_summary, aes(x=domain, y=carriage_type_prop))
       + geom_col(fill='white', color='grey20')
       #+ geom_text(mapping=aes(x=domain, y=max(df_dummy_summary %>% pull(carriage_type_prop)) * 1.05, label=carriage_type_sum, fill=unstable), size=7/.pt)
       + facet_wrap(~origin, ncol=1, scales='free_y')
       + theme_classic()
       + theme(axis.text.x = ggtext::element_markdown(angle=45, hjust=1, vjust=1))
       + scale_y_percent()
       + labs(x='', y="Active OTUs (%)")
       )


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
    ungroup() %>%
    dplyr::filter(carriage_type_sum>=10) %>%
    filter(partial_carriage==unstable)


col_pal <- ggsci::pal_nejm()(2)
names(col_pal) <- c('Archaea', 'Bacteria')
df_dummy_summary$color <- col_pal[df_dummy_summary$domain]
df_dummy_summary <- df_dummy_summary %>%
    mutate(phylum = stringr::str_c('<span style = \"color:', color, '\">', phylum, '<span>', sep=' '))


phylum_levels <- df_dummy_summary %>% group_by(phylum) %>%
    summarise(carriage_type_prop = sum(carriage_type_prop)) %>%
    arrange(desc(carriage_type_prop)) %>%
    ungroup() %>%
    pull(phylum)


df_dummy_summary <- df_dummy_summary %>%
    mutate(phylum = factor(phylum, levels = phylum_levels))


options(repr.plot.width=7, repr.plot.height=6, repr.plot.res=300)
gg <- (ggplot(data=df_dummy_summary, aes(x=phylum, y=carriage_type_prop))
       + geom_col(fill='white', color='grey20')
       #+ geom_text(aes(x=phylum, y=max(df_dummy_summary %>% pull(carriage_type_prop)) * 1.05, label=carriage_type_sum, fill=unstable), size=1.5)
       + facet_wrap(~origin, ncol=1, scales = 'free_y')
       + theme_classic()
       + theme(axis.text.x = ggtext::element_markdown(angle=45, hjust=1, vjust=1))
       + scale_y_percent()
       + labs(x='', y="Active OTUs (%)")
       )

fig.c <- gg





### PR include all RP gene shared between bac and arc

gene_tab_vec <- c('som-data/fig-data/contig_taxa/rp/all_info.all_gene.fix_domain.tsv')
active_gene_list_f <- 'som-data/fig-data/metat/add-hiseq-shared/all.sample.cov_0d9.active.list'
rec_tab_f <- 'som-data/mge_recombinase.tsv'
sanity_check_tab <- 'som-data/sample.metat.qc.tsv'
non_rp_gene_vec <- c('IS_Tn', 'Phage')
rp_gene_lst_f <- 'som-data/fig-data/contig_taxa/rp/rp.ko.bac_arc_shared.list'

rp_gene_vec <- readLines(rp_gene_lst_f)
select_gene_vec <- c(non_rp_gene_vec, rp_gene_vec)
df_sanity_check <- read_tsv(file=sanity_check_tab, col_types = cols()) %>%
    #dplyr::filter(read_forward_mapped >=0.8 | read_reverse_mapped >=0.8)
    dplyr::filter(read_reverse_mapped >=0.8)

df_gene_lst <- lapply(gene_tab_vec, read_tsv, col_types = cols())

df <- do.call(rbind, df_gene_lst) %>%
    dplyr::mutate(gene_name=if_else(gene_name=='RP', ko, gene_name)) %>% # replace RP gene_name w/ ko
    dplyr::select(-ko)

df_rec <- read_tsv(rec_tab_f, col_types = cols()) %>% 
    dplyr::select(all_of(c('recombinase','Sample','contig','contig_length','domain','phylum','class','order','family','genus','species','Year','Habitat','DepthAvg','origin2'))) %>%
    dplyr::rename(gene=recombinase, gene_name=origin2)

df <- rbind(df, df_rec) %>%
    dplyr::filter(!is.na(DepthAvg)) %>%
    dplyr::mutate(DepthLumping =  ifelse(DepthAvg >= 0 & DepthAvg < 10, "0-9", 
                           ifelse(DepthAvg >= 10 & DepthAvg < 20, "10-19", 
                           ifelse(DepthAvg >= 20 & DepthAvg < 30, "20-29",
                           ifelse(DepthAvg >= 30 & DepthAvg < 40, "30-39",
                           ifelse(DepthAvg >= 40 & DepthAvg < 50, "40-49",
                           ifelse(DepthAvg >= 50 & DepthAvg < 60, "50-59",
                           ifelse(DepthAvg >= 60 & DepthAvg < 70, "60-69",
                           ifelse(DepthAvg >= 70 & DepthAvg < 80, "70-79",
                           "80+")))))))))



active_gene_list <- readLines(active_gene_list_f)
df <- df %>% dplyr::filter(contig_length>=3000) %>%
    dplyr::filter(Sample %in% df_sanity_check$Sample) %>%
    dplyr::filter(domain!='d__Other') %>%
    dplyr::mutate(domain=stringr::str_remove(domain, 'd__')) %>%
    dplyr::mutate(phylum=stringr::str_remove(phylum, 'p__')) %>%
    dplyr::mutate(Habitat=factor(Habitat, levels = c('Palsa', 'Bog', 'Fen'))) %>%
    dplyr::filter(gene_name %in% select_gene_vec) %>%
    dplyr::mutate(active=if_else(gene %in% active_gene_list, 1, 0))


gene_vec <- df %>% dplyr::pull(gene_name) %>% unique()

##################
### by phylum; sep MGE type

df_by_taxa <- df %>% 
    dplyr::group_by(gene_name, domain, phylum) %>%
    dplyr::summarise(total_cnt=n(), active_cnt=sum(active)) %>%
    #dplyr::filter(domain=='Archaea') %>%
    dplyr::ungroup() %>%
    dplyr::mutate(gene_name=if_else(gene_name %in% rp_gene_vec, 'RP', gene_name)) %>%
    dplyr::group_by(gene_name, domain, phylum) %>% # only summarise for RP here
    dplyr::summarise(total_cnt=mean(total_cnt), active_cnt=mean(active_cnt)) %>% # take average for all RP genes
    dplyr::mutate(gene_name=factor(gene_name, levels = c(non_rp_gene_vec, 'RP'))) %>%
    dplyr::mutate(active_ratio=if_else(total_cnt==0, 0, active_cnt/total_cnt)) %>%
    dplyr::ungroup()


# fitlter small active_cnt and total_cnt based on RP, rare (active) members

df_taxa_base_ratio <- df_by_taxa %>% 
    filter(gene_name == 'RP') %>%
    filter(active_cnt >= 5 & total_cnt >= 5) %>%
    select(domain, phylum, base_active_ratio=active_ratio)


df_by_taxa <- df_by_taxa %>%
    #filter(gene_name != 'RP') %>%
    filter(active_cnt >= 5 & total_cnt >=5) %>%
    inner_join(df_taxa_base_ratio, by = c('domain', 'phylum')) %>%
    mutate(active_ratio_norm = active_ratio/base_active_ratio)


phylum_level_fig6b <- df_by_taxa %>%
    group_by(phylum) %>%
    summarise(active_ratio_norm = mean(active_ratio_norm)) %>%
    arrange(desc(active_ratio_norm)) %>%
    pull(phylum)

df_by_taxa <- df_by_taxa %>%
    mutate(phylum = factor(phylum, levels = phylum_level_fig6b))

phylum_level_fig6b <- df_by_taxa %>%
    group_by(phylum) %>%
    filter(gene_name %in% c('CE', 'Integron', 'IS_Tn', 'Phage')) %>%
    summarise(active_ratio_norm = sum(active_ratio_norm)) %>%
    arrange(desc(active_ratio_norm)) %>%
    pull(phylum)

df_by_taxa <- df_by_taxa %>%
    mutate(phylum = factor(phylum, levels = phylum_level_fig6b))

gg <- (ggplot(data=df_by_taxa %>% filter(gene_name %in% c('CE', 'Integron', 'IS_Tn', 'Phage')), aes(x=phylum, y=active_ratio_norm, fill=domain)) + geom_col() 
       + facet_wrap(~gene_name, ncol=1, scales='free_y')
        +theme_classic()
        +labs(x='', y='Normalized active ratio')
        + ggsci::scale_fill_nejm()
        #+ scale_y_continuous(breaks = c(0.25, 0.5, 0.75, 1))
        #+theme(legend.title=element_blank(), axis.text.x=element_blank())
        + theme(legend.title=element_blank(), axis.text.x=element_text(angle=45, hjust=1, vjust=1))
        )

fig.s6b <- gg

options(repr.plot.width=7.2, repr.plot.height=8.8, repr.plot.res=300)
layout <- '
AB#
DDD'

p <- (fig.a & scale_y_percent(accuracy = 0.1)) + (fig.b & scale_y_percent(accuracy = 0.1) & labs(x='', y="Active OTUs (%)")) + (fig.c & scale_y_percent(accuracy = 0.1)) + 
    plot_layout(design = layout, widths = c(2, 1, 1)) + plot_annotation(tag_levels = 'A')

figdir <- here::here('fig.outdir')
dir.create(figdir)

figfile <- here::here(figdir, 'fig.s09.metat.pdf')
ggsave(figfile, width = 7.2, height = 9, dpi = 300, device = 'pdf')

#figfile <- here::here(figdir, 'fig.s09.metat.png')
#ggsave(figfile, width = 7.2, height = 9, dpi = 300, device = 'png')
