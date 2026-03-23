source(here::here('setup.R'))

###
# fig.S6.A
###

tab_f <- 'som-data/fig-data/short_vs_long/long_vs_short.tsv'
allinfo_tab <- 'som-data/fig-data/short_vs_long/gene_mapping.allinfo.final.add_longonly.tsv'

mge_levels <- c('IS_Tn', 'Phage', 'CE', 'Integron', 'Ambiguous')
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
  "Other"
  )

phyla_colours <- ggsci::pal_d3(palette = 'category20', alpha=1)(15)[-4] %>% #remove red, not good w/ green
  head(13) %>% c("grey10")

phyla_colours[7] <- 'grey20' # change grey to black
names(phyla_colours) <- phyla_levels
color_pal <- phyla_colours

old_taxon_vec <- c('Actinobacteriota')
old2new_taxon <- c('Actinomycetota')
names(old2new_taxon) <- old_taxon_vec

df <- read_tsv(tab_f, col_types = cols())

df_info <- read_tsv(allinfo_tab, col_types=cols()) %>%
    dplyr::select(all_of(c('genome', 'classification'))) %>%
    dplyr::distinct() %>%
    tidyr::separate(classification, sep = ';', into = c('domain', 'phylum', 'class', 'order', 'family', 'genus','species')) %>%
    dplyr::mutate(phylum = stringr::str_remove(phylum, 'p__')) %>%
    dplyr::mutate(phylum = if_else(phylum %in% old_taxon_vec, old2new_taxon[phylum], phylum))


df <- df %>% left_join(df_info, by = c('name'='genome')) %>%
    dplyr::mutate(origin = factor(origin, levels=mge_levels)) %>%
    dplyr::filter(origin!='Ambiguous') %>%
    dplyr::filter(origin!='Integron')

options(repr.plot.width=4.2, repr.plot.height=2.8, repr.plot.res=300)
gg <- (ggplot(df, aes(x=n_long, y=n_short, color=phylum)) 
       + geom_abline(linetype='dashed')
       + geom_point()
       + facet_wrap(~origin, nrow=1)
       + scale_color_manual(name='', values = color_pal)
       + guides(color=guide_legend(nrow = 3))
       + theme_classic()
       + theme(legend.text = element_text(size=8), legend.position = 'bottom', legend.box.spacing = unit(0, 'pt'), legend.margin = margin(0,0.2,0,0))
       + labs(x='MGE recombinase number in hybrid MAG', y='MGE recombinase number\nin short read MAG')
       )

fig2a <- gg

gg_inset <- (ggplot(df %>% filter(origin=='CE'), aes(x=n_long, y=n_short, color=phylum))
             + geom_abline(linetype='dashed')
       + geom_point()
       + scale_color_manual(name='', values = color_pal, guide='none')
       + theme_classic()
       + labs(x='', y='')
             + theme(plot.margin = margin(0,0,0,0))
       )

fig2a <- fig2a + inset_element(gg_inset, 0.75, 0.01, 0.99, 0.99, align_to = 'panel', clip =T, ignore_tag = T)
fig.s4.a <- wrap_elements(full=fig2a)


###
# fig.S6.B
###

gene_f <- 'som-data/fig-data/short_vs_long/gene_mapping.allinfo.final.add_longonly.tsv'

df_gene <- read_tsv(gene_f, col_types=cols()) %>%
    dplyr::mutate(long_only=if_else(long_only==TRUE, 'missed', 'covered')) %>%
    tidyr::separate(classification, into = c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'), sep=';') %>%
    dplyr::mutate(genome = stringr::str_c(genome, phylum, sep='\n'))

df_gene2 <- df_gene %>% 
    dplyr::mutate(long_only=if_else(long_only=='missed' & breadth_cov_contig_sheared_sr_100ani >=1, stringr::str_c('missed', 'binning', sep='\n'), if_else(long_only=='missed', stringr::str_c('missed', 'assembly', sep='\n'), long_only)))

mag_f <- 'som-data/fig-data/short_vs_long/long_read_comparison_genomes.add_new_name.tsv'
mag_df <- read_tsv(file = mag_f, col_types = cols()) %>%
    mutate(new_name = stringr::str_c(phylum, phylum_idx, sep = ''))

mag_name <- 'cmr6.MA.201907_S_1_20to24_24'
mag_name_new <- mag_df  %>% filter(name == mag_name) %>% pull(new_name)

options(repr.plot.width=7.2, repr.plot.height=2, repr.plot.res=300)
comp <- list(c('covered', 'missed\nassembly'))
gg1 <- df_gene2 %>% dplyr::filter(stringr::str_detect(genome, pattern = mag_name)) %>% dplyr::select(OTU_90, long_only, origin2) %>% 
    dplyr::group_by(OTU_90) %>% 
    summarise(n=n(), long_only=stringr::str_c(unique(long_only), ';'), n_long_only=map_int(unique(long_only), length), origin2=stringr::str_c(unique(origin2), ';'), n_origin2=lengths(origin2)) %>%
    #dplyr::filter(n_long_only>1 | n_origin2 >1) #%>%
    mutate(long_only=str_remove_all(long_only, ';'), origin2=str_remove_all(origin2, ';')) %>%
    filter(origin2 %in% c('IS_Tn', 'Phage', 'CE')) %>%
    mutate(origin2 = factor(origin2, levels = c('IS_Tn', 'Phage', 'CE', 'Integron'))) %>%
    ggplot(aes(x=long_only, y=n)) + 
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) + 
    geom_hline(yintercept = 1, linetype='dashed', color='red') + 
    facet_wrap(~origin2, nrow = 1) + 
    #ggpubr::stat_compare_means(label = 'p.format', vjust = 1) + 
    #ggpubr::stat_compare_means(ref.group = 'covered', label = 'p.signif', hide.ns = T, vjust = 1) + 
    ggpubr::stat_compare_means(comparisons = comp, label = 'p.signif', hide.ns=F, vjust=1.1) + 
    labs(x='', y='Copy number (90% ANI)', title = mag_name_new) +
    theme_classic()

options(repr.plot.width=3, repr.plot.height=2, repr.plot.res=300)
gg2 <- df_gene2 %>% dplyr::filter(stringr::str_detect(genome, pattern = 'cmr6.MA.201907_S_1_20to24_24')) %>% dplyr::select(OTU_90, long_only, origin2) %>% 
    dplyr::group_by(OTU_90) %>% 
    summarise(n=n(), long_only=stringr::str_c(unique(long_only), ';'), n_long_only=map_int(unique(long_only), length), origin2=stringr::str_c(unique(origin2), ';'), n_origin2=lengths(origin2)) %>%
    mutate(long_only=str_remove_all(long_only, ';'), origin2=str_remove_all(origin2, ';')) %>%
    filter(origin2 %in% c('IS_Tn', 'Phage', 'CE', 'Integron')) %>%
    mutate(origin2 = factor(origin2, levels = c('IS_Tn', 'Phage', 'CE'))) %>%
    ggplot(aes(x=long_only, y=n)) + geom_boxplot() + 
    geom_hline(yintercept = 1, linetype='dashed', color='red') + 
    ggpubr::stat_compare_means(label='p.format', vjust = 1) + 
    labs(x='', y='Copy number (90% ANI)', title = mag_name_new) + 
    theme_classic()

options(repr.plot.width=7.2, repr.plot.height=2, repr.plot.res=300)
gg <- gg2 + gg1 + plot_layout(widths = c(1,4))

fig.s4.b <- gg1


###
# fig.S6.C
###

#cmr6.MA.201907_P_1_20to24_50

mag_name <- 'cmr6.MA.201907_P_1_20to24_50'
mag_name_new <- mag_df  %>% filter(name == mag_name) %>% pull(new_name)

options(repr.plot.width=7.2, repr.plot.height=2, repr.plot.res=300)
gg1 <- df_gene2 %>% dplyr::filter(stringr::str_detect(genome, pattern = 'cmr6.MA.201907_P_1_20to24_50')) %>% dplyr::select(OTU_90, long_only, origin2) %>% 
    dplyr::group_by(OTU_90) %>% 
    summarise(n=n(), long_only=stringr::str_c(unique(long_only), ';'), n_long_only=map_int(unique(long_only), length), origin2=stringr::str_c(unique(origin2), ';'), n_origin2=lengths(origin2)) %>%
    #dplyr::filter(n_long_only>1 | n_origin2 >1) #%>%
    mutate(long_only=str_remove_all(long_only, ';'), origin2=str_remove_all(origin2, ';')) %>%
    filter(origin2 %in% c('IS_Tn', 'Phage', 'CE')) %>%
    mutate(origin2 = factor(origin2, levels = c('IS_Tn', 'Phage', 'CE', 'Integron'))) %>%
    ggplot(aes(x=long_only, y=n)) + 
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) + 
    geom_hline(yintercept = 1, linetype='dashed', color='red') + 
    facet_wrap(~origin2) + 
    #ggpubr::stat_compare_means(label='p.format', vjust = 1) + 
    #ggpubr::stat_compare_means(label='p.signif', vjust = 1, ref.group = 'covered', hide.ns = T) +
    ggpubr::stat_compare_means(comparisons = comp, label = 'p.signif', hide.ns=F, vjust=1.1) + 
    labs(x='', y='Copy number (90% ANI)', title=mag_name_new) + 
    theme_classic()

options(repr.plot.width=3, repr.plot.height=2, repr.plot.res=300)
gg2 <- df_gene2 %>% dplyr::filter(stringr::str_detect(genome, pattern = 'cmr6.MA.201907_P_1_20to24_50')) %>% dplyr::select(OTU_90, long_only, origin2) %>% 
    dplyr::group_by(OTU_90) %>% 
    summarise(n=n(), long_only=stringr::str_c(unique(long_only), ';'), n_long_only=map_int(unique(long_only), length), origin2=stringr::str_c(unique(origin2), ';'), n_origin2=lengths(origin2)) %>%
    #dplyr::filter(n_long_only>1 | n_origin2 >1) #%>%
    mutate(long_only=str_remove_all(long_only, ';'), origin2=str_remove_all(origin2, ';')) %>%
    filter(origin2 %in% c('IS_Tn', 'Phage', 'CE')) %>%
    mutate(origin2 = factor(origin2, levels = c('IS_Tn', 'Phage', 'CE', 'Integron'))) %>%
    ggplot(aes(x=long_only, y=n)) + geom_boxplot() + geom_hline(yintercept = 1, linetype='dashed', color='red') + 
    ggpubr::stat_compare_means(label='p.format', vjust = 1) + 
    labs(x='', y='Copy number (90% ANI)', title=mag_name_new) + 
    theme_classic()

fig.s4.c <- gg1



###
# fig.S6.D
###


# cmr6.MA.201907_E_1_20to24_66

mag_name <- 'cmr6.MA.201907_E_1_20to24_66'
mag_name_new <- mag_df  %>% filter(name == mag_name) %>% pull(new_name)

options(repr.plot.width=7.2, repr.plot.height=2, repr.plot.res=300)
gg1 <- df_gene2 %>% dplyr::filter(stringr::str_detect(genome, pattern = 'cmr6.MA.201907_E_1_20to24_66')) %>% dplyr::select(OTU_90, long_only, origin2) %>% 
    dplyr::group_by(OTU_90) %>% 
    summarise(n=n(), long_only=stringr::str_c(unique(long_only), ';'), n_long_only=map_int(unique(long_only), length), origin2=stringr::str_c(unique(origin2), ';'), n_origin2=lengths(origin2)) %>%
    #dplyr::filter(n_long_only>1 | n_origin2 >1) #%>%
    mutate(long_only=str_remove_all(long_only, ';'), origin2=str_remove_all(origin2, ';')) %>%
    filter(origin2 %in% c('IS_Tn', 'Phage', 'CE')) %>%
    mutate(origin2 = factor(origin2, levels = c('IS_Tn', 'Phage', 'CE', 'Integron'))) %>%
    ggplot(aes(x=long_only, y=n)) + 
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5) + 
    geom_hline(yintercept = 1, linetype='dashed', color='red') + 
    facet_wrap(~origin2, ncol = 1) + 
    #ggpubr::stat_compare_means(label='p.format', vjust=1) + 
    #ggpubr::stat_compare_means(label='p.signif', vjust = 1, ref.group = 'covered', hide.ns = T) +
    ggpubr::stat_compare_means(comparisons = comp, label = 'p.signif', hide.ns=F, vjust=1.1) + 
    labs(x='', y='Copy number (90% ANI)', title = mag_name_new) + 
    theme_classic()


options(repr.plot.width=3, repr.plot.height=2, repr.plot.res=300)
gg2 <- df_gene2 %>% dplyr::filter(stringr::str_detect(genome, pattern = 'cmr6.MA.201907_E_1_20to24_66')) %>% dplyr::select(OTU_90, long_only, origin2) %>% 
    dplyr::group_by(OTU_90) %>% 
    summarise(n=n(), long_only=stringr::str_c(unique(long_only), ';'), n_long_only=map_int(unique(long_only), length), origin2=stringr::str_c(unique(origin2), ';'), n_origin2=lengths(origin2)) %>%
    #dplyr::filter(n_long_only>1 | n_origin2 >1) #%>%
    mutate(long_only=str_remove_all(long_only, ';'), origin2=str_remove_all(origin2, ';')) %>%
    filter(origin2 %in% c('IS_Tn', 'Phage', 'CE')) %>%
    mutate(origin2 = factor(origin2, levels = c('IS_Tn', 'Phage', 'CE', 'Integron'))) %>%
    ggplot(aes(x=long_only, y=n)) + geom_boxplot() + geom_hline(yintercept = 1, linetype='dashed', color='red') + 
    ggpubr::stat_compare_means(label='p.format', vjust=1) + 
    labs(x='', y='Copy number (90% ANI)') + 
    theme_classic()


fig.s4.d <- gg1


options(repr.plot.width=7.2, repr.plot.height=7.5, repr.plot.res=300)
layout <- '
AD
BB
CC
'

p <- fig.s4.a + fig.s4.b + fig.s4.c + fig.s4.d + plot_layout(design = layout, widths = c(4.5, 2.5), heights = c(3, 2, 2)) + plot_annotation(tag_levels = 'A')

figdir <- here::here('fig.outdir')
dir.create(figdir)

figfile <- here::here(figdir, 'fig.s06.short_vs_hybrid.example.pdf')
ggsave(figfile, p, width = 7.2, height = 7.5, dpi=300, device = 'pdf')

#figfile <- here::here(figdir, 'fig.s06.short_vs_hybrid.example.png')
#ggsave(figfile, p, width = 7.2, height = 7.5, dpi=300, device = 'png')
