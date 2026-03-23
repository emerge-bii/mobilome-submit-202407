source(here::here('setup.R'))

###
# fig.S4.A
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
       + labs(x='MGE recombinase number in hybrid MAG', y='MGE recombinase\nnumber in short read MAG')
       )

fig2a <- gg

gg_inset <- (ggplot(df %>% filter(origin=='CE'), aes(x=n_long, y=n_short, color=phylum))
             + geom_abline(linetype='dashed')
       + geom_point()
       + scale_color_manual(name='', values = color_pal, guide='none')
       + theme_classic()
       + labs(x=NULL, y=NULL, title=NULL)
       + theme(plot.margin = margin(0,0,0,0))
       )

fig2a <- fig2a + inset_element(gg_inset, 0.72, 0.01, 0.99, 0.99, align_to = 'panel', clip =T, ignore_tag = T)
fig2a <- wrap_elements(full=fig2a)


###
# fig.S4.C
###

### host centric view

rec_allinfo_f <- 'som-data/fig-data/short_vs_long/recombinase.allinfo.add_mag.tsv'
long_vs_short_f <- 'som-data/fig-data/short_vs_long/long_read_comparison_genomes.tsv'
mapping_info_f <- 'som-data/fig-data/short_vs_long/gene_mapping.allinfo.final.add_longonly.tsv'

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

#df_mag %>% write_tsv(file = 'som-data/fig-data/short_vs_long/long_read_comparison_genomes.add_new_name.tsv')

df_rec <- read_tsv(rec_allinfo_f, col_types = cols()) %>%
    dplyr::select(recombinase, origin, OTU_100, OTU_90, mag) %>%
    dplyr::mutate(origin=if_else(stringr::str_detect(origin, ";"), "Ambiguous", origin)) %>%
    dplyr::filter(!is.na(mag)) %>%
    dplyr::filter(mag %in% df_mag$name | mag %in% df_mag$name_short)


color_pal2 <- ggsci::pal_uchicago(palette = 'default')(2)
names(color_pal2) <- c('Short read recovered', 'Short read missed')

gg <- df_mag %>% dplyr::select(name, name_short) %>%
    dplyr::left_join(df_rec %>% select(name=mag, recombinase, origin, OTU_100) %>% tidyr::nest(.by=name, .key="data")) %>%
    dplyr::left_join(df_rec %>% select(name_short=mag, recombinase, origin, OTU_100) %>% tidyr::nest(.by=name_short, .key="data_short")) %>%
    dplyr::mutate(
        shared=map2(data, data_short, \(x, y) length(intersect(x$OTU_100, y$OTU_100))),
        long_only=map2(data, data_short, \(x, y) length(setdiff(x$OTU_100, y$OTU_100))),
        short_only=map2(data, data_short, \(x, y) length(setdiff(y$OTU_100, x$OTU_100)))
    ) %>%
    dplyr::select(name, name_short, shared, long_only, short_only) %>%
    tidyr::unnest(c(shared, long_only, short_only)) %>%
    #dplyr::mutate(recovery_rate=format(round(shared/(long_only+shared), 2), nsmall=2)) %>%
    dplyr::mutate(recovery_rate=(short_only+shared)/(short_only+long_only+shared)) %>%
    dplyr::mutate(recovery_rate=scales::percent(recovery_rate, accuracy = 1L)) %>%
    dplyr::mutate(rec_total = short_only+long_only+shared) %>%
    #dplyr::arrange(desc(recovery_rate)) %>%
    dplyr::arrange(desc(rec_total)) %>%
    dplyr::left_join(df_mag %>% select(name, mag_idx, mag_idx2)) %>%
    #dplyr::mutate(name = stringr::str_c(recovery_rate, mag_idx, sep = ' | ')) %>%
    dplyr::mutate(name = mag_idx2) %>%
    dplyr::mutate(name = factor(name, levels=unique(name))) %>%
    tidyr::pivot_longer(cols = c('shared', 'long_only', 'short_only'), names_to = 'groups', values_to = 'values') %>%
    #dplyr::mutate(groups=if_else(groups=='long_only', 'hybrid_only', groups)) %>%
    dplyr::mutate(groups=if_else(groups=='long_only', 'Short read missed', groups)) %>%
    dplyr::mutate(groups=if_else(groups=='shared' | groups=='short_only', 'Short read recovered', groups)) %>%
    #dplyr::mutate(recovery_rate = if_else(groups=='long_only', recovery_rate, NA)) %>%
    ggplot(aes(y=values, x=name, fill=groups)) + geom_col() + scale_fill_manual(values = color_pal2, guide='none') + 
        geom_text(aes(y=rec_total, label = recovery_rate), vjust=-0.2, size=3) +
        theme_classic() + theme(axis.text.x=ggtext::element_markdown(angle = 45, hjust = 1, vjust = 1), legend.title=element_blank()) + 
        ylim(0, 200) +
        labs(y='Unique recombinases', x='Individual short read - hybrid MAG pair')


options(repr.plot.width=6, repr.plot.height=2.5, repr.plot.res=300)
x_labs <- ggplot_build(gg)$layout$panel_params[[1]]$x$get_labels()
x_labs_new <- stringr::str_remove(x_labs, pattern = regex('MAG[0-9]+\\|'))

gg <- gg + scale_x_discrete(labels=x_labs_new)
fig2c <- gg


###
# fig.S4.B
###

### recombinsase type centric view

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
    dplyr::mutate(groups=if_else(groups=='long_only', 'Short read missed', groups)) %>%
    dplyr::mutate(groups=if_else(groups=='shared' | groups=='short_only', 'Short read recovered', groups)) %>%
    #dplyr::mutate(recovery_rate = if_else(groups=='long_only', recovery_rate, NA)) %>%
    dplyr::filter(origin!='Ambiguous') %>%
    dplyr::filter(origin!='Integron') %>%
    ggplot(aes(y=values, x=origin, fill=groups)) + geom_col() + scale_fill_manual(name='', values = color_pal2) + theme_classic() +
        geom_text(aes(label = recovery_rate, y=rec_total), vjust=-0.2, size=3) +
        ylim(c(0, 500)) +
        guides(fill=guide_legend(nrow = 2, reverse = F)) +
        theme(legend.title=element_blank(), legend.position='bottom', legend.box.spacing = unit(0, 'pt'), legend.margin = margin(0,0,0,0)) + 
        labs(y='Unique recombinases', x='')

fig2b <- gg


options(repr.plot.width=7.2, repr.plot.height=7.5, repr.plot.res=300)
layout <- '
AA
BC
DD
'

fig2d <- ggplot() + theme_void() + cowplot::draw_text('TODO: add flowchart')

p <- fig2d + fig2a + fig2b + fig2c + fig2d + plot_layout(design = layout, widths = c(5.4,1.8), heights = c(2.5,2.7,2.3)) + plot_annotation(tag_levels = 'A')

figdir <- here::here('fig.outdir')
dir.create(figdir)

figfile <- here::here(figdir, 'fig.s04.long_vs_short.17mag.pdf')
ggsave(figfile, p, width = 7.2, height = 7, dpi = 300, device = 'pdf')

#figfile <- here::here(figdir, 'fig.s04.long_vs_short.17mag.png')
#ggsave(figfile, p, width = 7.2, height = 7, dpi = 300, device = 'png')
