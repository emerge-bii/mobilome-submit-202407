source(here::here('setup.R'))

tab_f <- here::here('som-data/mge_recombinase.tsv')
contig_min_len <- 3000

df <- read_tsv(tab_f, col_types = cols()) %>% mutate(group = rec_subfamily, type = rec_family)


###
# fig.S2.B
###

# unique MGE recombinases across Habitats

pal <- ggpubr::get_palette(palette='npg', 4)
pal <- c(pal, 'grey70')

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

### MGE type distribution across habitats
### Add distribution of all identified recombinase w/o any length or habitat filtering

df1_uniq <- df %>% 
    dplyr::filter(Habitat %in% c('Palsa', 'Bog', 'Fen') & !is.na(OTU)) %>%
    dplyr::group_by(OTU) %>%
    dplyr::filter(row_number() == 1) %>% #Pick top row per group
    dplyr::ungroup() %>%
    group_by(Habitat, origin2) %>% summarise(count=n()) %>%
    mutate(type_sum=sum(count)) %>% mutate(perc=count/type_sum) %>%
    mutate(label_pos=1.02) %>%
    mutate(origin=origin2) %>%
    dplyr::mutate(Habitat=factor(Habitat, levels=c('Palsa', 'Bog', 'Fen'))) %>%
    dplyr::mutate(origin=factor(origin, levels=select_vec)) %>%
    dplyr::filter(!is.na(Habitat)) %>%
    dplyr::filter(Habitat!='Collapsed Palsa')


options(repr.plot.width=3.5, repr.plot.height=3.5, repr.plot.res=300)
gg <- (ggplot(data=df1_uniq, aes(x=Habitat, y=perc, fill=origin))
       + geom_col() + scale_fill_manual(values=pal, name='MGE')
       + geom_text(aes(y=label_pos, label=type_sum), angle=0, size=7/.pt)
       + theme_classic()
       + scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1))
       + labs(x='', y='Unique MGE recombinase (%)')
)

fig.s1b <- gg


###
# fig.S2.C
###

# MGE type by Habitat and Year

df1 <- df %>% 
    mutate(Year=as.factor(Year), origin2 = factor(origin2, levels = c('IS_Tn', 'Phage', 'CE', 'Integron', 'ambiguous'))) %>%
    filter(!is.na(Year)) %>%
    dplyr::mutate(Habitat=factor(Habitat, levels=c('Palsa', 'Bog', 'Fen'))) %>%
    dplyr::filter(!is.na(Habitat)) %>%
    group_by(Habitat, Year, origin2) %>% summarise(count=n()) %>%
    mutate(type_sum=sum(count)) %>% mutate(perc=count/type_sum) %>%
    mutate(label_pos=1.07)


options(repr.plot.width=6, repr.plot.height=4, repr.plot.res=300)
gg <- (ggplot(data=df1, aes(x=Year, y=perc, fill=origin2))
       + geom_col()
       + scale_fill_manual(values=pal, name='MGE')
       + geom_text(aes(y=label_pos, label=type_sum), size = 7/.pt, angle=0)
       + facet_wrap(~Habitat, ncol = 1)
       + scale_y_percent(accuracy = 1L, limits = c(0, 1.2), breaks = c(0, 0.25, 0.5, 0.75, 1))
       + theme_classic()
       + labs(x='', y='MGE recombinase (%)')
)

fig.s1c <- gg


###
# fig.S3.A
###

# MGE group (68 subfamily)

pal <- RColorBrewer::brewer.pal(12, 'Paired')[-c(5,6)]
pal <- c(pal, 'grey50')

# add distribution of all identified recombinase w/o any length or habitat filtering
# MGE group distribution across habitats

df_tmp <- df %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(count=n()) %>%
    dplyr::mutate(tmp_sum=sum(count)) %>%
    dplyr::mutate(perc=count/tmp_sum) %>%
    dplyr::arrange(desc(perc))


select_vec <- df_tmp %>% dplyr::pull(group) %>% head(n=10)
select_vec <- c(select_vec, 'Other')
names(pal) <- select_vec

df1 <- df %>% 
    dplyr::mutate(group2 = if_else(group %in% select_vec, group, 'Other')) %>%
    group_by(Habitat, group2) %>% summarise(count=n()) %>%
    mutate(type_sum=sum(count)) %>% mutate(perc=count/type_sum) %>%
    mutate(label_pos=1.02) %>%
    dplyr::mutate(Habitat=factor(Habitat, levels=c('Palsa', 'Bog', 'Fen'))) %>%
    dplyr::mutate(group2=factor(group2, levels=select_vec)) %>%
    dplyr::filter(!is.na(Habitat)) %>%
    dplyr::filter(Habitat!='Collapsed Palsa')


options(repr.plot.width=3.5, repr.plot.height=3.5, repr.plot.res=300)
gg <- (ggplot(data=df1, aes(x=Habitat, y=perc, fill=group2))
       + geom_col() + scale_fill_manual(values=pal, name='Subfamily')
       + geom_text(aes(y=label_pos, label=type_sum), angle=0, size=7/.pt)
       + scale_y_percent(accuracy = 1L, limits = c(0, 1.05))
       + theme_classic()
       + labs(x='', y='MGE recombinase (%)')
)

fig.s1d <- gg


###
# fig.S3.B
###

### subfamily by Habitat and Year

df1 <- df %>% 
    dplyr::mutate(group2 = if_else(group %in% select_vec, group, 'Other')) %>%
    mutate(Year=as.factor(Year)) %>%
    filter(!is.na(Year)) %>%
    dplyr::mutate(Habitat=factor(Habitat, levels=c('Palsa', 'Bog', 'Fen'))) %>%
    dplyr::filter(!is.na(Habitat)) %>%
    group_by(Habitat, Year, group2) %>% summarise(count=n()) %>%
    mutate(type_sum=sum(count)) %>% mutate(perc=count/type_sum) %>%
    mutate(label_pos=1.07) %>%
    dplyr::mutate(group2=factor(group2, levels=select_vec))


options(repr.plot.width=6, repr.plot.height=4, repr.plot.res=300)
gg <- (ggplot(data=df1, aes(x=Year, y=perc, fill=group2))
       + geom_col()
       + scale_fill_manual(values=pal, name='Subfamily')
       + geom_text(aes(y=label_pos, label=type_sum), size = 7/.pt, angle=0)
       + facet_wrap(~Habitat, ncol = 1)
       + scale_y_percent(accuracy = 1L, limits = c(0, 1.2), breaks = c(0, 0.25, 0.5, 0.75, 1))
       + theme_classic()
       + labs(x='', y='MGE recombinase (%)')
)

fig.s1e <- gg


###
# fig.S2.A
###

## alluvial between mge, family, and subfamily

require(ggalluvial)
df1 <- df %>% 
    dplyr::mutate(group2 = if_else(group %in% select_vec, group, 'Other')) %>% # NOTE: select_vec from subfamily
    dplyr::mutate(origin2 = factor(origin2, levels=c('IS_Tn', 'Phage', 'CE', 'Integron', 'ambiguous'))) %>%
    dplyr::mutate(type = factor(type, levels = c('Cas', 'DDE', 'HUH', 'Tyr', 'Ser'))) %>%
    dplyr::select(origin2, type, group2) %>%
    dplyr::count(across(everything()))

pal <- c(ggsci::pal_npg()(4), 'grey70')
names(pal) <- c('IS_Tn', 'Phage', 'CE', 'Integron', 'ambiguous')
options(repr.plot.width = 4.7, repr.plot.height = 4, repr.plot.res = 300)
gg <- (ggplot(data=df1, aes(axis1=origin2, axis2=type, axis3=group2, y=n))
       + scale_x_discrete(limits=c('origin2', 'type', 'group2'), labels=c('MGE', 'Recombinase\nfamily', 'Recombinase\nsubfamily'))
       #+ scale_x_discrete(labels=c('MGE', 'Recombinase\nfamily', 'Recombinase\nsubfamily'))
       + ggalluvial::geom_alluvium(aes(fill=origin2))
       #+ ggalluvial::geom_stratum(aes(fill=factor(origin2, levels=c('IS_Tn', 'Phage', 'CE', 'Integron', 'ambiguous'))))
       + ggalluvial::geom_stratum(aes(fill=origin2), width=0.4)
       #+ ggalluvial::geom_stratum()
       + geom_text(stat='stratum', aes(label=after_stat(stratum)), size=7/.pt, nudge_x = 0)
       + scale_fill_manual(values = pal, name='', guide='none')
       + scale_y_comma()
       + labs(y='MGE recombinase count')
       + theme_minimal())

fig.s1a <- gg


options(repr.plot.width = 7.2, repr.plot.height = 7, repr.plot.res=300)
layout <- '
AAL
BCC
'

p <- list(
    A = fig.s1a,
    B = fig.s1b,
    C = fig.s1c,
    L = guide_area()
) %>%
    wrap_plots(heights = c(1, 1), widths = c(1,1.5), design = layout, guides = 'collect') + plot_annotation(tag_levels = 'A')

figdir <- here::here('fig.outdir')
dir.create(figdir)
figfile <- here::here(figdir, 'fig.s02.rec_overview_by_mge.pdf')
ggsave(figfile, p, width = 7.2, height = 7, dpi = 300, device = 'pdf')

#figfile <- here::here(figdir, 'fig.s02.rec_overview_by_mge.png')
#ggsave(figfile, p, width = 7.2, height = 7, dpi = 300, device = 'png')


### fig.s2
options(repr.plot.width = 7.2, repr.plot.height = 5, repr.plot.res=300)
layout <- '
AB
LL
'

gg <- fig.s1d & theme(legend.margin = margin(0,0,0,0), legend.box.margin = margin(0,0,0,0), legend.key.size = unit(10, 'pt'), legend.position = 'bottom')

leg <- cowplot::get_legend(gg)

p <- list(
    A = fig.s1d & labs(tag = 'A') & theme(legend.position = 'none'),
    B = fig.s1e & labs(y='', tag = 'B') & theme(legend.position = 'none'),
    C = leg
) %>%
    wrap_plots(widths = c(1, 2), heights = c(5, 1), design = layout)


figdir <- here::here('fig.outdir')
figfile <- here::here(figdir, 'fig.s03.rec_overview_by_subfam.pdf')
ggsave(figfile, p, width = 7.2, height = 5, dpi = 300, device = 'pdf')

#figfile <- here::here(figdir, 'fig.s03.rec_overview_by_subfam.png')
#ggsave(figfile, p, width = 7.2, height = 5, dpi = 300, device = 'png')
