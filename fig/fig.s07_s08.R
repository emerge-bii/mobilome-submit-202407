source(here::here('setup.R'))

tab_f <- 'som-data/mge_recombinase.tsv'
scg_allinfo_f <- 'som-data/fig-data/contig_taxa/rp/all_info.all_gene.fix_domain.tsv'
contig_min_len <- 3000

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
  ) #%>% rev()

phyla_colours <- ggsci::pal_d3(palette = 'category20', alpha=1)(20)[-4] %>% #remove red, not good w/ green
  head(17) %>% c("grey50")
phyla_colours[7] <- 'grey20' # change grey to black; was postion 8 but becomes 7 after removing red (4)
phyla_colours_lines <- phyla_colours

df <- read_tsv(tab_f, col_types = cols())
df_scg_per_gene <- read_tsv(scg_allinfo_f, col_types=cols())

### add distribution of all identified recombinase w/o any length or habitat filtering
### MGE type distribution across habitats

pal <- ggpubr::get_palette(palette='npg', 4)
pal <- c(pal, 'grey50')

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
    #dplyr::mutate(origin=if_else(origin %in% c('Cell', 'IS_Tn;CE', 'IS_Tn;Phage', 'Phage;CE'), 'Other', origin)) %>%
    group_by(Habitat, origin2) %>% summarise(count=n()) %>%
    mutate(type_sum=sum(count)) %>% mutate(perc=count/type_sum) %>%
    mutate(label_pos=1.05) %>%
    mutate(origin=origin2) %>%
    dplyr::mutate(Habitat=factor(Habitat, levels=c('Palsa', 'Bog', 'Fen'))) %>%
    dplyr::mutate(origin=factor(origin, levels=select_vec)) %>%
    dplyr::filter(!is.na(Habitat)) %>%
    dplyr::filter(Habitat!='Collapsed Palsa')


# MGE type distribution across domain

taxon_level <- 'domain'
df2 <- df %>% 
    dplyr::filter(contig_length>= contig_min_len) %>%
    dplyr::filter((!is.na(Habitat)) & Habitat!='Collapsed Palsa') %>%
    dplyr::mutate(Habitat=factor(Habitat, levels=c('Palsa', 'Bog', 'Fen'))) %>%
    dplyr::group_by(!!rlang::sym(taxon_level), origin2) %>% 
    dplyr::summarise(count=n()) %>%
    dplyr::mutate(tmp_sum=sum(count)) %>%
    dplyr::mutate(perc=count/tmp_sum) %>% dplyr::ungroup() %>%
    dplyr::mutate(origin=origin2) %>%
    dplyr::mutate(origin=factor(origin, levels=select_vec)) %>%
    dplyr::mutate(!!rlang::sym(taxon_level):=stringr::str_remove(!!rlang::sym(taxon_level), 'd__')) %>%
    dplyr::mutate(!!rlang::sym(taxon_level):=if_else(!!rlang::sym(taxon_level)=='Other', 'Unknown', !!rlang::sym(taxon_level))) %>%
    dplyr::mutate(label_pos=1.05)


df3 <- df2 %>% 
    dplyr::group_by(!!rlang::sym(taxon_level)) %>%
    dplyr::summarise(count=sum(count)) %>%
    dplyr::mutate(tmp_sum=sum(count)) %>%
    dplyr::mutate(perc=count/tmp_sum) %>%
    #dplyr::filter(!!rlang::sym(taxon_level)!='Other') %>%
    #dplyr::filter(perc >= 0.01) %>%
    dplyr::arrange(desc(perc))


select_taxon_vec <- df3 %>% dplyr::pull(!!rlang::sym(taxon_level)) %>% head(n=15)

df4 <- df2 %>%
    dplyr::filter(!!rlang::sym(taxon_level) %in% select_taxon_vec) %>%
    dplyr::mutate(!!rlang::sym(taxon_level):=factor(!!rlang::sym(taxon_level), levels=select_taxon_vec))

options(repr.plot.width=5, repr.plot.height=2, repr.plot.res=300)
gg <- (ggplot(data=df4, aes(x=!!rlang::sym(taxon_level), y=perc, fill=origin)) 
       + geom_col() + scale_fill_manual(values=pal, name="MGE type")
       + geom_text(aes(y=label_pos, label=tmp_sum), size=7/.pt, angle=0)
       + scale_y_percent(limits = c(0, 1.1))
       + theme_classic()
       + labs(x='', y='MGE recombinase (%)') + coord_flip()
)

fig.s5a <- gg



### MGE type distribution across phylum

taxon_level <- 'phylum'
df2 <- df %>% 
    dplyr::filter(contig_length>= contig_min_len) %>%
    dplyr::filter((!is.na(Habitat)) & Habitat!='Collapsed Palsa') %>%
    dplyr::mutate(Habitat=factor(Habitat, levels=c('Palsa', 'Bog', 'Fen'))) %>%
    dplyr::group_by(!!rlang::sym(taxon_level), origin2) %>% 
    dplyr::summarise(count=n()) %>%
    dplyr::mutate(tmp_sum=sum(count)) %>%
    dplyr::mutate(perc=count/tmp_sum) %>% dplyr::ungroup() %>%
    dplyr::mutate(origin=origin2) %>%
    dplyr::mutate(origin=factor(origin, levels=select_vec)) %>%
    dplyr::mutate(!!rlang::sym(taxon_level):=stringr::str_remove(!!rlang::sym(taxon_level), 'p__')) %>%
    dplyr::mutate(!!rlang::sym(taxon_level):=if_else(!!rlang::sym(taxon_level)=='Other', 'Unknown', !!rlang::sym(taxon_level))) %>%
    dplyr::mutate(label_pos=1.05)


df3 <- df2 %>% 
    dplyr::group_by(!!rlang::sym(taxon_level)) %>%
    dplyr::summarise(count=sum(count)) %>%
    dplyr::mutate(tmp_sum=sum(count)) %>%
    dplyr::mutate(perc=count/tmp_sum) %>%
    dplyr::filter(!(!!rlang::sym(taxon_level) %in% c('Other', 'Unknown'))) %>%
    #dplyr::filter(perc >= 0.01) %>%
    dplyr::arrange(desc(perc))


select_taxon_vec <- df3 %>% dplyr::pull(!!rlang::sym(taxon_level)) %>% head(n=15)
df4 <- df2 %>%
    dplyr::filter(!!rlang::sym(taxon_level) %in% select_taxon_vec) %>%
    dplyr::mutate(!!rlang::sym(taxon_level):=factor(!!rlang::sym(taxon_level), levels=select_taxon_vec))


options(repr.plot.width=6, repr.plot.height=3, repr.plot.res=300)
gg <- (ggplot(data=df4, aes(x=!!rlang::sym(taxon_level), y=perc, fill=origin)) 
       + geom_col() + scale_fill_manual(values=pal, name='')
       + geom_text(aes(y=label_pos, label=tmp_sum), size=7/.pt, angle=0)
       #+ theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
       + scale_y_percent(limits = c(0, 1.1))
       + theme_classic()
       + labs(x='', y='MGE recombinase (%)') + coord_flip()
)

fig.s5b <- gg


# change the line below to adjust for different data

df2 <- df %>% 
  dplyr::filter(contig_length>=contig_min_len & (!is.na(Habitat))) %>%
  #dplyr::mutate(origin=if_else(origin %in% c('Cell', 'IS_Tn;CE', 'IS_Tn;Phage', 'Phage;CE'), 'Other', origin)) %>%
  group_by(origin2, phylum, Sample, Habitat, Year, DepthAvg) %>% summarise(count=n()) %>% rename(origin=origin2) %>% ungroup()


# ### host using rplB only
# df_host <- df_scg_per_gene %>% dplyr::filter(gene_name=='rplB') %>%
#     dplyr::filter(contig_length>=contig_min_len & (!is.na(Habitat))) %>%
#     group_by(phylum, Sample, Habitat, Year, DepthAvg) %>% summarise(count=n()) %>% ungroup() %>%
#     dplyr::mutate(origin='Host')


### host using all RP genes

df_host <- df_scg_per_gene %>%
    dplyr::filter(contig_length>=contig_min_len & (!is.na(Habitat))) %>%
    group_by(phylum, Sample, Habitat, Year, DepthAvg, gene_name, ko) %>% # gene_name are all "RP" for RP genes
    summarise(count=n()) %>%
    summarise(count=mean(count)) %>% # average among RP genes
    ungroup() %>%
    dplyr::mutate(origin='Host')


df2 <- dplyr::bind_rows(df2, df_host) %>%
  dplyr::mutate(Habitat=factor(Habitat, levels=c('Palsa', 'Bog', 'Fen'))) %>%
  dplyr::filter(!is.na(Habitat) & Habitat!='Collapsed Palsa' & (!is.na(origin))) %>%
  dplyr::filter(!is.na(DepthAvg)) %>%
  dplyr::mutate(Depth =  ifelse(DepthAvg >= 1 & DepthAvg < 9, "0-9", 
       ifelse(DepthAvg >= 10 & DepthAvg < 20, "10-19", 
       ifelse(DepthAvg >= 20 & DepthAvg < 30, "20-29",
       ifelse(DepthAvg >= 30 & DepthAvg < 40, "30-39",
       ifelse(DepthAvg >= 40 & DepthAvg < 50, "40-49",
       ifelse(DepthAvg >= 50 & DepthAvg < 60, "50-59",
       ifelse(DepthAvg >= 60 & DepthAvg < 70, "60-69",
       ifelse(DepthAvg >= 70 & DepthAvg < 80, "70-79",
       "80+")))))))))


### stack bar of phylum whole dataset

df2 <- df2 %>% 
    mutate(origin=if_else(origin=='IS', 'IS_Tn', origin))
df3 <- df2 %>% group_by(origin, phylum) %>% summarise(count=sum(count)) %>% mutate(type_sum=sum(count)) %>%
  dplyr::mutate(perc=count/type_sum)

df4 <- df3 %>% 
    group_by(origin) %>% 
    dplyr::mutate(phylum=stringr::str_remove(phylum, 'p__')) %>%
    dplyr::mutate(phylum=if_else(phylum %in% phyla_levels, phylum, 'Other')) %>%
    dplyr::mutate(phylum=factor(phylum, levels=c(phyla_levels))) %>%
    dplyr::filter(origin!='ambiguous') %>%
    dplyr::mutate(label_pos=1.05)


options(repr.plot.width=5, repr.plot.height=5, repr.plot.res=300)
gg <- (ggplot(data=df4, aes(x=factor(origin, levels=c('Host', 'Integron', 'CE', 'Phage', 'IS_Tn')), y=perc, fill=phylum)) 
       + geom_col() + scale_fill_manual(values=phyla_colours_lines, breaks = phyla_levels)
       #+ geom_text(aes(y=label_pos, label=type_sum), size=7/.pt, angle=0)
       + scale_y_percent(limits = c(0, 1.1))
       + theme_classic()
       + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
       + labs(x='', y='Proportional count')
)

fig.s5c <- gg


########## 

### recombinase per genome; enrichment by phylum


# ### only using rplB
# df_scg <- df_scg_per_gene %>% dplyr::filter(gene_name=='rplB') %>%
#     dplyr::filter(contig_length>=contig_min_len & (!is.na(Habitat))) %>%
#     group_by(Sample, phylum) %>% summarise(scg_cnt=n()) %>% ungroup()

### using all RP genes

df_scg <- df_scg_per_gene %>%
    dplyr::filter(contig_length>=contig_min_len & (!is.na(Habitat))) %>%
    group_by(Sample, phylum, gene_name, ko) %>%
    summarise(scg_cnt=n()) %>% 
    summarise(scg_cnt = mean(scg_cnt)) %>%
    ungroup()


df_add_scg <- df %>% dplyr::filter(contig_length>=contig_min_len & (!is.na(Habitat))) %>%
  dplyr::mutate(origin=if_else(origin %in% c('Cell', 'IS_Tn;CE', 'IS_Tn;Phage', 'Phage;CE'), 'Other', origin)) %>%
  group_by(origin2, Sample, Habitat, Year, DepthAvg, phylum) %>% summarise(mge_cnt=n()) %>% ungroup() %>%
  left_join(df_scg, by=c('Sample', 'phylum')) %>% filter(!is.na(scg_cnt)) %>%
  dplyr::mutate(Habitat=factor(Habitat, levels=c('Palsa', 'Bog', 'Fen')), Year=as.factor(Year)) %>%
  dplyr::filter(!is.na(Habitat)) %>%
  #dplyr::filter(origin2 != 'ambiguous') %>%
  dplyr::mutate(phylum=stringr::str_replace(phylum, 'p__', '')) %>% dplyr::filter(phylum!='Other') %>%
  dplyr::filter(!is.na(DepthAvg)) %>%
  dplyr::mutate(Depth =  ifelse(DepthAvg >= 1 & DepthAvg < 9, "0-9", 
       ifelse(DepthAvg >= 10 & DepthAvg < 20, "10-19", 
       ifelse(DepthAvg >= 20 & DepthAvg < 30, "20-29",
       ifelse(DepthAvg >= 30 & DepthAvg < 40, "30-39",
       ifelse(DepthAvg >= 40 & DepthAvg < 50, "40-49",
       ifelse(DepthAvg >= 50 & DepthAvg < 60, "50-59",
       ifelse(DepthAvg >= 60 & DepthAvg < 70, "60-69",
       ifelse(DepthAvg >= 70 & DepthAvg < 80, "70-79",
       NA)))))))))


df_add_scg7 <- df_add_scg %>% group_by(origin2, phylum) %>%
    summarise(mge_cnt=sum(mge_cnt), scg_cnt=sum(scg_cnt)) %>% 
    dplyr::filter(phylum %in% select_taxon_vec) %>%
    dplyr::mutate(phylum=factor(phylum, levels=select_taxon_vec %>% rev)) %>%
    mutate(mge_per_cell=mge_cnt/scg_cnt) %>%
    mutate(mge_per_cell_mean_aross_phylum = mean(mge_per_cell)) %>%
    ungroup()


df_add_scg8 <- df_add_scg %>% group_by(origin2) %>%
    summarise(mge_cnt=sum(mge_cnt), scg_cnt=sum(scg_cnt)) %>% 
    mutate(mge_per_cell=mge_cnt/scg_cnt)


options(repr.plot.width=8, repr.plot.height=4, repr.plot.res=300)
gg <- (ggplot(df_add_scg7, aes(x=phylum, y=mge_per_cell, fill=phylum))
       #+ geom_boxplot(position = position_dodge(preserve = 'single'))
       + geom_col()
       + geom_hline(data=df_add_scg8, mapping=aes(yintercept=mge_per_cell), color='red', linetype='dashed')
       + geom_hline(data=. %>% group_by(origin2) %>% filter(row_number() == 1), mapping=aes(yintercept = mge_per_cell_mean_aross_phylum), color='blue', linetype='dashed')
       + scale_fill_manual(values=phyla_colours_lines, breaks = phyla_levels, guide='none')
       #+ scale_color_manual(values=color_pal_emerge, guide='none')
       #+ scale_color_manual(values=ggpubr::get_palette('npg', 7), guide='none')
       #+ theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
       #+ facet_wrap(~origin2, nrow=1, scales='free_x')
       + facet_wrap(~origin2, nrow=1)
       + labs(x='', y='MGE recombinase number per genome')
       #+ ggpubr::stat_compare_means(mapping=aes(group=Habitat), label='p.signif', hide.ns=F)
       + theme_classic()
       + coord_flip()
)


###
# fig.S7.D
###

gg <- (ggplot(df_add_scg7 %>% filter(origin2 != 'ambiguous'), aes(x=phylum, y=mge_per_cell, fill=phylum))
       #+ geom_boxplot(position = position_dodge(preserve = 'single'))
       + geom_col()
       + geom_hline(data=df_add_scg8 %>% filter(origin2 != 'ambiguous'), mapping=aes(yintercept=mge_per_cell), color='red', linetype='dashed')
       + geom_hline(data=. %>% filter(origin2 != 'ambiguous') %>% group_by(origin2) %>% filter(row_number() == 1), mapping=aes(yintercept = mge_per_cell_mean_aross_phylum), color='blue', linetype='dashed')
       + scale_fill_manual(values=phyla_colours_lines, breaks = phyla_levels, guide='none')
       #+ theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
       #+ facet_wrap(~origin2, nrow=1, scales='free_x')
       + facet_wrap(~origin2, nrow=1, scales='free_x')
       + labs(x='', y='MGE recombinases per genome (contig)')
       #+ ggpubr::stat_compare_means(mapping=aes(group=Habitat), label='p.signif', hide.ns=F)
       + theme_classic()
       + coord_flip()
)

fig.s5d <- gg



###
# fig.S7.A.v2
###

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
  ) #%>% rev()


phyla_colours <- ggsci::pal_d3(palette = 'category20', alpha=1)(20)[-4] %>% #remove red, not good w/ green
  head(17) %>% c("grey50")
phyla_colours[7] <- 'grey20' # change grey to black; was postion 8 but becomes 7 after removing red (4)
phyla_colours_lines <- phyla_colours

mge_levels <- c("IS_Tn", "Phage", "CE", "Integron", "ambiguous")
mge_colours <- ggsci::pal_npg()(4) %>% c('grey50')

mag_mge_plot_data <- df_add_scg7 %>%
    ungroup() %>%
    mutate(
        phylum = as.character(phylum),
        phylum = ifelse(phylum %in% phyla_levels, phylum, "Other"),
        origin2 = factor(origin2, levels = mge_levels)
    ) %>%
    group_by(phylum) %>%
    mutate(n_total = sum(mge_cnt), per_genome_total = sum(mge_per_cell), n=mge_cnt, per_genome=mge_per_cell, n_genomes = scg_cnt, n_genomes_total=sum(scg_cnt)) %>%
    filter(phylum != 'Other') %>%
    ungroup()


### y as absolute number

mag_mge_plot_data2 <- mag_mge_plot_data %>%
  arrange(desc(per_genome_total)) %>%
  mutate(phylum = as.character(phylum))


phyla_levels2 <- unique(mag_mge_plot_data2$phylum)
mag_mge_plot_data2$phylum <- factor(mag_mge_plot_data2$phylum, levels = phyla_levels2)

mag_mge_plot <-
  mag_mge_plot_data2 %>%
  ggplot(aes(phylum, n)) +
  geom_col(aes(fill = origin2)) +
  geom_text(aes(y = n_total, label = per_genome_total %>% round(0)), hjust = 1.2, data = . %>% select(phylum, n_total, per_genome_total) %>% distinct()) +
  coord_flip() +
  scale_fill_manual("", breaks = mge_levels, values = mge_colours) +
  scale_y_reverse(labels = c('0', '100K', '200K', '300K', '400K'), breaks = c(0, 1e+5, 2e+5, 3e+5, 4e+5), limits = c(0, 4.1e+5) %>% rev) +
  scale_x_discrete(position = 'top') +
  guides(fill=guide_legend(nrow = 3, byrow = T)) +
  theme_classic() +
  theme(legend.position = 'bottom') +
  #ylim(0, 5.5e+5) + # NOTE: this is redudant w/ scale_y_continuous(limits=limits)
  labs(x="", y="MGE recombinase number")

fig.s5.v2.a <- mag_mge_plot


###
# fig.S7.B.v2
###

tile_layers <- list(
    geom_tile(),
    #xlab("MGE type"),
    #ylab("Phylum"),
    labs(x=element_blank(), y=element_blank()),
    theme_classic()
    #cowplot::theme_cowplot()
)

options(repr.plot.width=1.5, repr.plot.height=3.5, repr.plot.res=300)
alpha_per_genomes_plot <- mag_mge_plot_data2 %>%
    select(origin2, phylum, per_genome) %>%
    complete(origin2, phylum, fill = list(per_genome=0)) %>%
    ggplot(aes(origin2, phylum, fill = per_genome)) +
    tile_layers +
    #scale_fill_distiller("MGE recombinase \nper genome", palette = "YlGn", direction = 1) +
    scale_fill_distiller('MGE recombinase\nper genome', palette = "YlGn", direction = 1) +
    #guides(fill = guide_legend(title.position = 'bottom')) +
    guides(x=guide_axis(angle=45), fill=guide_colorbar(title.position = 'top', frame.colour = 'black', ticks.colour = 'black')) +
    theme(
      legend.position = "bottom",
      legend.justification = "centre",
      axis.title.y = element_blank(),
      axis.text.y = element_blank(),
      )

###
# fig.S7.C.v2
###

tmp_df <- mag_mge_plot_data %>%
    select(phylum, n_total, n_genomes_total) %>%
    distinct()

average_mge_per_genome <- (tmp_df %>% pull(n_total) %>% sum) / (tmp_df %>% pull(n_genomes_total) %>% sum)

genomes_comp_plot <- tmp_df %>%
  ggplot(aes(n_genomes_total, n_total, colour = phylum, label = phylum)) +
  geom_point() +
  geom_abline(slope = average_mge_per_genome, intercept = 0, linetype = "dashed") +
  ggrepel::geom_text_repel(size = 7/.pt, data = . %>% filter(phylum != "other"), force = 3, force_pull = 0.5, min.segment.length = 0.2) +
  scale_color_manual(values = phyla_colours_lines, breaks = phyla_levels, guide = "none") +
  scale_x_continuous(labels = scales::unit_format(unit='K', scale = 1e-3)) +
  scale_y_continuous(labels = scales::unit_format(unit='K', scale = 1e-3)) +
  xlab("Genome number") +
  ylab("MGE recombinase number") +
  theme_classic()


### log log

genomes_comp_plot_loglog <- tmp_df %>%
  ggplot(aes(n_genomes_total, n_total, colour = phylum, label = phylum)) +
  geom_point() +
  geom_abline(slope = average_mge_per_genome, intercept = 0, linetype = "dashed") +
  ggrepel::geom_text_repel(size = 7/.pt, data = . %>% filter(phylum != "other"), force = 3, force_pull = 0.5, min.segment.length = 0.2) +
  scale_color_manual(values = phyla_colours_lines, breaks = phyla_levels, guide = "none") +
  scale_x_log10() +
  scale_y_log10() +
  xlab("Genome number") +
  ylab("MGE recombinase number") +
  theme_classic()


### inset

genomes_comp_plot_inset <- tmp_df %>%
  ggplot(aes(n_genomes_total, n_total, colour = phylum, label = phylum)) +
  geom_abline(slope = average_mge_per_genome, intercept = 0, linetype = "dashed") +
  geom_point() +
  ggrepel::geom_text_repel(size = 7/.pt, data = . %>% filter(phylum != "other"), force = 3, force_pull = 0.5, min.segment.length = 0.2) +
  scale_color_manual(values = phyla_colours_lines, breaks = phyla_levels, guide = "none") +
  scale_x_continuous(labels = scales::unit_format(unit = 'K', scale = 1e-3), limits = c(0, 3000)) +
  scale_y_continuous(labels = scales::unit_format(unit = 'K', scale = 1e-3), limits = c(0, 16000)) +
  xlab("") +
  ylab("") +
  theme_classic()



### Add MAG host lineage data to compare with contig data here

### dirctories
v2_mags_directory <- here::here("som-data", "fig-data", "emerge_mags_v2")
recombinase_directory <- here::here("som-data")
v3_contig_tracking_directory <- here::here("som-data", "fig-data", "contig_tracking_v3")
recombinase_clustering_directory <- here::here("som-data", "fig-data", "recombinase_clustering_v1")

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
    d <- read_tsv(here(recombinase_directory, "mge_recombinase.tsv"), show_col_types = FALSE)
    return(d)
}

recombinase_contig_info <- read_recombinase_contig_info()

### filtered contigs with length minimal length cutoff
list_filtered_contigs <- function(recombinase_contig_info = read_recombinase_contig_info()) {
    d <- recombinase_contig_info %>%
        filter(contig_length >= contig_length_cutoff) %>%
        select(contig)

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
    mutate(seq_model_simple = if_else(stringr::str_detect(seq_model, 'NovaSeq'), 'JGI', 'Cronin'))

mag_metadata_df_ori <- read_tsv(mag_metadata_f, col_types = cols()) %>% rename(mag_name=MAG) 

mag_metadata_df <- mag_metadata_df_ori %>%
    filter(stringr::str_detect(SampleID__, 'MainAutochamber')) %>%
    filter(folder %in% c('JGI', 'Cronin_v1', 'Cronin_v2')) %>%
    mutate(seq_model_simple = if_else(folder %in% c('Cronin_v1', 'Cronin_v2'), 'Cronin', 'JGI')) %>%
    left_join(sample_metadata_df, by = c('SampleID__', 'seq_model_simple')) %>%
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

n_contig_binned <- df_contig_tracking_filt %>% nrow

mge_to_mags_checkm2_filt <- mge_to_mags_checkm2 %>%
  inner_join(df_contig_tracking_filt)

n_rec_binned <- mge_to_mags_checkm2_filt %>% nrow
# total recombinase encoding contigs >= 3kbp
n_contig_total <- recombinase_contig_info %>% 
    filter(contig_length>=3000) %>% 
    select(contig, contig_length) %>% 
    distinct() %>% nrow
# some stats
cat('[INFO] # of contigs >=3kb binned: ', n_contig_binned, '\n')
cat('[INFO] # of recombinase binned: ', n_rec_binned, '\n')
cat('[INFO] # of contits >=3kb: ', n_contig_total, '\n')
cat('[INFO] % of contigs >=3kb binned: ', scales::percent(n_contig_binned/n_contig_total, accuracy = .1), '\n')

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

type_proportions <- mag_contigs %>%
  group_by(phylum, origin2) %>%
  summarise(n = n()) %>%
  mutate(prop = n / sum(n)) %>%
  mutate(n_total = sum(n), prop_total =  sum(prop))

phylum_species_counts <- mag_metadata_df %>%
    separate(Classification, sep =';', into = c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species')) %>%
    group_by(phylum) %>%
    summarise(
        n_species = unique(mag_cluster) %>% length(),
        n_genomes = n()
    )

###
# fig.S8.B
###

average_mge_per_genome <- (mag_contigs %>% count() %>% pull(n)) / (phylum_species_counts %>% summarise(n_genomes = sum(n_genomes)) %>% pull(n_genomes))


genomes_comp_plot_mag <- mag_contigs %>%
  count(phylum) %>%
  left_join(phylum_species_counts) %>%
  mutate(
    phylum = str_remove(phylum, "p__"),
    phylum = ifelse(phylum %in% phyla_levels, phylum, "other"),
    phylum = factor(phylum, levels = phyla_levels)
    ) %>%
  ggplot(aes(n_genomes, n, colour = phylum, label = phylum)) +
  geom_point() +
  geom_abline(slope = average_mge_per_genome, intercept = 0, linetype = "dashed") +
  ggrepel::geom_text_repel(size = 7/.pt, data = . %>% filter(phylum != "other"), force = 3, force_pull = 0.5, min.segment.length = 0.2) +
  scale_color_manual(values = phyla_colours_lines, breaks = phyla_levels, guide = "none") +
  scale_x_continuous(labels = scales::unit_format(unit='K', scale = 1e-3)) +
  scale_y_continuous(labels = scales::unit_format(unit='K', scale = 1e-3)) +
  xlab("Genome number") +
  ylab("MGE recombinase number") +
  theme_classic()


### inset 
genomes_comp_plot_mag_inset <- mag_contigs %>%
  count(phylum) %>%
  left_join(phylum_species_counts) %>%
  mutate(
    phylum = str_remove(phylum, "p__"),
    phylum = ifelse(phylum %in% phyla_levels, phylum, "Other"),
    phylum = factor(phylum, levels = phyla_levels)
    ) %>%
  filter(phylum != 'Other') %>%
  ggplot(aes(n_genomes, n, colour = phylum, label = phylum)) +
  geom_abline(slope = average_mge_per_genome, intercept = 0, linetype = "dashed") +
  geom_point() +
  ggrepel::geom_text_repel(size = 7/.pt, data = . %>% filter(phylum != "Other"), force = 3, force_pull = 0.5, min.segment.length = 0.2) +
  scale_color_manual(values = phyla_colours_lines, breaks = phyla_levels, guide = "none") +
  scale_x_continuous(labels = scales::unit_format(unit='K', scale = 1e-3), limits = c(0, 400)) +
  scale_y_continuous(labels = scales::unit_format(unit='K', scale = 1e-3), limits = c(0, 5500)) +
  xlab("") +
  ylab("") +
  theme_classic()


### Fig.S8.C

tmp_df <- mag_contigs %>%
  group_by(origin2, phylum) %>%
  summarise(n_per_origin_per_phylum = n()) %>%
  left_join(phylum_species_counts) %>%
  mutate(mge_per_cell_global = sum(n_per_origin_per_phylum) / sum(n_genomes)) %>%
  mutate(
    phylum = str_remove(phylum, "p__"),
    phylum = ifelse(phylum %in% phyla_levels, phylum, "Other"),
    phylum = factor(phylum, levels = phyla_levels %>% rev)
    ) %>%
  filter(phylum != 'Other') %>%
  mutate(mge_per_cell = n_per_origin_per_phylum / n_genomes) %>%
  mutate(mge_per_cell_mean_across_phylum = mean(mge_per_cell)) #%>%
  #filter(origin2 != 'ambiguous')

tmp_df2 <- mag_contigs %>%
  group_by(origin2, phylum, MAG) %>%
  summarise(n_per_mag = n()) #%>%
  #filter(origin2 != 'ambiguous')


### boxplot
gg <- (ggplot(tmp_df2, aes(x=phylum, y=n_per_mag, fill=phylum))
       + geom_boxplot(position = position_dodge(preserve = 'single'))
       #+ geom_col()
       + geom_hline(data= tmp_df %>% group_by(origin2) %>% filter(row_number() == 1), mapping=aes(yintercept=mge_per_cell_global), color='red', linetype='dashed')
       + geom_hline(data= tmp_df %>% group_by(origin2) %>% filter(row_number() == 1), mapping=aes(yintercept = mge_per_cell_mean_across_phylum), color='blue', linetype='dashed')
       + scale_fill_manual(values=phyla_colours_lines, breaks = phyla_levels, guide='none')
       #+ theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
       #+ facet_wrap(~origin2, nrow=1, scales='free_x')
       + facet_wrap(~origin2, nrow=1, scales='free_x')
       + labs(x='', y='MGE recombinase number per genome')
       #+ ggpubr::stat_compare_means(mapping=aes(group=Habitat), label='p.signif', hide.ns=F)
       + theme_classic()
       + coord_flip()
)

### bar chart to be consistent with contig fig
gg <- (ggplot(tmp_df %>% filter(origin2 != 'ambiguous'), aes(x=phylum, y=mge_per_cell, fill=phylum))
       #+ geom_boxplot(position = position_dodge(preserve = 'single'))
       + geom_col()
       + geom_hline(data=. %>% filter(origin2 != 'ambiguous') %>% group_by(origin2) %>% filter(row_number() == 1), mapping=aes(yintercept=mge_per_cell_global), color='red', linetype='dashed')
       + geom_hline(data=. %>% filter(origin2 != 'ambiguous') %>% group_by(origin2) %>% filter(row_number() == 1), mapping=aes(yintercept = mge_per_cell_mean_across_phylum), color='blue', linetype='dashed')
       + scale_fill_manual(values=phyla_colours_lines, breaks = phyla_levels, guide='none')
       #+ theme(axis.text.x=element_text(angle=45, hjust=1, vjust=1))
       #+ facet_wrap(~origin2, nrow=1, scales='free_x')
       + facet_wrap(~origin2, nrow=1, scales='free_x')
       + labs(x='', y='MGE recombinases per genome (MAG)')
       #+ ggpubr::stat_compare_means(mapping=aes(group=Habitat), label='p.signif', hide.ns=F)
       + theme_classic()
       + coord_flip()
)

fig.s8.c <- wrap_elements(full=gg)


### fig.S7

options(repr.plot.width=7.2, repr.plot.height=9, repr.plot.res=300)
upper <- (mag_mge_plot & theme(legend.key.size = unit(8, 'pt'))) + (alpha_per_genomes_plot& theme(legend.key.size = unit(10, 'pt'))) + guide_area() + plot_layout(widths = c(1.8,0.8,1.4), guides = 'collect') + plot_annotation(tag_levels = 'A')
upper <- wrap_elements(full=upper)
mid <- (fig.s5c & theme(legend.key.size = unit(8, 'pt')) & guides(fill=guide_legend(ncol=2)) & labs(tag='C'))
mid <- wrap_elements(full=mid)
lower <- (fig.s5d & labs(tag='D'))
lower <- wrap_elements(full=lower)
p <- upper / mid / lower

figdir <- here::here('fig.outdir')
dir.create(figdir)

figfile <- here::here(figdir, 'fig.s07.host_lineage.contig.pdf')
ggsave(figfile, width = 7.2, height = 9, dpi = 300, device = 'pdf')

figfile <- here::here(figdir, 'fig.s07.host_lineage.contig.png')
ggsave(figfile, width = 7.2, height = 9, dpi = 300, device = 'png')


###
# fig.S08
###
options(repr.plot.width=7.2, repr.plot.height=9, repr.plot.res=300)
layout <- '
AA
BC
DE
'
p <- (fig.s8.c & labs(tag='A')) + (genomes_comp_plot & labs(y='MGE recombinase number (contig)', tag='B')) + genomes_comp_plot_inset + 
    (genomes_comp_plot_mag & labs(y='MGE recombinase number (MAG)', tag='C')) + genomes_comp_plot_mag_inset + plot_layout(design = layout, widths=c(2, 4.2))

figdir <- here::here('fig.outdir')
dir.create(figdir)

figfile <- here::here(figdir, 'fig.s08.host_lineage.contig_vs_mag.v2.pdf')
ggsave(figfile, width = 7.2, height = 9, dpi = 300, device = 'pdf')

figfile <- here::here(figdir, 'fig.s08.host_lineage.contig_vs_mag.v2.png')
ggsave(figfile, width = 7.2, height = 9, dpi = 300, device = 'png')



