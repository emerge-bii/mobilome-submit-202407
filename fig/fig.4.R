source(here::here('setup.R'))


setwd('som-data/fig-data')
#############
# Fig.4.A
############

allinfo_tab <- '../metat.gene_cnt.tsv'
sanity_check_tab <- '../sample.metat.qc.tsv'
rp_gene_lst_f <- 'contig_taxa/rp/rp.ko.bac_arc_shared.list'
rp_gene_vec <- readLines(rp_gene_lst_f)
mge_levels <- c("IS_Tn", "Phage", "CE", "Integron", "ambiguous")
mge_colours <- ggsci::pal_npg()(4) %>% c('grey50')
names(mge_colours) <- mge_levels
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
                           "80+")))))))))
df_sanity_check <- read_tsv(file=sanity_check_tab, col_types = cols()) %>%
    #dplyr::filter(read_forward_mapped >=0.8 | read_reverse_mapped >=0.8)
    dplyr::filter(read_reverse_mapped >=0.8)
sample_vec <- df_sanity_check$Sample
df_pre <- df %>% dplyr::filter(sample %in% sample_vec) # could be used for rec sep analysis later
df <- df_pre %>%
    dplyr::mutate(ko=if_else(stringr::str_starts(ko, 'rec__'), 'recombinase', ko)) %>% # this line for combining sep rec
    dplyr::mutate(anno=if_else(ko=='recombinase', 'recombinase', anno)) %>% 
    dplyr::mutate(anno=if_else(is.na(anno), ko, anno)) %>%
    dplyr::filter(cnt_metag>=20) %>%
    # convert gene names for better plots
    dplyr::mutate(anno=if_else(anno=='recombinase', 'MGE\nrecombinase',
                              if_else(stringr::str_detect(anno, 'glnA'), 'glnA',
                                     if_else(stringr::str_detect(anno, 'hupB'), 'hupB',
                                            if_else(stringr::str_detect(anno, 'rpoE'), 'rpoE',
                                                   if_else(stringr::str_detect(anno, 'cspA'), 'cspA', anno))))))


df_ribo_sum_per_sample <- df %>% 
    #filter(stringr::str_detect(anno, regex('ribosomal protein', ignore_case = T))) %>%
    filter(ko %in% rp_gene_vec) %>%
    dplyr::select(one_of(c('cnt_metat', 'cnt_metag', 'ko', 'anno', 'sample', 'Year', 'Habitat', 'DepthAvg', 'DepthLumping'))) %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(cnt_metat=mean(cnt_metat), cnt_metag=mean(cnt_metag), ko='RP', anno='RP', across(Year:DepthLumping, ~head(.x, n=1))) %>%
    dplyr::mutate(ratio=cnt_metat/cnt_metag, rp_ratio=cnt_metat/cnt_metag, rp_ave_cnt_metat=cnt_metat)
df_nonribo <- df %>% 
    #filter(!stringr::str_detect(anno, regex('ribosomal protein', ignore_case = T))) %>%
    dplyr::filter(!ko %in% rp_gene_vec) %>%
    dplyr::select(c('cnt_metat', 'cnt_metag', 'ko', 'anno', 'sample', 'Year', 'Habitat', 'DepthAvg', 'DepthLumping', 'ratio'))
df_ribo_sum_per_sample2 <- df_ribo_sum_per_sample %>% dplyr::select(-c('rp_ratio', 'rp_ave_cnt_metat'))
df2 <- rbind(df_ribo_sum_per_sample2, df_nonribo)

# overall metaT gene count
df_tmp2 <- df2 %>% dplyr::group_by(ko, anno) %>% 
    dplyr::summarise(cnt_metat=sum(cnt_metat)) %>% ungroup() %>% 
    dplyr::arrange(desc(cnt_metat)) %>%
    #dplyr::filter(!is.na(anno)) %>%
    #dplyr::filter(!stringr::str_detect(anno, regex('ribosomal protein', ignore_case = T))) %>%
    dplyr::filter(!ko %in% rp_gene_vec) %>%
    #dplyr::slice_head(n=5) %>%
    dplyr::filter(row_number()<=5 | anno=='RP')

col_purple <- RColorBrewer::brewer.pal(n = 3, name='BuPu') %>% rev %>% head(n=1)
df_tmp2 <- df_tmp2 %>% 
    #mutate(fill = if_else(anno == 'MGE\nrecombinase', 'grey50', 'NA'))
    mutate(fill = if_else(anno == 'MGE\nrecombinase', col_purple, 'NA'))
color_pal <- df_tmp2$fill
names(color_pal) <- df_tmp2$anno
gg_a1 <- (ggplot(data=df_tmp2, aes(x=factor(anno, levels=anno %>% rev), y=cnt_metat, fill=anno)) + geom_col(color='grey20') 
          + theme_classic() 
          + scale_x_discrete(labels=c(df_tmp2$anno[1:5], 'RP'))
          + scale_fill_manual(values = color_pal)
          + theme(legend.position='none')
          + labs(x="", y="Number of genes active")
          + theme(axis.text.x = element_blank())
         )

df_ribo_sum_per_sample3 <- df_ribo_sum_per_sample %>% 
    dplyr::select(c('sample', 'rp_ratio', 'rp_ave_cnt_metat'))
df2 <- df2 %>% inner_join(df_ribo_sum_per_sample3, by='sample') %>%
    mutate(norm_ratio=ratio/rp_ratio, norm_cnt_metat=cnt_metat/rp_ave_cnt_metat)
df_rec <- df2 %>% dplyr::mutate(anno2=if_else(ko %in% rp_gene_vec, 'RP', 
                                             if_else(anno=='recombinase', 'recombinase', anno))) %>%
    dplyr::filter(anno2=='recombinase'|anno2=='RP' |ko %in% df_tmp2$ko | stringr::str_detect(anno2, regex('lexA|recA', ignore_case = F))) %>%
    dplyr::mutate(anno2=if_else(stringr::str_detect(anno2, 'lexA'), 'lexA', 
                         if_else(stringr::str_detect(anno2, 'recA'), 'recA', anno2)))

df_rec_pivot_long <- df_rec %>%
    tidyr::pivot_longer(cols=c('ratio', 'norm_ratio', 'norm_cnt_metat'), names_to = 'name', values_to='value')

### fig6A
gg_a2 <- (ggplot(data=df_rec_pivot_long %>% filter(name %in% c('ratio') & (!stringr::str_detect(anno2, regex('lexA|recA', ignore_case=F)))) %>% mutate(anno2=factor(anno2, levels=c(df_tmp2$anno %>% rev))), 
                 aes(x=anno2, y=value, fill=anno2))
       + geom_boxplot(outlier.size = 0.8, outlier.shape = 21, outlier.alpha = 0.6)
       + ggpubr::stat_compare_means(label='p.signif', hide.ns=T, ref.group='RP', vjust=1)
       + scale_y_continuous(trans='log10', labels = scales::percent, limits = c(0.001, 1.1))
       + scale_fill_manual(values = color_pal, guide='none')
       + theme_classic()
       + labs(x='', y='Percentage of genes active')
       + annotation_logticks(sides='l')
       + theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1))
       )
fig6a1 <- gg_a1
fig6a2 <- gg_a2
fig6a <- gg_a1 / gg_a2 #+ plot_layout(axes='collect')
fig6a <- wrap_elements(full=fig6a)





##################
# Fig.4.B
##########

###
# fig6b_left
###

gg <- (ggplot(data=df_rec_pivot_long %>% filter(name %in% c('norm_ratio') & (anno2 %in% c('MGE\nrecombinase'))), aes(x=Habitat, y=value, color=Habitat)) 
       + geom_boxplot(outlier.size = 0.8, outlier.shape = 21, outlier.alpha = 0.8) + facet_wrap(~anno2, nrow = 1, scales = 'free_y')
       + scale_color_manual(values=color_pal_emerge, guide='none')
       + ggpubr::stat_compare_means(label='p.signif', hide.ns=T, ref.group='Palsa', vjust=2.5)
       + guides(x = guide_axis(angle = 45))
       + theme_classic()
       + theme(strip.background = element_rect(fill=alpha(col_purple, 0.4)))
       + labs(x='', y='Normalized active ratio')
       )
fig6b_left <- gg


###
#  fig6c_right
###
df <- df_pre %>%
    dplyr::mutate(anno=if_else(stringr::str_starts(string = ko, 'rec__'), ko, anno)) %>%
    dplyr::filter(cnt_metag>=20)
df_ribo <- df %>% filter(stringr::str_detect(anno, regex('ribosomal protein', ignore_case = T))) %>%
    group_by(sample) %>%
    summarise(rp_ratio=sum(cnt_metat)/sum(cnt_metag), rp_ave_cnt_metat=sum(cnt_metat)/n())
df2 <- df %>% inner_join(df_ribo, by='sample') %>%
    mutate(norm_ratio=ratio/rp_ratio, norm_cnt_metat=cnt_metat/rp_ave_cnt_metat)
df_rec <- df2 %>% dplyr::filter(stringr::str_starts(string = ko, 'rec__')) %>%
    dplyr::mutate(ko=stringr::str_remove(string = ko, pattern = '^rec__'))
df_rec_pivot_long <- tidyr::pivot_longer(df_rec, 
                                         #cols=c('rank_metat', 'rank_ratio', 'cnt_metat', 'cnt_metag', 'ratio'), 
                                         cols=c('norm_cnt_metat', 'norm_ratio', 'ratio'), 
                                         names_to = 'name', values_to='value')
comparisons_lst <- list(c('Palsa', 'Fen'))
options(repr.plot.width=2, repr.plot.height=2.5, repr.plot.res=300)
gg <- (ggplot(data=df_rec_pivot_long %>% dplyr::filter(name=='norm_ratio' & ko %in% c('IS_Tn')) %>% dplyr::mutate(ko=if_else(ko=='IS_Tn', 'IS_Tn\n', ko)), 
              aes(x=Habitat, y=value, color=Habitat)) 
       #+ geom_point() 
       + geom_boxplot(outlier.size = 0.8, outlier.shape = 21, outlier.alpha = 0.8) + facet_wrap(~ko, nrow = 1)
       + scale_color_manual(values=color_pal_emerge, guide='none') 
       #+ ggpubr::stat_compare_means(mapping=aes(group=Habitat), label='p.signif', hide.ns=F) 
       #+ ggpubr::stat_compare_means(label='p.signif', hide.ns=T, ref.group='Palsa', vjust = 2.5)
       + ggpubr::stat_compare_means(label='p.signif', hide.ns=T, comparisons = comparisons_lst, vjust=0.5)
       + theme_classic()
       + theme(strip.background = element_rect(fill=alpha(mge_colours['IS_Tn'], 0.4)))
       + guides(x=guide_axis(angle = 45))
       + labs(x='', y=element_blank())
       )
fig6b_mid <- gg

gg <- (ggplot(data=df_rec_pivot_long %>% dplyr::filter(name=='norm_ratio' & ko %in% c('Phage')) %>% dplyr::mutate(ko=if_else(ko=='Phage', 'Phage\n', ko)), 
              aes(x=Habitat, y=value, color=Habitat)) 
       #+ geom_point() 
       + geom_boxplot(outlier.size = 0.8, outlier.shape = 21, outlier.alpha = 0.8) + facet_wrap(~ko, nrow = 1)
       + scale_color_manual(values=color_pal_emerge, guide='none') 
       #+ ggpubr::stat_compare_means(mapping=aes(group=Habitat), label='p.signif', hide.ns=F) 
       + ggpubr::stat_compare_means(label='p.signif', hide.ns=T, ref.group='Palsa', vjust = 0.5)
       + theme_classic()
       + theme(strip.background = element_rect(fill=alpha(mge_colours['Phage'], 0.4)))
       + guides(x=guide_axis(angle = 45))
       + labs(x='', y=element_blank())
       )

fig6b_right <- gg
fig6b <- fig6b_left + fig6b_mid + fig6b_right
fig6b <- wrap_elements(full=fig6b)


###################
# Fig.4.C
#############


### PR include all RP gene shared between bac and arc
gene_tab_vec <- c('contig_taxa/rp/all_info.all_gene.fix_domain.tsv')
active_gene_list_f <- 'metat/add-hiseq-shared/all.sample.cov_0d9.active.list'
rec_tab_f <- '../mge_recombinase.tsv'
sanity_check_tab <- '../sample.metat.qc.tsv'
non_rp_gene_vec <- c('IS_Tn', 'Phage')
rp_gene_lst_f <- 'contig_taxa/rp/rp.ko.bac_arc_shared.list'

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
       #+ theme(strip.background = element_rect(fill=alpha(col_purple, 1)))
        )

fig6c <- gg
fig6c <- wrap_elements(full=fig6c)


#########
# Fig.4.D
###########

###
# combined pw file, large (~700MB)
# remove ambigous recombinases
# define unstable OTU as w/ at least one low_ani_perc < 1
# w/o OTU count label at top
###

pw_iden_f <- '../genomic_neighborhood.pairwise_ani.tsv'
rec_allinfo_f <- '../mge_recombinase.tsv'
mge_levels <- c("IS_Tn", "Phage", "CE", "Integron", "ambiguous")
mge_colours <- ggsci::pal_npg()(4) %>% c('grey50')
names(mge_colours) <- mge_levels
# variable for neighborhood types
stable <- 'stable'     # vertical transfer
unstable <- 'unstable' # horizontal
df_pw <- read.table(pw_iden_f, sep='\t', header=T)
colnames(df_pw) <- c('seqname1', 'seqname2', 'low_ani', 'high_ani', 'OTU')
df_rec <- read_tsv(rec_allinfo_f, col_types=cols())
df_info <- df_rec %>%
    dplyr::rename(mem=recombinase) %>%
    dplyr::filter(!is.na(OTU)) %>%
    dplyr::select(all_of(c('mem', 'domain', 'phylum', 'Year', 'Habitat', 'DepthAvg', 'origin'))) %>%
    dplyr::mutate(Habitat=factor(Habitat, levels=c('Palsa', 'Bog', 'Fen'))) %>%
    dplyr::mutate(DepthLumping =  ifelse(DepthAvg >= 0 & DepthAvg < 10, "0-9", 
                           ifelse(DepthAvg >= 10 & DepthAvg < 20, "10-19", 
                           ifelse(DepthAvg >= 20 & DepthAvg < 30, "20-29",
                           ifelse(DepthAvg >= 30 & DepthAvg < 40, "30-39",
                           ifelse(DepthAvg >= 40 & DepthAvg < 50, "40-49",
                           ifelse(DepthAvg >= 50 & DepthAvg < 60, "50-59",
                           ifelse(DepthAvg >= 60 & DepthAvg < 70, "60-69",
                           ifelse(DepthAvg >= 70 & DepthAvg < 80, "70-79",
                           "80+")))))))))

# remove non-MGE recombinases
df_pw <- df_pw %>% dplyr::filter(seqname1 %in% df_info$mem & seqname2 %in% df_info$mem)
df_merged <- df_pw %>% dplyr::left_join(df_info, by=c('seqname1'='mem')) %>%
    dplyr::left_join(df_info, by=c('seqname2'='mem'), suffix=c('1', '2')) %>%
    dplyr::mutate(Year_diff=factor(abs(Year1-Year2))) %>%
    dplyr::filter(origin1==origin2) %>%
    #dplyr::filter(Habitat1==Habitat2 & Habitat1=="Bog") %>%
    dplyr::mutate(origin=origin1) %>%
    dplyr::filter(origin %in% c('CE', 'Integron', 'IS_Tn', 'Phage'))

df_sum <- df_merged %>% dplyr::group_by(OTU) %>% dplyr::summarise(low_ani_cnt=sum(low_ani<1), total=n()) %>% 
    dplyr::mutate(low_ani_perc=low_ani_cnt/total*100)
### change the number here to set minimal OTU size;
### "total" is the total pairs; 1 pair -> OTU size of 2; 10 pairs -> OTU size of 5;
df_sum <- df_sum %>% dplyr::filter(total>=1)
unstable_otu_vec <- df_sum %>% dplyr::filter(low_ani_cnt>=1) %>% dplyr::pull(OTU)
df_merged <- df_merged %>% dplyr::filter(OTU %in% df_sum$OTU)

###
# all data per MGE type; fig7a
###
df_all <- df_merged %>%
    dplyr::group_by(origin, OTU) %>%
    dplyr::summarise(low_ani_cnt=sum(low_ani<1), total=n()) %>%
    dplyr::mutate(low_ani_perc=low_ani_cnt/total*100) %>%
    #dplyr::mutate(transfer_type=if_else(OTU %in% unstable_otu_vec, unstable, stable)) %>%
    dplyr::mutate(transfer_type=if_else(low_ani_cnt>=1, unstable, stable)) %>%
    dplyr::mutate(transfer_type=factor(transfer_type, levels=c(unstable, stable)))
df_cnt <- df_all %>% group_by(origin) %>% summarise(n=n())
df_all_summary <- df_all %>% 
    group_by(origin, transfer_type) %>% 
    summarise(transfer_type_cnt=n()) %>% 
    mutate(transfer_type_prop=transfer_type_cnt/sum(transfer_type_cnt)) %>%
    filter(transfer_type==unstable) %>%
    mutate(origin = factor(origin, levels=mge_levels[mge_levels!='ambiguous'])) %>%
    ungroup()

gg <- (ggplot(data=df_all_summary, aes(x=origin, y=transfer_type_prop, fill=origin))
       + geom_col()
       #+ geom_text(stat='count', aes(label=..count..))
       #+ geom_text(data=df_cnt, mapping=aes(x=origin, y=max(df_all_summary %>% pull(transfer_type_prop)) * 1.05, label=n, fill=unstable), size=1.5)
       #+ facet_wrap(~origin, ncol=1, scale='free_y')
       + theme_classic()
       + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
       + scale_fill_manual(values = mge_colours, name=element_blank(), guide='none')
       + scale_y_continuous(labels = scales::percent_format(accuracy = 0.1L))
       + labs(x='', y='Percentage of unstable OTUs')
       )
fig7a <- gg


##########################
# Fig.4.D.partial_carriage
#############################

### partial carriage

tab_f <- 'partial_carriage/all_recombinase.coord.min_up_down_20000.rec_vs_contig_depth_cov.add_allinfo.rm_other_rec.dedup.add_cluster.tsv'
stable <- 'N'
unstable <- 'Y'
cohen_d_cutoff <- 1
mge_levels <- c("IS_Tn", "Phage", "CE", "Integron", "ambiguous")
mge_colours <- ggsci::pal_npg()(4) %>% c('grey50')
names(mge_colours) <- mge_levels

df <- read_tsv(tab_f, col_types = cols()) %>%
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
                           "80-89")))))))))

###
# all data per MGE type; fig8a
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
    filter(partial_carriage==unstable) %>%
    mutate(origin = factor(origin, levels=mge_levels[mge_levels!='ambiguous'])) %>%
    ungroup()

df_cnt <- df_all %>% dplyr::group_by(origin) %>% dplyr::summarise(n=n())
gg <- (ggplot(data=df_all_summary, aes(x=origin, y=carriage_type_prop, fill=origin))
       + geom_col()
       #+ geom_text(data=df_cnt, mapping=aes(x=origin, y=max(df_all_summary %>% pull(carriage_type_prop)) * 1.05, label=n, fill=unstable), size=1.5)
       + theme_classic()
       + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
       + scale_fill_manual(values = mge_colours, name=element_blank(), guide='none')
       + scale_y_continuous(labels=scales::percent_format(accuracy = 0.1L))
       + labs(x='', y="Percentage of OTUs w/ partial carriage by host")
       )
fig8a <- gg

################
# Fig.4.D.metat
#################
### metaT on recomibnase OTU stability or activity
### cov_0d9
rec_tab_f <- '../mge_recombinase.tsv'
active_gene_list_f <- 'metat/add-hiseq-shared/all.sample.cov_0d9.active.list'
sanity_check_tab <- '../sample.metat.qc.tsv'

mge_levels <- c("IS_Tn", "Phage", "CE", "Integron", "ambiguous")
mge_colours <- ggsci::pal_npg()(4) %>% c('grey50')
names(mge_colours) <- mge_levels
stable <- 'N'
unstable <- 'Y'
df_sanity_check <- read_tsv(file=sanity_check_tab, col_types = cols()) %>%
    dplyr::filter(read_reverse_mapped >=0.8)
sample_vec <- df_sanity_check$Sample

df <- read_tsv(rec_tab_f, col_types=cols()) %>%
    dplyr::rename(mem = recombinase) %>%
    dplyr::filter(!is.na(OTU)) %>%
    dplyr::filter(Sample %in% sample_vec)
active_gene_vec <- readLines(active_gene_list_f)
active_rec_vec <- intersect(active_gene_vec, df$mem)
# partial is meant to be "active"; used here for re-using the code
df <- df %>% 
    dplyr::filter(origin != 'Cell') %>%
    dplyr::mutate(origin = if_else(stringr::str_detect(origin, ';'), 'ambiguous', origin)) %>%
    dplyr::mutate(Habitat = factor(Habitat, levels = c('Palsa', 'Bog', 'Fen')), Year=as.factor(Year)) %>%
    dplyr::mutate(partial=if_else(mem %in% active_rec_vec, 1, 0)) %>%
    dplyr::filter(origin != 'ambiguous') %>%
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
    filter(partial_carriage==unstable) %>%
    mutate(origin = factor(origin, levels = mge_levels)) %>%
    ungroup()

df_cnt <- df_all %>% dplyr::group_by(origin) %>% dplyr::summarise(n=n())
gg <- (ggplot(data=df_all_summary, aes(x=origin, y=carriage_type_prop, fill=origin))
       + geom_col()
       #+ geom_text(data=df_cnt, mapping=aes(x=origin, y=max(df_all_summary %>% pull(carriage_type_prop)) * 1.05, label=n, fill=unstable), size=1.5)
       + theme_classic()
       + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
       + scale_fill_manual(values = mge_colours, name=element_blank(), guide='none')
       + scale_y_continuous(labels=scales::percent_format(accuracy = 0.1L))
       + labs(x='', y="Percentage of active OTUs")
       )
fig9a <- gg


layout <- '
###
ABC
'
fig6_extra <- list(
    (fig9a & labs(y='Transcription') & coord_flip() & theme(axis.text.y = element_text(margin = margin('l', 0)))),
    (fig7a & labs(y='Genomic context change') & coord_flip() & theme(axis.text.y = element_blank())),
    (fig8a & labs(y='Partial carriage by population') & coord_flip() & theme(axis.text.y = element_blank()))
) %>% wrap_plots(design=layout, widths = c(1, 1, 1)) + plot_layout(tag_level = 'new')
fig6_extra <- wrap_elements(full=fig6_extra)

options(repr.plot.width=7.2, repr.plot.height=7, repr.plot.res=300)
layout <- '
AB
AC
DD
'
fig6_add_extra <- fig6a + fig6b + fig6c + fig6_extra + 
    plot_layout(widths = c(2.5,4), heights = c(1.1,1.3,1), design=layout, guides = 'collect') + 
    plot_annotation(tag_levels = 'A')
p <- fig6_add_extra

figdir <- here::here('fig.outdir')
dir.create(figdir)

figfile <- here::here(figdir, 'fig.4.metat.add_extra.pdf')
ggsave(figfile, p, width = 7.2, height = 8, dpi = 300, device = 'pdf')

#figfile <- here::here(figdir, 'fig.4.metat.add_extra.png')
#ggsave(figfile, p, width = 7.2, height = 8, dpi = 300, device = 'png')
