source(here::here('setup.R'))

###
# 300 bp neighborhood cutoff
# define unstable OTU as w/ at least one low_ani_perc < 1
###

pw_iden_f <- 'som-data/genomic_neighborhood.pairwise_ani.tsv'
rec_allinfo_f <- 'som-data/mge_recombinase.tsv'

# variable for neighborhood types
stable <- 'stable'     # vertical transfer
unstable <- 'unstable' # horizontal

df_pw <- read_tsv(pw_iden_f, col_types = cols())
colnames(df_pw) <- c('seqname1', 'seqname2', 'low_ani', 'high_ani', 'OTU')

df_rec <- read_tsv(rec_allinfo_f, col_types=cols())

df_info <- df_rec %>%
    dplyr::rename(mem = recombinase) %>%
    dplyr::filter(Habitat %in% c('Palsa', 'Bog', 'Fen') & !is.na(OTU)) %>%
    dplyr::select(all_of(c('mem', 'domain', 'phylum', 'Year', 'Habitat', 'DepthAvg', 'origin'))) %>%
    dplyr::mutate(Habitat=factor(Habitat, levels=c('Palsa', 'Bog', 'Fen'))) %>%
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

df_pw <- df_pw %>% dplyr::filter(seqname1 %in% df_info$mem & seqname2 %in% df_info$mem)

df_merged <- df_pw %>% dplyr::left_join(df_info, by=c('seqname1'='mem')) %>%
    dplyr::left_join(df_info, by=c('seqname2'='mem'), suffix=c('1', '2')) %>%
    dplyr::mutate(Year_diff=factor(abs(Year1-Year2))) %>%
    dplyr::filter(origin1==origin2) %>%
    dplyr::mutate(origin=origin1) %>%
    dplyr::filter(origin %in% c('CE', 'Integron', 'IS_Tn', 'Phage'))


df_sum <- df_merged %>% dplyr::group_by(OTU) %>% dplyr::summarise(low_ani_cnt=sum(low_ani<1), total=n()) %>% 
    dplyr::mutate(low_ani_perc=low_ani_cnt/total*100)

### change the number here to set minimal OTU size;
### "total" is the total pairs; 1 pair -> OTU size of 2; 10 pairs -> OTU size of 5;

df_sum <- df_sum %>% dplyr::filter(total>=1)
unstable_otu_vec <- df_sum %>% dplyr::filter(low_ani_cnt>=1) %>% dplyr::pull(OTU)
df_merged <- df_merged %>% dplyr::filter(OTU %in% df_sum$OTU)

df_merged_300 <- df_merged

###
# all data per MGE type
###

df_all <- df_merged %>%
    dplyr::group_by(origin, OTU) %>%
    dplyr::summarise(low_ani_cnt=sum(low_ani<1), total=n()) %>%
    dplyr::mutate(low_ani_perc=low_ani_cnt/total*100) %>%
    dplyr::mutate(transfer_type=if_else(low_ani_cnt>=1, unstable, stable)) %>%
    dplyr::mutate(transfer_type=factor(transfer_type, levels=c(unstable, stable)))


df_cnt <- df_all %>% group_by(origin) %>% summarise(n=n())
df_all_summary <- df_all %>% 
    group_by(origin, transfer_type) %>% 
    summarise(transfer_type_cnt=n()) %>% 
    mutate(transfer_type_prop=transfer_type_cnt/sum(transfer_type_cnt)) %>%
    filter(transfer_type==unstable)


options(repr.plot.width=2.5, repr.plot.height=3, repr.plot.res=300)
gg <- (ggplot(data=df_all_summary, aes(x=origin, y=transfer_type_prop, fill=transfer_type))
       + geom_col()
       + theme_classic()
       + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
       + scale_fill_brewer(palette = 'Dark2', name=element_blank(), guide='none')
       + scale_y_continuous(labels = scales::percent)
       + labs(x='', y='Percentage of unstable OTUs')
       )

#fig7a <- gg

df_300 <- df_all_summary %>% mutate(cutoff = 300)


###
# across years
###

df_single_year <- df_merged %>% 
    dplyr::filter(Year1!=2010 & Year2!=2010) %>%
    dplyr::filter(Year1==Year2 & origin1==origin2) %>% dplyr::select(all_of(c('origin1', 'Year1', 'low_ani', 'high_ani', 'OTU'))) %>%
    dplyr::mutate(Year1=factor(Year1)) %>%
    dplyr::rename(origin=origin1, Year=Year1) %>%
    dplyr::group_by(OTU, Year, origin) %>%
    dplyr::summarise(low_ani_cnt=sum(low_ani<1), total=n()) %>%
    dplyr::mutate(low_ani_perc=low_ani_cnt/total*100) %>%
    #dplyr::mutate(transfer_type=if_else(OTU %in% unstable_otu_vec, unstable, stable)) %>%
    dplyr::mutate(transfer_type=if_else(low_ani_cnt>=1, unstable, stable)) %>%
    dplyr::mutate(transfer_type=factor(transfer_type, levels=c(unstable, stable)))


df_single_year_summary <- df_single_year %>% 
    group_by(origin, Year, transfer_type) %>% 
    summarise(transfer_type_cnt=n()) %>%
    mutate(transfer_type_prop=transfer_type_cnt/sum(transfer_type_cnt)) %>%
    filter(transfer_type==unstable)

df_cnt <- df_single_year %>% group_by(Year, origin) %>% summarise(n=n())
options(repr.plot.width=3, repr.plot.height=3, repr.plot.res=300)
gg <- (ggplot(data=df_single_year_summary %>% filter(origin %in% c('IS_Tn', 'Phage')), aes(x=Year, y=transfer_type_prop, fill=transfer_type))
       + geom_col(color='grey20')
       #+ geom_bar(position='fill', aes(fill=transfer_type))
       #+ geom_text(stat='count', aes(label=..count..))
       #+ geom_text(data=df_cnt %>% filter(origin %in% c('IS_Tn', 'Phage')), mapping=aes(x=Year, y=max(df_single_year_summary %>% pull(transfer_type_prop))*1.05, fill='unstable', label=n), size=1.5)
       + facet_wrap(~origin, ncol=1, scale='free_y')
       + theme_classic()
       + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
       #+ scale_fill_brewer(palette = 'Dark2', name=element_blank(), guide='none')
       + scale_fill_manual(values=c('white'), name=element_blank(), guide='none')
       + scale_y_continuous(labels = scales::percent)
       + labs(x='Year', y="Unstable OTUs (%)")
       #+ labs(x='', y="")
       )

fig7b <- gg


options(repr.plot.width=3, repr.plot.height=6, repr.plot.res=300)
gg <- (ggplot(data=df_single_year_summary, aes(x=Year, y=transfer_type_prop, fill=transfer_type))
       + geom_col()
       #+ geom_bar(position='fill', aes(fill=transfer_type))
       #+ geom_text(stat='count', aes(label=..count..))
       #+ geom_text(data=df_cnt, mapping=aes(x=Year, y=max(df_single_year_summary %>% pull(transfer_type_prop))*1.05, fill='unstable', label=n), size=1.5)
       + facet_wrap(~origin, ncol=1, scale='free_y')
       + theme_classic()
       + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
       + scale_fill_brewer(palette = 'Dark2', name=element_blank(), guide='none')
       + scale_y_continuous(labels = scales::percent)
       + labs(x='Year', y="Unstable OTUs (%)")
       )


###
# compared to 2011
###

df_year_diff <- df_merged %>% 
    dplyr::filter(Year1==2011 | Year2==2011) %>%
    dplyr::filter(Year1!=2010 & Year2!=2010) %>%
    dplyr::select(all_of(c('origin1', 'Year_diff', 'low_ani', 'high_ani', 'OTU'))) %>%
    dplyr::rename(origin=origin1) %>%
    dplyr::group_by(OTU, Year_diff, origin) %>%
    dplyr::summarise(low_ani_cnt=sum(low_ani<1), total=n()) %>%
    dplyr::mutate(low_ani_perc=low_ani_cnt/total*100) %>%
    #dplyr::mutate(transfer_type=if_else(OTU %in% unstable_otu_vec, unstable, stable)) %>%
    dplyr::mutate(transfer_type=if_else(low_ani_cnt>=1, unstable, stable)) %>%
    dplyr::mutate(transfer_type=factor(transfer_type, levels = c(unstable, stable)))


df_year_diff_summary <- df_year_diff %>% 
    group_by(origin, Year_diff, transfer_type) %>% 
    summarise(transfer_type_cnt=n()) %>%
    mutate(transfer_type_prop=transfer_type_cnt/sum(transfer_type_cnt)) %>%
    filter(transfer_type==unstable)


df_cnt <- df_year_diff %>% group_by(Year_diff, origin) %>% summarise(n=n())
options(repr.plot.width=3, repr.plot.height=3, repr.plot.res=300)
gg <- (ggplot(data=df_year_diff_summary %>% filter(origin %in% c('IS_Tn', 'Phage')), aes(x=Year_diff, y=transfer_type_prop, fill=transfer_type))
       + geom_col()
       #+ geom_bar(position='fill', aes(fill=transfer_type))
       #+ geom_text(stat='count', aes(label=..count..))
       #+ geom_text(data=df_cnt %>% filter(origin %in% c('IS_Tn', 'Phage')), mapping=aes(x=Year_diff, y=max(df_year_diff_summary %>% pull(transfer_type_prop))*1.05, fill='unstable', label=n), size=1.5)
       + facet_wrap(~origin, ncol=1, scale='free_y')
       + theme_classic()
       + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
       + scale_fill_brewer(palette = 'Dark2', name=element_blank(), guide='none')
       + scale_y_continuous(labels = scales::percent)
       + labs(x='Years after 2011', y="Unstable OTUs (%)")
       )


options(repr.plot.width=3, repr.plot.height=6, repr.plot.res=300)
gg <- (ggplot(data=df_year_diff_summary, aes(x=Year_diff, y=transfer_type_prop, fill=transfer_type))
       + geom_col()
       #+ geom_bar(position='fill', aes(fill=transfer_type))
       #+ geom_text(stat='count', aes(label=..count..))
       #+ geom_text(data=df_cnt, mapping=aes(x=Year_diff, y=max(df_year_diff_summary %>% pull(transfer_type_prop))*1.05, fill='unstable', label=n), size=1.5)
       + facet_wrap(~origin, ncol=1, scale='free_y')
       + theme_classic()
       + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
       + scale_fill_brewer(palette = 'Dark2', name=element_blank(), guide='none')
       + scale_y_continuous(labels = scales::percent)
       + labs(x='Years after 2011', y="Unstable OTUs (%)")
       )


###
# look at habitat; fig7c
###

df_single_habitat <- df_merged %>% dplyr::filter(Habitat1==Habitat2) %>% dplyr::select(all_of(c('origin1', 'Habitat1', 'low_ani', 'high_ani', 'OTU'))) %>%
    dplyr::rename(Habitat=Habitat1, origin=origin1) %>%
    dplyr::group_by(OTU, Habitat, origin) %>%
    dplyr::summarise(low_ani_cnt=sum(low_ani<1), total=n()) %>%
    dplyr::mutate(low_ani_perc=low_ani_cnt/total*100) %>%
    #dplyr::mutate(transfer_type=if_else(OTU %in% unstable_otu_vec, unstable, stable)) %>%
    dplyr::mutate(transfer_type=if_else(low_ani_cnt>=1, unstable, stable)) %>%
    dplyr::mutate(transfer_type=factor(transfer_type, levels=c(unstable, stable)))


df_single_habitat_summary <- df_single_habitat %>% 
    group_by(origin, Habitat, transfer_type) %>% 
    summarise(transfer_type_cnt=n()) %>%
    mutate(transfer_type_prop=transfer_type_cnt/sum(transfer_type_cnt)) %>%
    filter(transfer_type==unstable)

df_cnt <- df_single_habitat %>% group_by(Habitat, origin) %>% summarise(n=n())
options(repr.plot.width=1.5, repr.plot.height=3, repr.plot.res=300)
gg <- (ggplot(data=df_single_habitat_summary %>% filter(origin %in% c('IS_Tn', 'Phage')), aes(x=Habitat, y=transfer_type_prop, fill=transfer_type))
       + geom_col(color='grey20')
       #+ geom_bar(position='fill', aes(fill=transfer_type))
       #+ geom_text(stat='count', aes(label=..count..))
       #+ geom_text(data=df_cnt %>% filter(origin %in% c('IS_Tn', 'Phage')), mapping=aes(x=Habitat, y=max(df_single_habitat_summary %>% pull(transfer_type_prop))*1.05, fill='unstable', label=n), size=1.5)
       + facet_wrap(~origin, ncol=1, scale='free_y')
       + theme_classic()
       + theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))
       #+ scale_fill_brewer(palette = 'Dark2', name=element_blank(), guide='none')
       + scale_fill_manual(values = c('white'), name=element_blank(), guide='none')
       + scale_y_continuous(labels = scales::percent)
       #+ labs(x='Habitat', y="Proportion of unstable OTUs")
       + labs(x='', y="")
       )

fig7c <- gg


###
# look at phylum; fig7d
###

df_dummy <- df_merged %>% dplyr::filter(phylum1==phylum2, origin1==origin2) %>% 
    dplyr::select(all_of(c('origin1', 'domain1', 'phylum1', 'low_ani', 'high_ani', 'OTU'))) %>%
    dplyr::rename(origin=origin1, domain=domain1, phylum=phylum1) %>%
    dplyr::group_by(OTU, domain, phylum, origin) %>%
    dplyr::summarise(low_ani_cnt=sum(low_ani<1), total=n()) %>%
    dplyr::mutate(low_ani_perc=low_ani_cnt/total*100) %>%
    dplyr::mutate(transfer_type=if_else(OTU %in% unstable_otu_vec, unstable, stable)) %>%
    dplyr::mutate(transfer_type=factor(transfer_type, levels=c(unstable, stable))) %>%
    dplyr::mutate(domain=stringr::str_remove_all(domain, pattern = 'd__')) %>%
    dplyr::mutate(phylum=stringr::str_remove_all(phylum, pattern = 'p__')) %>%
    dplyr::filter(phylum != 'Other')


df_dummy_summary <- df_dummy %>% 
    group_by(origin, domain, phylum, transfer_type) %>% 
    summarise(transfer_type_cnt=n()) %>%
    mutate(transfer_type_prop=transfer_type_cnt/sum(transfer_type_cnt), transfer_type_sum=sum(transfer_type_cnt)) %>%
    dplyr::filter(transfer_type_sum>=10) %>% #otu_per_phylum_origin
    filter(transfer_type==unstable)

col_pal <- ggsci::pal_nejm()(2)
names(col_pal) <- c('Archaea', 'Bacteria')
df_dummy_summary$color <- col_pal[df_dummy_summary$domain]
df_dummy_summary <- df_dummy_summary %>%
    mutate(phylum2 = stringr::str_c('<span style = \"color:', color, '\">', phylum, '<span>', sep=' '))

phylum2_level <- df_dummy_summary %>% 
    group_by(domain, phylum2, transfer_type) %>%
    summarise(transfer_type_prop=sum(transfer_type_prop)) %>% 
    arrange(desc(transfer_type_prop)) %>% 
    pull(phylum2)

options(repr.plot.width=7.2, repr.plot.height=3, repr.plot.res=300)
gg <- (ggplot(data=df_dummy_summary %>% filter(origin %in% c('IS_Tn', 'Phage')), 
              aes(x=factor(phylum2, levels=phylum2_level), y=transfer_type_prop, fill=transfer_type))
       + geom_col(color='grey20')
       #+ geom_bar(position='fill', aes(fill=transfer_type))
       #+ geom_text(stat='count', aes(label=..count..))
       #+ geom_text(mapping=aes(x=phylum2, y=max(df_dummy_summary %>% pull(transfer_type_prop))*1.05, fill=unstable, label=transfer_type_sum), size=1.5)
       + facet_wrap(~origin, ncol=1, scale='free_y')
       + theme_classic()
       + theme(axis.text.x = ggtext::element_markdown(angle=45, hjust=1, vjust=1))
       #+ theme(axis.text.x = ggtext::element_markdown(angle=90, hjust=1, vjust=0.5))
       #+ scale_fill_brewer(palette = 'Dark2', name=element_blank(), guide='none')
       + scale_fill_manual(values = c('white'), name=element_blank(), guide='none')
       + scale_y_continuous(labels=scales::percent)
       + labs(x='', y="Unstable OTUs (%)")
       )

fig7d <- gg


###
# 600bp neighborhood
###

pw_iden_f <- 'som-data/fig-data/stability/nh_600bp/all.gene_neighbor.blastn.pw_iden_within_cluster.sep_ends_sum.tsv'
rec_allinfo_f <- 'som-data/mge_recombinase.tsv'

# variable for neighborhood types
stable <- 'stable'     # vertical transfer
unstable <- 'unstable' # horizontal

df_pw <- read.table(pw_iden_f, sep='\t', header=T)
colnames(df_pw) <- c('seqname1', 'seqname2', 'low_ani', 'high_ani', 'OTU')

df_rec <- read_tsv(rec_allinfo_f, col_types=cols())

df_info <- df_rec %>%
    dplyr::rename(mem = recombinase) %>%
    dplyr::filter(Habitat %in% c("Palsa", "Bog", "Fen") & !is.na(OTU)) %>%
    dplyr::select(all_of(c('mem', 'domain', 'phylum', 'Year', 'Habitat', 'DepthAvg', 'origin'))) %>%
    dplyr::mutate(Habitat=factor(Habitat, levels=c('Palsa', 'Bog', 'Fen'))) %>%
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


# remove non-MGE recombinases
df_pw <- df_pw %>% dplyr::filter(seqname1 %in% df_info$mem & seqname2 %in% df_info$mem)

df_merged <- df_pw %>% dplyr::left_join(df_info, by=c('seqname1'='mem')) %>%
    dplyr::left_join(df_info, by=c('seqname2'='mem'), suffix=c('1', '2')) %>%
    dplyr::mutate(Year_diff=factor(abs(Year1-Year2))) %>%
    dplyr::filter(origin1==origin2) %>%
    dplyr::mutate(origin=origin1) %>%
    dplyr::filter(origin %in% c('CE', 'Integron', 'IS_Tn', 'Phage'))


df_sum <- df_merged %>% dplyr::group_by(OTU) %>% dplyr::summarise(low_ani_cnt=sum(low_ani<1), total=n()) %>% 
    dplyr::mutate(low_ani_perc=low_ani_cnt/total*100)

### change the number here to set minimal OTU size;
### "total" is the total pairs; 1 pair -> OTU size of 2; 10 pairs -> OTU size of 5;

df_sum <- df_sum %>% dplyr::filter(total>=1)
unstable_otu_vec <- df_sum %>% dplyr::filter(low_ani_cnt>=1) %>% dplyr::pull(OTU)
df_merged <- df_merged %>% dplyr::filter(OTU %in% df_sum$OTU)

df_merged_600 <- df_merged

###
# all data per MGE type
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
    filter(transfer_type==unstable)


options(repr.plot.width=2.5, repr.plot.height=3, repr.plot.res=300)
gg <- (ggplot(data=df_all_summary, aes(x=origin, y=transfer_type_prop, fill=transfer_type))
       + geom_col()
       #+ geom_text(stat='count', aes(label=..count..))
       #+ geom_text(data=df_cnt, mapping=aes(x=origin, y=max(df_all_summary %>% pull(transfer_type_prop)) * 1.05, label=n, fill=unstable), size=1.5)
       #+ facet_wrap(~origin, ncol=1, scale='free_y')
       + theme_classic()
       + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
       + scale_fill_brewer(palette = 'Dark2', name=element_blank(), guide='none')
       + scale_y_continuous(labels = scales::percent)
       + labs(x='', y='Percentage of unstable OTUs')
       )

#fig7a <- gg

df_600 <- df_all_summary %>% mutate(cutoff=600)



###
# 1kb neighborhood cutoff
###

pw_iden_f <- 'som-data/fig-data/stability/nh_1kbp/all.gene_neighbor.blastn.pw_iden_within_cluster.sep_ends_sum.tsv'
rec_allinfo_f <- 'som-data/mge_recombinase.tsv'

# variable for neighborhood types
stable <- 'stable'     # vertical transfer
unstable <- 'unstable' # horizontal

df_pw <- read.table(pw_iden_f, sep='\t', header=T)
colnames(df_pw) <- c('seqname1', 'seqname2', 'low_ani', 'high_ani', 'OTU')

df_rec <- read_tsv(rec_allinfo_f, col_types=cols())

df_info <- df_rec %>%
    dplyr::rename(mem = recombinase) %>%
    dplyr::filter(Habitat %in% c("Palsa", "Bog", "Fen") & !is.na(OTU)) %>%
    dplyr::select(all_of(c('mem', 'domain', 'phylum', 'Year', 'Habitat', 'DepthAvg', 'origin'))) %>%
    dplyr::mutate(Habitat=factor(Habitat, levels=c('Palsa', 'Bog', 'Fen'))) %>%
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

total_otu_vec_1000 <- df_sum %>% dplyr::pull(OTU)
total_seqname_vec_1000 <- union(df_merged$seqname1, df_merged$seqname2)

###
# all data per MGE type
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
    filter(transfer_type==unstable)

options(repr.plot.width=2.5, repr.plot.height=3, repr.plot.res=300)
gg <- (ggplot(data=df_all_summary, aes(x=origin, y=transfer_type_prop, fill=transfer_type))
       + geom_col()
       #+ geom_text(stat='count', aes(label=..count..))
       #+ geom_text(data=df_cnt, mapping=aes(x=origin, y=max(df_all_summary %>% pull(transfer_type_prop)) * 1.05, label=n, fill=unstable), size=1.5)
       #+ facet_wrap(~origin, ncol=1, scale='free_y')
       + theme_classic()
       + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
       + scale_fill_brewer(palette = 'Dark2', name=element_blank(), guide='none')
       + scale_y_continuous(labels = scales::percent)
       + labs(x='', y='Percentage of unstable OTUs')
       )

#fig7a <- gg

df_1000 <- df_all_summary %>% mutate(cutoff = 1000)


###
# 600bp neighborhood; screened by total_otu_vec_1000 or total_seqname_vec_1000
###

#df_merged <- df_merged_600 %>% dplyr::filter(OTU %in% total_otu_vec_1000)
df_merged <- df_merged_600 %>% dplyr::filter(seqname1 %in% total_seqname_vec_1000 & seqname2 %in% total_seqname_vec_1000)

###
# all data per MGE type
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
    filter(transfer_type==unstable)


options(repr.plot.width=2.5, repr.plot.height=3, repr.plot.res=300)
gg <- (ggplot(data=df_all_summary, aes(x=origin, y=transfer_type_prop, fill=transfer_type))
       + geom_col()
       #+ geom_text(stat='count', aes(label=..count..))
       #+ geom_text(data=df_cnt, mapping=aes(x=origin, y=max(df_all_summary %>% pull(transfer_type_prop)) * 1.05, label=n, fill=unstable), size=1.5)
       #+ facet_wrap(~origin, ncol=1, scale='free_y')
       + theme_classic()
       + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
       + scale_fill_brewer(palette = 'Dark2', name=element_blank(), guide='none')
       + scale_y_continuous(labels = scales::percent)
       + labs(x='', y='Percentage of unstable OTUs')
       )

#fig7a <- gg

df_600_filt <- df_all_summary %>% mutate(cutoff=600)


###
# 300bp neighborhood; screened by total_otu_vec_1000 or total_seqname_vec_1000
###

#df_merged <- df_merged_300 %>% dplyr::filter(OTU %in% total_otu_vec_1000)
df_merged <- df_merged_300 %>% dplyr::filter(seqname1 %in% total_seqname_vec_1000 & seqname2 %in% total_seqname_vec_1000)

###
# all data per MGE type
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
    filter(transfer_type==unstable)


options(repr.plot.width=2.5, repr.plot.height=3, repr.plot.res=300)
gg <- (ggplot(data=df_all_summary, aes(x=origin, y=transfer_type_prop, fill=transfer_type))
       + geom_col()
       #+ geom_text(stat='count', aes(label=..count..))
       #+ geom_text(data=df_cnt, mapping=aes(x=origin, y=max(df_all_summary %>% pull(transfer_type_prop)) * 1.05, label=n, fill=unstable), size=1.5)
       #+ facet_wrap(~origin, ncol=1, scale='free_y')
       + theme_classic()
       + theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
       + scale_fill_brewer(palette = 'Dark2', name=element_blank(), guide='none')
       + scale_y_continuous(labels = scales::percent)
       + labs(x='', y='Percentage of unstable OTUs')
       )

#fig7a <- gg

df_300_filt <- df_all_summary %>% mutate(cutoff=300)


df_all_summary_merged <- rbind(df_300, df_600, df_1000)
df_all_summary_filt_merged <- merged <- rbind(df_300_filt, df_600_filt, df_1000)

options(repr.plot.width = 3, repr.plot.height = 4, repr.plot.res=300)
gg_sensitivity <- ggplot(data=df_all_summary_merged, aes(x=origin, y=transfer_type_prop)) +
    geom_col(fill='white', color='grey20') + facet_wrap(~cutoff, ncol = 1) +
    scale_y_continuous(labels=scales::percent) +
    labs(x='', y="Unstable OTU (%)") +
    theme_classic() +
    guides(x = guide_axis(angle = 45))

options(repr.plot.width = 3, repr.plot.height = 4, repr.plot.res=300)
gg_sensitivity_filt <- ggplot(data=df_all_summary_filt_merged, aes(x=origin, y=transfer_type_prop)) +
    geom_col(fill='white', color='grey20') + facet_wrap(~cutoff, ncol = 1) +
    scale_y_continuous(labels=scales::percent) +
    labs(x='', y="Unstable OTU (%)") +
    theme_classic() +
    guides(x = guide_axis(angle = 45))

options(repr.plot.width=7.2, repr.plot.height=6, repr.plot.res=300)
layout <-'
ABC
DDD
'
p <- gg_sensitivity_filt + (fig7b & labs(x='', y="")) + fig7c + (fig7d & labs(x='', y="Unstable OTUs (%)")) + plot_layout(design = layout, heights = c(1.2, 1)) + plot_annotation(tag_levels = 'A')

figdir <- here::here('fig.outdir')
dir.create(figdir)

figfile <- here::here(figdir, 'fig.s10.nh_stability.pdf')
ggsave(figfile, p, width = 7.2, height = 6, dpi = 300, device = 'pdf')

#figfile <- here::here(figdir, 'fig.s10.nh_stability.png')
#ggsave(figfile, p, width = 7.2, height = 6, dpi = 300, device = 'png')

