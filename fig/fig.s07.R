source(here::here('setup.R'))


### coverM (minimap2-sr) with sheared contigs as coverage and with shared contig 100ani as separating "missed assembly" and "missed_binning"
### NOTE: adding 100ani or not has minimal impact on coverage results;

###
# fig.S5.ABD
###

gene_f <- 'som-data/fig-data/short_vs_long/gene_mapping.allinfo.final.add_longonly.tsv'

df_gene <- read_tsv(gene_f, col_types=cols()) %>%
    #dplyr::filter((contig_len>=0) & ((contig_len-end_bp)>=500) & (start_bp>=500)) %>%
    #dplyr::filter(!is.na(breadth_cov_contig_bam_by_sam)) %>%
    dplyr::mutate(long_only=if_else(long_only==TRUE, 'missed', 'covered')) %>%
    tidyr::separate(classification, into = c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'), sep=';') %>%
    dplyr::mutate(genome = stringr::str_c(genome, phylum, sep='\n'))

### separate those missed by assembly and binning

df_gene2 <- df_gene %>% 
    dplyr::mutate(long_only=if_else(long_only=='missed' & breadth_cov_contig_sheared_sr_100ani >=1, stringr::str_c('missed', 'binning', sep='\n'), 
                                    if_else(long_only=='missed', stringr::str_c('missed', 'assembly', sep='\n'), long_only)))


comp <- list(c('covered', 'missed\nassembly'), c('covered', 'missed\nbinning'), c('missed\nassembly', 'missed\nbinning'))
options(repr.plot.width=12, repr.plot.height=3, repr.plot.res=300)
gg1_1 <- (ggplot(data=df_gene2, aes(x=long_only, y=depth_cov_read)) + geom_boxplot(outlier.shape = 21, outlier.alpha = 0.5, outlier.size = 0.5) 
        + geom_hline(yintercept = 10, color='red', linetype='longdash')
        #+ stat_summary(fun = mean, geom='point', shape=4, color='red')
        + theme_classic()
          +  scale_y_log10()
          + annotation_logticks(sides = 'l')
        #  + guides(x=guide_axis(angle = 45))
        #+ ggpubr::stat_compare_means(label='p.signif', comparisons = comp, vjust=1.5)
        #+ ggpubr::stat_compare_means(label='p.signif', ref.group='covered')
        + labs(x='', y = 'Read depth\ncoverage')
        )

comp <- list(c('covered', 'missed\nassembly'))
gg1_2 <- (ggplot(data=df_gene2, aes(x=long_only, y=breadth_cov_read)) + geom_boxplot(outlier.shape = 21, outlier.alpha = 0.5, outlier.size = 0.5) 
        #+ geom_hline(yintercept = 1, color='red', linetype='dashed')
        #+ stat_summary(fun=mean, geom='point', shape=4, color='red')
        + theme_classic()
        #  + guides(x=guide_axis(angle = 45))
        + ggpubr::stat_compare_means(label='p.signif', comparisons = comp, vjust=0.5)
        #+ ggpubr::stat_compare_means(label='p.signif', ref.group='covered')
        + scale_y_continuous(labels = scales::percent_format(accuracy = 1))
        + labs(x='', y = 'Read breadth\ncoverage')
        )

comp <- list(c('covered', 'missed\nassembly'), c('covered', 'missed\nbinning'))
gg2_1 <- (ggplot(data=df_gene2, aes(x=long_only, y=depth_cov_contig_sheared_sr)) + geom_boxplot(outlier.shape = 21, outlier.alpha = 0.5, outlier.size = 0.5) 
        + theme_classic()
        + ggpubr::stat_compare_means(label='p.signif', comparisons = comp, vjust=0.5)
        #+ ggpubr::stat_compare_means(label='p.signif', ref.group='covered')
        + labs(x='', y = 'Contig depth\ncoverage')
        )


comp <- list(c('covered', 'missed\nassembly'))
gg2_2 <- (ggplot(data=df_gene2, aes(x=long_only, y=breadth_cov_contig_sheared_sr)) + geom_boxplot(outlier.shape = 21, outlier.alpha = 0.5, outlier.size = 0.5) 
        + theme_classic()
        + ggpubr::stat_compare_means(label='p.signif', comparisons = comp, vjust=0.5)
        #+ ggpubr::stat_compare_means(label='p.signif', ref.group='covered')
        + scale_y_continuous(labels = scales::percent_format(accuracy = 1))
        #  + guides(x=guide_axis(angle = 45))
        + labs(x='', y = 'Contig breadth\ncoverage')
        )

comp <- list(c('covered', 'missed\nassembly'))
gg3 <- (ggplot(data=df_gene2, aes(x=long_only, y=pi)) 
        + geom_boxplot(outlier.shape = 21, outlier.alpha = 0.5, outlier.size = 0.5) 
        #+ ggbeeswarm::geom_beeswarm(shape=21, size=1, alpha=0.6)
        #+ stat_summary(fun = mean, geom = 'point', shape=4, color='red')
        + theme_classic()
        #+ guides(x=guide_axis(angle = 45))
        + ggpubr::stat_compare_means(label='p.signif', comparisons = comp, vjust=0.5)
        #+ ggpubr::stat_compare_means(label='p.signif', ref.group='covered')
        + labs(x='', y='Nucleotide\ndiversity')
        )


fig.s3a <- gg1_1 | gg1_2
fig.s3a <- wrap_elements(full=fig.s3a)

fig.s3b <- gg2_1 | gg2_2
fig.s3b <- wrap_elements(full=fig.s3b)

fig.s3d <- gg3


###
# fig.S5.C
###

### genomic context variation

gene_f <- 'som-data/fig-data/short_vs_long/gene_mapping.allinfo.final.add_longonly.tsv'
clust_f <- 'som-data/fig-data/short_vs_long/stability/gene.prefilter.hmmtblout.tophit.cdhit1d0.s_eq_1.faa.clstr.tsv'
nh_f <- 'som-data/fig-data/short_vs_long/stability/all.gene_neighbor.blastn.pw_iden_within_cluster.sep_ends_sum.unstable_cnt.add_nh_per_otu.tsv'

df_gene <- read_tsv(gene_f, col_types=cols()) %>%
    dplyr::mutate(long_only=if_else(long_only==TRUE, 'missed', 'covered'))

df_clust <- read_tsv(clust_f, col_types=cols())
df_nh <-read_tsv(nh_f, col_types=cols())

df <- df_gene %>% 
    tidyr::separate(col = classification, 
                    into = c('domain', 'phylum', 'class', 'order', 'family', 'genus', 'species'),
                    sep = ';') %>%
    dplyr::mutate(genome=stringr::str_c(genome, phylum, sep='\n')) %>%
    dplyr::arrange(phylum) %>%
    dplyr::mutate(genome=factor(genome, levels = unique(genome))) %>%
    dplyr::mutate(long_only = if_else(breadth_cov_contig_sheared_sr>=1 & long_only=='missed', stringr::str_c('missed', 'binning', sep='\n'),
                                     if_else(long_only=='missed', stringr::str_c('missed', 'assembly', sep='\n'), long_only))) %>%
    dplyr::left_join(df_clust, by = c('recombinase'='mem')) %>%
    dplyr::inner_join(df_nh, by = c('OTU')) %>%
    dplyr::filter(total_cnt>=5)


comp <- list(c('covered', 'missed\nassembly', 'missed\nbinning'))

options(repr.plot.width=3, repr.plot.height=3, repr.plot.res=300)
gg <- (ggplot(data=df, aes(x=long_only, y=uniq_cnt)) 
       + geom_boxplot(outlier.shape = 21, outlier.alpha = 0.5, outlier.size = 0.5)
       + labs(x="", y="Genomic context\nnumber")
       + theme_classic()
       + ggpubr::stat_compare_means(label='p.signif', comparisons = comp, vjust=0.5) 
      )

fig.s3c <- gg

###
# fig.S5.E
###

### it's manually made

options(repr.plot.width=6, repr.plot.height=8.5, repr.plot.res=300)
layout <-'
AA
BB
CD
EE'

fig.s3e <- ggplot() + theme_void()

p <- fig.s3a + fig.s3b + fig.s3c + fig.s3d + fig.s3e + plot_layout(design = layout, heights = c(2,2,1.5,2)) + plot_annotation(tag_levels = 'A')

figdir <- here::here('fig.outdir')
figfile <- here::here(figdir, 'fig.s05.short_vs_hybrid.pdf')
ggsave(figfile, p, width = 6, height = 6.5, dpi = 300, device = 'pdf')

#figfile <- here::here(figdir, 'fig.s05.short_vs_hybrid.png')
#ggsave(figfile, p, width = 6, height = 6.5, dpi = 300, device = 'png' )
