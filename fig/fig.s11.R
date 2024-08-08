source(here::here('setup.R'))


### unstable neighborhood examples
require(gggenes)

tab <- 'som-data/fig-data/stability/example/IS.unstable.OTU998629.dram.10000_10000.4plot.add_allinfo.ge1k.curate.tsv'

df <- read.table(tab, sep='\t', header=T) %>%
        #dplyr::filter(contig_length>=10000) %>%
	dplyr::filter(!is.na(DepthAvg)) %>%
        dplyr::mutate(DepthLumping =  ifelse(DepthAvg >= 0 & DepthAvg < 10, "0-9", 
                               ifelse(DepthAvg >= 10 & DepthAvg < 20, "10-19", 
                               ifelse(DepthAvg >= 20 & DepthAvg < 30, "20-29",
                               ifelse(DepthAvg >= 30 & DepthAvg < 40, "30-39",
                               ifelse(DepthAvg >= 40 & DepthAvg < 50, "40-49",
                               ifelse(DepthAvg >= 50 & DepthAvg < 60, "50-59",
                               ifelse(DepthAvg >= 60 & DepthAvg < 70, "60-69",
                               ifelse(DepthAvg >= 70 & DepthAvg < 80, "70-79",
                               "80+"))))))))) %>%
        #dplyr::mutate(genome=stringr::str_c(Habitat, Year, DepthLumping, contig_length, phylum, sep='|'))
        dplyr::mutate(genome=stringr::str_c(Habitat, Year, contig_length, sep='|'), phylum=stringr::str_remove(phylum, 'p__')) %>%
        dplyr::mutate(genome=stringr::str_c(genome, phylum, sep='\n'))


dummies <- make_alignment_dummies(
  df,
  aes(xmin = start, xmax = end, y = genome, id = gene),
  #on = "integrase",
    on = 'transposase'
)

#colors <- c(RColorBrewer::brewer.pal(8, 'Dark2'), RColorBrewer::brewer.pal(10, 'Paired'))
old_colors <- c(RColorBrewer::brewer.pal(8, 'Dark2')[1:7], RColorBrewer::brewer.pal(10, 'Paired'))
gene_levels <- df %>% pull(gene) %>% factor() %>% levels()
select_gene_vec <- c('int', 'intI', 'relaxase', 'transposase', 'transposase2',  'unknown')
other_gene_levels <- gene_levels[!gene_levels %in% select_gene_vec]
other_colors <- old_colors[1:length(other_gene_levels)]
new_gene_levels <- c('transposase', 'transposase2', 'unknown') %>% c(other_gene_levels)
colors <- c('grey10', 'grey70', 'grey100') %>% c(other_colors)
names(colors) <- new_gene_levels
df <- df %>% mutate(gene = factor(gene, levels = new_gene_levels))

options(repr.plot.width=7, repr.plot.height=2, repr.plot.res=300)
gg <- (ggplot(data=df, aes(xmin=start, xmax=end, y=genome, fill=gene, forward=strand))
       + geom_gene_arrow()
       + geom_blank(data=dummies, aes(forward=1))
       + facet_wrap(~ genome, scales='free', ncol=1)
       + scale_fill_manual(name='', values=colors)
       + labs(y="")
       + theme_genes()
       #+ theme(legend.position = 'top', legend.key.size = unit(7, 'pt'))
       + theme(legend.key.size = unit(10, 'pt'))
       )

fig.s8.a <- gg


gg <- (ggplot(data=df, aes(xmin=start, xmax=end, y=genome, fill=gene, forward=strand))
       + geom_gene_arrow()
       #+ facet_wrap(~ contig, ncol=1)
       + scale_fill_manual(values=colors)
       + labs(y="")
       + theme_genes())

tab <- 'som-data/fig-data/stability/example/Phage.unstable.OTU427098.dram.10000_10000.4plot.add_allinfo.ge1k.curate.tsv'

df <- read.table(tab, sep='\t', header=T) %>%
        #dplyr::filter(contig_length>=10000) %>%
	dplyr::filter(!is.na(DepthAvg)) %>%
        dplyr::mutate(DepthLumping =  ifelse(DepthAvg >= 0 & DepthAvg < 10, "0-9", 
                               ifelse(DepthAvg >= 10 & DepthAvg < 20, "10-19", 
                               ifelse(DepthAvg >= 20 & DepthAvg < 30, "20-29",
                               ifelse(DepthAvg >= 30 & DepthAvg < 40, "30-39",
                               ifelse(DepthAvg >= 40 & DepthAvg < 50, "40-49",
                               ifelse(DepthAvg >= 50 & DepthAvg < 60, "50-59",
                               ifelse(DepthAvg >= 60 & DepthAvg < 70, "60-69",
                               ifelse(DepthAvg >= 70 & DepthAvg < 80, "70-79",
                               "80+"))))))))) %>%
        #dplyr::mutate(genome=stringr::str_c(Habitat, Year, DepthLumping, contig_length, phylum, sep='|')) %>%
        dplyr::mutate(genome=stringr::str_c(Habitat, Year, contig_length, sep='|'), phylum=stringr::str_remove(phylum, 'p__')) %>%
        dplyr::mutate(genome=stringr::str_c(genome, phylum, sep='\n'))


dummies <- make_alignment_dummies(
  df,
  aes(xmin = start, xmax = end, y = genome, id = gene),
  #on = "integrase",
    on = 'int'
)

#colors <- c(RColorBrewer::brewer.pal(8, 'Dark2'), RColorBrewer::brewer.pal(10, 'Paired'))
old_colors <- c(RColorBrewer::brewer.pal(8, 'Dark2')[1:7], RColorBrewer::brewer.pal(10, 'Paired'))
gene_levels <- df %>% pull(gene) %>% factor() %>% levels()
select_gene_vec <- c('int', 'intI', 'relaxase', 'transposase', 'transposase2',  'unknown')
other_gene_levels <- gene_levels[!gene_levels %in% select_gene_vec]
other_colors <- old_colors[1:length(other_gene_levels)]
new_gene_levels <- c('int', 'unknown') %>% c(other_gene_levels)
colors <- c('grey10', 'grey100') %>% c(other_colors)

names(colors) <- new_gene_levels
df <- df %>% mutate(gene = factor(gene, levels = new_gene_levels))

options(repr.plot.width=7, repr.plot.height=3, repr.plot.res=300)
gg <- (ggplot(data=df, aes(xmin=start, xmax=end, y=genome, fill=gene, forward=strand))
       + geom_gene_arrow()
       + geom_blank(data=dummies, aes(forward=1))
       + facet_wrap(~ genome, scales='free', ncol=1)
       + scale_fill_manual(name='', values=colors)
       + labs(y="")
       + theme_genes()
       + theme(legend.key.size = unit(10, 'pt'))
      ) 

fig.s8.b <- gg

gg <- (ggplot(data=df, aes(xmin=start, xmax=end, y=genome, fill=gene, forward=strand))
       + geom_gene_arrow()
       #+ facet_wrap(~ contig, ncol=1)
       + scale_fill_manual(values=colors)
       + labs(y="")
       + theme_genes())



tab <- 'som-data/fig-data/stability/example/CE.unstable.OTU112905.dram.10000_10000.4plot.add_allinfo.ge1k.curate.tsv'

df <- read.table(tab, sep='\t', header=T) %>%
        #dplyr::filter(contig_length>=10000) %>%
	dplyr::filter(!is.na(DepthAvg)) %>%
        dplyr::mutate(DepthLumping =  ifelse(DepthAvg >= 0 & DepthAvg < 10, "0-9", 
                               ifelse(DepthAvg >= 10 & DepthAvg < 20, "10-19", 
                               ifelse(DepthAvg >= 20 & DepthAvg < 30, "20-29",
                               ifelse(DepthAvg >= 30 & DepthAvg < 40, "30-39",
                               ifelse(DepthAvg >= 40 & DepthAvg < 50, "40-49",
                               ifelse(DepthAvg >= 50 & DepthAvg < 60, "50-59",
                               ifelse(DepthAvg >= 60 & DepthAvg < 70, "60-69",
                               ifelse(DepthAvg >= 70 & DepthAvg < 80, "70-79",
                               "80+"))))))))) %>%
        #dplyr::mutate(genome=stringr::str_c(Habitat, Year, DepthLumping, contig_length, phylum, sep='|')) %>%
        dplyr::mutate(genome=stringr::str_c(Habitat, Year, contig_length, sep='|'), phylum=stringr::str_remove(phylum, 'p__')) %>%
        dplyr::mutate(genome=stringr::str_c(genome, phylum, sep='\n'))


dummies <- make_alignment_dummies(
  df,
  aes(xmin = start, xmax = end, y = genome, id = gene),
  #on = "integrase",
    on = 'relaxase'
)

old_colors <- c(RColorBrewer::brewer.pal(8, 'Dark2')[1:7], RColorBrewer::brewer.pal(10, 'Paired'))
gene_levels <- df %>% pull(gene) %>% factor() %>% levels()
select_gene_vec <- c('int', 'intI', 'relaxase', 'transposase', 'transposase2',  'unknown')
other_gene_levels <- gene_levels[!gene_levels %in% select_gene_vec]
other_colors <- old_colors[1:length(other_gene_levels)]
new_gene_levels <- c('relaxase', 'transposase', 'unknown') %>% c(other_gene_levels)
colors <- c('grey10', 'grey70', 'grey100') %>% c(other_colors)
names(colors) <- new_gene_levels

df <- df %>% mutate(gene = factor(gene, levels = new_gene_levels))

options(repr.plot.width=7.2, repr.plot.height=2, repr.plot.res=300)
gg <- (ggplot(data=df, aes(xmin=start, xmax=end, y=genome, fill=gene, forward=strand))
       + geom_gene_arrow()
       + geom_blank(data=dummies, aes(forward=1))
       + facet_wrap(~ genome, scales='free', ncol=1)
       + scale_fill_manual(name='', values=colors)
       + labs(y="")
       + theme_genes()
       + theme(legend.key.size = unit(10, 'pt'))
      ) 

fig.s8.c <- gg

gg <- (ggplot(data=df, aes(xmin=start, xmax=end, y=genome, fill=gene, forward=strand))
       + geom_gene_arrow()
       #+ facet_wrap(~ contig, ncol=1)
       + scale_fill_manual(values=colors)
       + labs(y="")
       + theme_genes())


tab <- 'som-data/fig-data/stability/example/Integron.unstable.OTU799083.dram.10000_10000.4plot.add_allinfo.ge1k.curate.tsv'

df <- read.table(tab, sep='\t', header=T) %>%
        #dplyr::filter(contig_length>=10000) %>%
	dplyr::filter(!is.na(DepthAvg)) %>%
        dplyr::mutate(DepthLumping =  ifelse(DepthAvg >= 0 & DepthAvg < 10, "0-9", 
                               ifelse(DepthAvg >= 10 & DepthAvg < 20, "10-19", 
                               ifelse(DepthAvg >= 20 & DepthAvg < 30, "20-29",
                               ifelse(DepthAvg >= 30 & DepthAvg < 40, "30-39",
                               ifelse(DepthAvg >= 40 & DepthAvg < 50, "40-49",
                               ifelse(DepthAvg >= 50 & DepthAvg < 60, "50-59",
                               ifelse(DepthAvg >= 60 & DepthAvg < 70, "60-69",
                               ifelse(DepthAvg >= 70 & DepthAvg < 80, "70-79",
                               "80+"))))))))) %>%
        #dplyr::mutate(genome=stringr::str_c(Habitat, Year, DepthLumping, contig_length, phylum, sep='|')) %>%
        dplyr::mutate(genome=stringr::str_c(Habitat, Year, contig_length, sep='|'), phylum=stringr::str_remove(phylum, 'p__')) %>%
        dplyr::mutate(genome=stringr::str_c(genome, phylum, sep='\n'))


dummies <- make_alignment_dummies(
  df,
  aes(xmin = start, xmax = end, y = genome, id = gene),
  #on = "integrase",
    on = 'intI'
)

old_colors <- c(RColorBrewer::brewer.pal(8, 'Dark2')[1:7], RColorBrewer::brewer.pal(10, 'Paired'))
gene_levels <- df %>% pull(gene) %>% factor() %>% levels()
select_gene_vec <- c('int', 'intI', 'relaxase', 'transposase', 'transposase2',  'unknown')
other_gene_levels <- gene_levels[!gene_levels %in% select_gene_vec]
other_colors <- old_colors[1:length(other_gene_levels)]
new_gene_levels <- c('intI', 'transposase', 'unknown') %>% c(other_gene_levels)
colors <- c('grey10', 'grey70', 'grey100') %>% c(other_colors)
names(colors) <- new_gene_levels
df <- df %>% mutate(gene = factor(gene, levels = new_gene_levels))

options(repr.plot.width=7.2, repr.plot.height=2, repr.plot.res=300)
gg <- (ggplot(data=df, aes(xmin=start, xmax=end, y=genome, fill=gene, forward=strand))
       + geom_gene_arrow()
       + geom_blank(data=dummies, aes(forward=1))
       + facet_wrap(~ genome, scales='free', ncol=1)
       + scale_fill_manual(name='', values=colors)
       + labs(y="")
       + theme_genes()
       + theme(legend.key.size = unit(10, 'pt'))
      ) 

fig.s8.d <- gg


gg <- (ggplot(data=df, aes(xmin=start, xmax=end, y=genome, fill=gene, forward=strand))
       + geom_gene_arrow()
       #+ facet_wrap(~ contig, ncol=1)
       + scale_fill_manual(values=colors)
       + labs(y="")
       + theme_genes())


options(repr.plot.width=7.2, repr.plot.height=8, repr.plot.res=300)
p <- fig.s8.a + fig.s8.b + fig.s8.c + fig.s8.d + plot_layout(heights = c(2,3,1.5,1.5)) + plot_annotation(tag_levels = 'A')

figdir <- 'fig.outdir'
dir.create(figdir)

figfile <- here::here(figdir, 'fig.s11.nh_stability.example.pdf')
ggsave(figfile, p, width = 7.2, height = 8, dpi = 300, device = 'pdf')

#figfile <- here::here(figdir, 'fig.s11.nh_stability.example.png')
#ggsave(figfile, p, width = 7.2, height = 8, dpi = 300, device = 'png')


