#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggthemes)


# Missense variants
missense <- fread('esm-missense.csv')

missense[, pos := gsub('[A-z]', '', mut_name)]
missense[, pos := as.numeric(pos)]

veps <- missense[, .('mean_esm'=mean(esm_score)), by=pos]

g.missense <- ggplot(veps[pos>1], aes(x=pos,y=1, fill=mean_esm)) + geom_tile() +
    scale_x_continuous(breaks=seq(0,90,10)) +
    labs(title='mean esm score across all missense variants') +
    ylab('') +
    theme_few() +
    scale_fill_gradientn(colors = c("#a020f0", "#a020f0", "#ffff00"))

ggsave(g.missense, file='esm-missense.png', width=20, height=5, units='cm')
ggsave(g.missense, file='esm-missense.pdf', width=20, height=5, units='cm')


# Deletions
dat <- fread('esm-deletions.csv')
setnames(dat, 'start_pos','start')

dat[, mut_length := nchar(mut_seq)]
dat[, del_length := 87-mut_length]
dat[, end_value := start + del_length]


g.dels <- ggplot(dat[start>1], aes(x = start, y = end_value, fill = esm_score)) +
  geom_tile() +  
  scale_fill_gradientn(colors = c("#a020f0", "#a020f0", "#ffff00")) +
  scale_x_continuous(breaks = seq(0,90,10)) +
  scale_y_continuous(breaks =seq(0,90,10)) +
  theme_few() +
  labs(title='esm score for deletion mutants')

ggsave(g.dels, file='esm-deletions.png', width=20, height=20, units='cm')
ggsave(g.dels, file='esm-deletions.pdf', width=20, height=20, units='cm')