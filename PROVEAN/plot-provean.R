#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(ggthemes)

dat <- fread('provean-out.txt', skip='# VARIATION')
setnames(dat, c('var','score'))

dat[, var := gsub('del','', var)]
dat[, var := gsub('[A-Z]','', var)]

dat[, c('del_start','del_stop') := tstrsplit(var, split='_')]
dat[is.na(del_stop), del_stop := del_start]
dat[, del_start := as.numeric(del_start)]
dat[, del_stop := as.numeric(del_stop)]



g.dels <- ggplot(dat[del_start>1], aes(x = del_start, y = del_stop, fill = score)) +
  geom_tile() +  
  scale_fill_gradientn(colors = c("#a020f0", "#a020f0", "#ffff00")) +
  scale_x_continuous(breaks = seq(0,90,10)) +
  scale_y_continuous(breaks =seq(0,90,10)) +
  theme_few() +
  labs(title='PROVEAN score for deletion mutants')

ggsave(g.dels, file='provean-deletions.png', width=20, height=20, units='cm')
ggsave(g.dels, file='provean-deletions.pdf', width=20, height=20, units='cm')