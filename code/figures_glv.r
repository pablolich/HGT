setwd("~/Desktop/phd/HGT/code")
require(ggplot2)
require(tidyverse)
library(ggpubr)
library(gridExtra)
library(tikzDevice)

data = read.csv('../data/abundances.csv')
data_plot = data %>% 
  pivot_longer(cols = c("ab", "ab_hgt")) %>% 
  group_by(X, name) %>% 
  mutate(mean = mean(value), sd = sd(value)) %>% 
  filter(sd > 0)

assembly = 
  ggplot(data = data_plot,
       aes(x = X, y = mean))+
  geom_line(aes(color = name), size = 0.8)+
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd,
                    group = name,
                    fill = name), 
              alpha = 0.3)+
  theme(aspect.ratio = 1,
        legend.position = c(0.3, 0.8),
        legend.spacing.y = unit(0,"cm"),
        legend.title = element_text(size = 21, hjust = 0.5, margin = margin(r = 7)),
        legend.text = element_text(size = 21),
        legend.key.width = unit(10, 'mm'),
        legend.key = element_rect(fill = NA, colour = NA),
        legend.background = element_rect(fill = NA),
        legend.margin = margin(b =  0),
        panel.background = element_blank(), 
        panel.border = element_rect(color = 'black', fill=NA),
        axis.title = element_text(size = 21), 
        axis.text = element_text(size = 15))+
  scale_colour_manual(values = c('#E1341E','#1ECBE1'),
                      labels = c('Mutations', 'Mutations + HGT'))+
  scale_fill_manual(values = c('#E1341E','#1ECBE1'),
                    labels = c('Mutations', 'Mutations + HGT'))+
  labs(x = 'Number of mutations', 
       y = 'Community richness',
       color = "",
       fill = "")

data_eigen = read_csv('../data/eigen.csv', col_names = F) 
data_eig = data_eigen %>% filter(X1 > -10)

eigen = ggplot(data = data_eig, aes(x = X1)) + 
  geom_histogram(aes(y=..density..), bins = 30,
                 color = 'black',
                 fill = '#E1341E')+
  theme(aspect.ratio = 1/2,
        legend.position = c(0.2, 0.8),
        legend.spacing.y = unit(0,"cm"),
        legend.title = element_text(size = 21, hjust = 0.5, margin = margin(r = 7)),
        legend.text = element_text(size = 21),
        legend.key.width = unit(10, 'mm'),
        legend.key = element_rect(fill = NA, colour = NA),
        legend.background = element_rect(fill = NA),
        legend.margin = margin(b =  0),
        panel.background = element_blank(), 
        panel.border = element_rect(color = 'black', fill=NA),
        axis.title = element_text(size = 21), 
        axis.text = element_text(size = 15))+
  labs(x = '', 
       y = '',
       colour = "",
       fill = "")



data_eigen_hgt = read_csv('../data/eigen_hgt.csv', col_names = F) %>% 
  filter(X1 > -10)
eigen_hgt = ggplot(data = data_eigen_hgt, aes(x = X1)) + 
  geom_histogram(aes(y = ..density..), bins = 30, 
                 fill = '#1ECBE1', color = 'black')+
  theme(aspect.ratio = 1/2,
        legend.position = c(0.2, 0.8),
        legend.spacing.y = unit(0,"cm"),
        legend.title = element_text(size = 21, hjust = 0.5, margin = margin(r = 7)),
        legend.text = element_text(size = 21),
        legend.key.width = unit(10, 'mm'),
        legend.key = element_rect(fill = NA, colour = NA),
        legend.background = element_rect(fill = NA),
        legend.margin = margin(b =  0),
        panel.background = element_blank(), 
        panel.border = element_rect(color = 'black', fill=NA),
        axis.title = element_text(size = 21), 
        axis.text = element_text(size = 15),
        axis.title.y = element_text(size = 21, hjust = 1.7, vjust = 2.7))+
  labs(x = 'Eigenvalues', 
       y = 'Density',
       colour = "",
       fill = "")


#Plot alltogether
pdf('../sandbox/figure_glv.pdf', width = 12.7, height = 5.71)
lay <- rbind(c(1, 1, 2, 2),
             c(1, 1, 3, 3))
grid.arrange(assembly, eigen, eigen_hgt, 
             layout_matrix = lay)
dev.off()
