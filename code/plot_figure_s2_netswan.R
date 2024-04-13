library(tidyverse)
library(NetSwan)
library(SpiecEasi)
library(igraph)
library(here)

source(here("code/helpful_functions.R"))
load(here("data/processed/se_results.RData"))

ps <- readRDS(here("data/processed/ps.rds"))

se_SC_net <- adj2igraph(se_SC_refit)
se_NBH_net <- adj2igraph(se_NBH_refit %>% symBeta())

# Swan combinatory took a long time so these results were run once and then saved

# swan_NBH_wild <- swan_combinatory(se_NBH_net, 10)
# swan_SC_wild <- swan_combinatory(se_SC_net, 10)
# save(swan_NBH_wild, swan_SC_wild, file=here("data/processed/netswan_wild.RData"))

load(here("data/processed/netswan_wild.RData"))
colnames(swan_NBH_wild) <- colnames(swan_SC_wild) <- c("fraction", "betweenness", "degree", "cascade", "random")
swan_wild <- pivot_longer(swan_NBH_wild %>% as_tibble(), cols = betweenness:random, names_to = "attack") %>% add_column(fishType=rep("New Bedford Harbor wild", n=nrow(.))) %>% bind_rows(pivot_longer(swan_SC_wild %>% as_tibble(), cols=betweenness:random, names_to = "attack") %>% add_column(fishType=rep("Scorton Creek wild", n=nrow(.))))
# ggplot(swan_wild, aes(x=fraction, y=value, color=site)) + geom_point() + facet_wrap(~attack)

p<-ggplot(swan_wild %>% filter(attack=="degree"), aes(x=fraction, y=value, color=fishType)) + geom_point(size=4) + scale_color_manual(values = c(cols.fishType[1], cols.fishType[3])) + labs(x="Fraction of nodes removed", y=expression("Loss of connectivity ->")) + theme_light()
p <- p + theme(axis.title = element_text(size=15), legend.title = element_blank(), legend.position = c(.7,.2), legend.box.background = element_rect(linetype = 3))
ggsave(p, filename = here("figures/figure_s2.pdf"), width=4, height=4, units="in")

sc_blocks <- cohesive_blocks(se_SC_net)
nb_blocks <- cohesive_blocks(se_NBH_net)
