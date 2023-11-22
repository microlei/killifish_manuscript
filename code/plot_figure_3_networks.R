# code for plotting the network

library(SpiecEasi)
library(igraph)
library(tidygraph)
library(ggraph)
library(graphlayouts)
library(RColorBrewer)

load(here("data/processed/se_results.RData"))

ps <- readRDS(here("data/processed/ps.rds"))
ps_NBH_wild <- ps %>% subset_samples(fishType=="New Bedford Harbor wild") %>% prune_taxa(taxa_sums(.)>50,.)
ps_SC_wild <- ps %>% subset_samples(fishType=="Scorton Creek wild") %>% prune_taxa(taxa_sums(.)>50,.)
ps_NBH_F2 <- ps %>% subset_samples(fishType=="New Bedford Harbor F2") %>% prune_zeros()
ps_SC_F2 <- ps %>% subset_samples(fishType== "Scorton Creek F2") %>% prune_zeros()

taxonomy <- tax_table(ps) %>% as.data.frame() %>% rownames_to_column(var="ASV")

se_SC_net <- adj2igraph(se_SC_refit, vertex.attr = list(name=taxa_names(ps_SC_wild)))
se_NBH_net <- adj2igraph(se_NBH_refit %>% symBeta(), vertex.attr = list(name=taxa_names(ps_NBH_wild)))
se_SC_F2_net <- adj2igraph(se_SC_F2_refit, vertex.attr = list(name=taxa_names(ps_SC_F2)))
se_NBH_F2_net <- adj2igraph(se_NBH_F2_refit, vertex.attr = list(name=taxa_names(ps_NBH_F2)))

# I used the top_ind from plot_networks_nbh_subsampled.R because it is more complete
#
# top_ind <- tribble(
#   ~"Topological Indices", ~"NBH wild", ~"SC wild", ~"NBH F2", ~"SC F2",
#   "No. ASVs (nodes)", length(V(se_NBH_net)), length(V(se_SC_net)), length(V(se_NBH_F2_net)), length(V(se_SC_F2_net)),
#   "No. interactions (edges)", length(E(se_NBH_net)), length(E(se_SC_net)), length(E(se_NBH_F2_net)), length(E(se_SC_F2_net)),
#   "Mean degree", round(mean(degree(se_NBH_net)),0), round(mean(degree(se_SC_net)),0), round(mean(degree(se_NBH_F2_net)),0), round(mean(degree(se_SC_F2_net)),0),
#   "SD degree", round(sd(degree(se_NBH_net)),0), round(sd(degree(se_SC_net)),0), round(sd(degree(se_NBH_F2_net)),0), round(sd(degree(se_SC_F2_net)),0),
#   "Transitivity/Cluster coefficient", transitivity(se_NBH_net), transitivity(se_SC_net), transitivity(se_NBH_F2_net), transitivity(se_SC_F2_net),
#   "Avg. path length", mean_distance(se_NBH_net), mean_distance(se_SC_net), mean_distance(se_NBH_F2_net), mean_distance(se_SC_F2_net),
#   "Density", edge_density(se_NBH_net), edge_density(se_SC_net), edge_density(se_NBH_F2_net), edge_density(se_SC_F2_net),
#   "Modularity (fast-greedy)", modularity(cluster_fast_greedy(se_NBH_net)), modularity(cluster_fast_greedy(se_SC_net)), modularity(cluster_fast_greedy(se_NBH_F2_net)), modularity(cluster_fast_greedy(se_SC_F2_net)),
#   "Degree centralization", centr_degree(se_NBH_net)$centralization, centr_degree(se_SC_net)$centralization, centr_degree(se_NBH_F2_net)$centralization, centr_degree(se_SC_F2_net)$centralization,
#   "Number of components, excluding singletons", sum(components(se_NBH_net)$csize>1), sum(components(se_SC_net)$csize>1), sum(components(se_NBH_F2_net)$csize>1), sum(components(se_SC_F2_net)$csize>1))
# 
# write_csv(top_ind, file=here("figures/top_ind.csv"))

se_SC_tidynet <- as_tbl_graph(se_SC_net) %>% activate(nodes) %>% mutate(trans_local = local_transitivity(), degree_local = local_ave_degree()) %>% 
  left_join(taxonomy, by=c("name"="ASV"))
se_NBH_tidynet <- as_tbl_graph(se_NBH_net) %>% activate(nodes) %>% mutate(trans_local = local_transitivity(), degree_local = local_ave_degree()) %>% 
  left_join(taxonomy, by=c("name"="ASV"))
se_SC_F2_tidynet <- as_tbl_graph(se_SC_F2_net) %>% activate(nodes) %>% mutate(trans_local = local_transitivity(), degree_local = local_ave_degree()) %>% 
  left_join(taxonomy, by=c("name"="ASV"))
se_NBH_F2_tidynet <- as_tbl_graph(se_NBH_F2_net) %>% activate(nodes) %>% mutate(trans_local = local_transitivity(), degree_local = local_ave_degree()) %>% 
  left_join(taxonomy, by=c("name"="ASV"))

graph_color_phyla <- tibble(breaks=unique(as_tibble(se_NBH_tidynet)$Phylum), values=c("#d25032","#576bd8","#5bc24f","#9459d2","#babb3a","#d265cf","#58952c","#c94097","#6cb66c","#e03d6b","#58c99e","#8d53a5","#958e2c","#a495e1","#dc913b","#536db1","#bdaf6c","#52a3d8","#a84b3d","#42c0c7","#b24669","#3d7f46","#dd87bb","#61712c","#95527d","#348a6e","#e4877b","#956a32"))

# Make layout using igraph's "nicely" algo
sc_layout <- create_layout(se_SC_tidynet, layout="stress")
sc_wild_graph <- ggraph(sc_layout) + geom_edge_link(edge_colour="grey66") + geom_node_point(aes(color=Phylum), size=.3) + scale_color_manual(breaks=graph_color_phyla$breaks, values=graph_color_phyla$values) + labs(title = "Scorton Creek wild") + theme(legend.position = "none")

# remove an outlying node to improve plotting
nbh_layout <- create_layout(se_NBH_tidynet %>% filter(name!="ASV702"), layout="stress")
nbh_wild_graph <- ggraph(nbh_layout) + geom_edge_link(edge_color="grey66", edge_alpha=0.15) + geom_node_point(aes(color=Phylum), size=.3) + scale_color_manual(breaks=graph_color_phyla$breaks, values=graph_color_phyla$values) + labs(title="New Bedford Harbor wild") + theme(legend.position = "none")

sc_F2_layout <- create_layout(se_SC_F2_tidynet, layout="stress")
sc_F2_graph <- ggraph(sc_F2_layout) + geom_edge_link(edge_colour="grey66") + geom_node_point(aes(color=Phylum), size=.3) + scale_color_manual(breaks=graph_color_phyla$breaks, values=graph_color_phyla$values) +labs(title="Scorton Creek F2") + theme(legend.position = "none")

nbh_F2_layout <- create_layout(se_NBH_F2_tidynet, layout="stress")
nbh_F2_graph <- ggraph(nbh_F2_layout) + geom_edge_link(edge_color="grey66", edge_alpha=0.15) + geom_node_point(aes(color=Phylum), size=.3) + scale_color_manual(breaks=graph_color_phyla$breaks, values=graph_color_phyla$values) + labs(title="New Bedford Harbor F2") + theme(legend.position = "none")

all_net <- plot_grid(sc_wild_graph, nbh_wild_graph, sc_F2_graph, nbh_F2_graph, labels = c("A", "B", "C", "D")) # this step takes a WHILE. be patient. Do NOT try to render it in RStudio!

ggsave(all_net,filename=here("figures/plot_networks.png"), width=5, height=5, units="in")


# plot degree distributions (no separate labels for each plot)

deg_dist <- tibble(fishType="Scorton Creek wild", deg = degree.distribution(se_SC_net)) %>% 
  bind_rows(tibble(fishType="Scorton Creek F2", deg = degree.distribution(se_SC_F2_net))) %>% 
  bind_rows(tibble(fishType="New Bedford Harbor wild", deg = degree.distribution(se_NBH_net))) %>% 
  bind_rows(tibble(fishType="New Bedford Harbor F2", deg = degree.distribution(se_NBH_F2_net))) %>% 
  group_by(fishType) %>% mutate(x = 0:(n()-1))

deg_dist$fishType_f <- factor(deg_dist$fishType, levels = c("Scorton Creek wild", "New Bedford Harbor wild", "Scorton Creek F2", "New Bedford Harbor F2"))

label_f <- labeller(fishType_f=labels.short)

p_deg_dist <- ggplot(deg_dist, aes(x=x, y=deg, color=fishType)) + geom_point(shape=1)+geom_line()+ facet_wrap(~fishType_f, scales="free", labeller = label_f) + scale_color_manual(values=cols.fishType, labels=labels.fishType) +
  labs(x="Degree", y="Frequency") + theme_light()+ theme(legend.position = "none", strip.text = element_text(size=10, color="black"), strip.background = element_rect(fill="white"))
ggsave(p_deg_dist, filename=here("figures/plot_degree_dist.png"), height=3, width =4, units = "in")

# plot degree distributions (separate labels for each plot)

p_sc_w <- tibble(fishType="Scorton Creek wild", deg = degree.distribution(se_SC_net)) %>% mutate(x = 0:(n()-1)) %>% 
  ggplot(aes(x=x, y=deg, color=fishType)) + geom_point(shape=1) + geom_line() + scale_color_manual(values = cols.fishType) + 
  labs(x="Degree", y="Frequency", title = "Scorton Creek wild") + theme_light() + theme(legend.position = "none", plot.title = element_text(hjust=0.5, face = "bold", size = 10))
p_nbh_w <- tibble(fishType="New Bedford Harbor wild", deg = degree.distribution(se_NBH_net)) %>% mutate(x = 0:(n()-1)) %>% 
  ggplot(aes(x=x, y=deg, color=fishType)) + geom_point(shape=1) + geom_line() + scale_color_manual(values = cols.fishType) + 
  labs(x="Degree", y="Frequency", title = "New Bedford Harbor wild") + theme_light() + theme(legend.position = "none", plot.title = element_text(hjust=0.5, face = "bold", size = 10))
p_sc_f2 <- tibble(fishType="Scorton Creek F2", deg = degree.distribution(se_SC_F2_net)) %>% mutate(x = 0:(n()-1)) %>% 
  ggplot(aes(x=x, y=deg, color=fishType)) + geom_point(shape=1) + geom_line() + scale_color_manual(values = cols.fishType) + 
  labs(x="Degree", y="Frequency", title = "Scorton Creek F2") + theme_light() + theme(legend.position = "none", plot.title = element_text(hjust=0.5, face = "bold", size = 10))
p_nbh_f2 <- tibble(fishType="New Bedford Harbor F2", deg = degree.distribution(se_NBH_F2_net)) %>% mutate(x = 0:(n()-1)) %>% 
  ggplot(aes(x=x, y=deg, color=fishType)) + geom_point(shape=1) + geom_line() + scale_color_manual(values = cols.fishType) + 
  labs(x="Degree", y="Frequency", title = "New Bedford Harbor F2") + theme_light() + theme(legend.position = "none", plot.title = element_text(hjust=0.5, face = "bold", size = 10))

p_deg_dist <- plot_grid(p_sc_w, p_nbh_w, p_sc_f2, p_nbh_f2, ncol=2, labels= c("E", "F", "G", "H"))

# Combine plots
network_plot <- plot_grid(all_net, p_deg_dist, ncol = 1)
ggsave(network_plot, filename=here("figures/plot_networks_combined.png"), width=5, height=8, units = "in")

