# plotting the degree network and degree distribution of the subsampled NBH wild 
library(SpiecEasi)
library(tidygraph)
library(igraph)
library(ggraph)
library(here)
source(here("code/helpful_functions.R"))

load(here("data/processed/se_tests_refits.RData"))
ps <- readRDS(here("data/processed/ps.rds"))
ps_sc_pruned <- ps %>% subset_samples(fishType=="Scorton Creek wild") %>% prune_taxa(taxa_sums(.)>50,.)
ps_nbh_top <- ps  %>% subset_samples(fishType=="New Bedford Harbor wild")  %>% prune_taxa(names(sort(taxa_sums(.), TRUE))[1:ntaxa(ps_sc_pruned)],.)

se_NBH_top_net <- adj2igraph(se_test_nbh_top_gl %>% symBeta(), vertex.attr = list(name=taxa_names(ps_nbh_top)))
V(se_NBH_top_net)$degree <- degree(se_NBH_top_net)

taxonomy <- tax_table(ps) %>% as.data.frame() %>% rownames_to_column(var="ASV")

load(here("data/processed/se_results.RData"))

ps_NBH_wild <- ps %>% subset_samples(fishType=="New Bedford Harbor wild") %>% prune_taxa(taxa_sums(.)>50,.)
ps_SC_wild <- ps %>% subset_samples(fishType=="Scorton Creek wild") %>% prune_taxa(taxa_sums(.)>50,.)
ps_NBH_F2 <- ps %>% subset_samples(fishType=="New Bedford Harbor F2") %>% prune_zeros()
ps_SC_F2 <- ps %>% subset_samples(fishType== "Scorton Creek F2") %>% prune_zeros()

se_SC_net <- adj2igraph(se_SC_refit, vertex.attr = list(name=taxa_names(ps_SC_wild)))
se_NBH_net <- adj2igraph(se_NBH_refit %>% symBeta(), vertex.attr = list(name=taxa_names(ps_NBH_wild)))
se_SC_F2_net <- adj2igraph(se_SC_F2_refit, vertex.attr = list(name=taxa_names(ps_SC_F2)))
se_NBH_F2_net <- adj2igraph(se_NBH_F2_refit, vertex.attr = list(name=taxa_names(ps_NBH_F2)))

# this is just to make the phylum coloring consistent with the main graph
se_NBH_tidynet <- as_tbl_graph(se_NBH_net) %>% activate(nodes) %>% mutate(trans_local = local_transitivity(), degree_local = local_ave_degree()) %>% 
  left_join(taxonomy, by=c("name"="ASV"))
graph_color_phyla <- tibble(breaks=unique(as_tibble(se_NBH_tidynet)$Phylum), values=c("#d25032","#576bd8","#5bc24f","#9459d2","#babb3a","#d265cf","#58952c","#c94097","#6cb66c","#e03d6b","#58c99e","#8d53a5","#958e2c","#a495e1","#dc913b","#536db1","#bdaf6c","#52a3d8","#a84b3d","#42c0c7","#b24669","#3d7f46","#dd87bb","#61712c","#95527d","#348a6e","#e4877b","#956a32"))

nbh_top_ind <-tribble(
  ~"Topological Indices", ~"NBH wild", ~"SC wild", ~"NBH F2", ~"SC F2", ~"NBH wild subsampled",
  "No. ASVs (nodes)", length(V(se_NBH_net)), length(V(se_SC_net)), length(V(se_NBH_F2_net)), length(V(se_SC_F2_net)), length(V(se_NBH_top_net)), 
  "No. interactions (edges)", length(E(se_NBH_net)), length(E(se_SC_net)), length(E(se_NBH_F2_net)), length(E(se_SC_F2_net)), length(E(se_NBH_top_net)),
  "Mean degree", round(mean(degree(se_NBH_net)),0), round(mean(degree(se_SC_net)),0), round(mean(degree(se_NBH_F2_net)),0), round(mean(degree(se_SC_F2_net)),0), round(mean(degree(se_NBH_top_net)),0),
  "SD degree", round(sd(degree(se_NBH_net)),0), round(sd(degree(se_SC_net)),0), round(sd(degree(se_NBH_F2_net)),0), round(sd(degree(se_SC_F2_net)),0), round(sd(degree(se_SC_F2_net)),0),
  "Transitivity/Cluster coefficient", transitivity(se_NBH_net), transitivity(se_SC_net), transitivity(se_NBH_F2_net), transitivity(se_SC_F2_net), transitivity(se_NBH_top_net),
  "Avg. path length", mean_distance(se_NBH_net), mean_distance(se_SC_net), mean_distance(se_NBH_F2_net), mean_distance(se_SC_F2_net), mean_distance(se_NBH_top_net),
  "Density", edge_density(se_NBH_net), edge_density(se_SC_net), edge_density(se_NBH_F2_net), edge_density(se_SC_F2_net), edge_density(se_NBH_top_net),
  "Modularity (fast-greedy)", modularity(cluster_fast_greedy(se_NBH_net)), modularity(cluster_fast_greedy(se_SC_net)), modularity(cluster_fast_greedy(se_NBH_F2_net)), modularity(cluster_fast_greedy(se_SC_F2_net)), modularity(cluster_fast_greedy(se_NBH_top_net)),
  "Degree centralization", centr_degree(se_NBH_net)$centralization, centr_degree(se_SC_net)$centralization, centr_degree(se_NBH_F2_net)$centralization, centr_degree(se_SC_F2_net)$centralization, centr_degree(se_NBH_top_net)$centralization,
  "Number of components, excluding singletons", sum(components(se_NBH_net)$csize>1), sum(components(se_SC_net)$csize>1), sum(components(se_NBH_F2_net)$csize>1), sum(components(se_SC_F2_net)$csize>1), sum(components(se_NBH_top_net)$csize>1))

write_csv(nbh_top_ind, file=here("figures/top_ind.csv"))

se_NBH_top_tidynet <- as_tbl_graph(se_NBH_top_net) %>% activate(nodes) %>% mutate(trans_local = local_transitivity(), degree_local = local_ave_degree()) %>% 
  left_join(taxonomy, by=c("name"="ASV"))

nbh_top_layout <- create_layout(se_NBH_top_tidynet, layout="stress")
nbh_wild_top_graph <- ggraph(nbh_top_layout) + geom_edge_link(edge_color="grey66", edge_alpha=0.15) + geom_node_point(aes(color=Phylum), size=.3) + scale_color_manual(breaks=graph_color_phyla$breaks, values=graph_color_phyla$values) + labs(title="New Bedford Harbor wild subsampled") + theme(legend.position = "none")

nbh_wild_top_graph
ggsave(nbh_wild_top_graph, filename=here("figures/plot_network_subsampled.png"))

nbh_deg_dist <- tibble(fishType = "New Bedford Harbor wild subsample", deg = degree.distribution(se_NBH_top_net))%>% mutate(x = 0:(n()-1))

p_nbh_deg_dist <- ggplot(nbh_deg_dist, aes(x=x, y=deg, color=fishType)) + geom_point(shape=1)+geom_line() + labs(x="Degree", y="Frequency") + theme_light() + theme(legend.position="none")

plot_grid(nbh_wild_top_graph, p_nbh_deg_dist, labels = c("A", "B"))
ggsave(filename=here("figures/plot_network_subsampled.png"), height=3, width=6, units="in")
