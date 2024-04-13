library(here)
source(here("code/plot_diversity_metrics.R"))
source(here("code/plot_pca_plain.R"))

# the diversity plot
p_full <- plot_grid(p_alpha, p_beta, labels= c("A", "B"))/p_legend + plot_layout(heights=c(1,.1))

# the pca plot
pca_legend <- get_legend(pca.plot + theme(legend.position = "bottom") + guides(color = guide_legend(nrow=2)))
pca_full <- plot_grid(wild_subplot, F2_subplot, ncol=1, labels = c("C", "D")) %>% plot_grid(pca.plot + theme(legend.position = "none"), rel_widths = c(.5,1), labels=c("", "E"))/pca_legend + plot_layout(heights=c(1,.1))

# the combined plot
plot_grid(p_full, pca_full, ncol=1)

ggsave(filename=here("figures/figure_1.pdf"), device="pdf", width=8, height=11, units = "in")
