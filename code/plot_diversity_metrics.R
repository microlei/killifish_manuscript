# Code to plot the alpha and beta diversity

library(tidyverse)
library(phyloseq)
library(vegan)
library(cowplot)
library(patchwork)
library(broom)
library(here)

source(here("code/helpful_functions.R"))

ps <- readRDS(here("data/processed/ps.rds"))
# sample_data(ps)$site_short <- ifelse(sample_data(ps)$site=="Scorton Creek", "SC", "NBH")
# sample_data(ps)$site_short <- factor(sample_data(ps)$site_short, levels=c("NBH", "SC"))
# sample_data(ps)$wild_or_F2 <- factor(sample_data(ps)$wild_or_F2, levels=c("wild", "F2", "water"))
# 
p <- ps %>% subset_samples(sampleType=="gut") %>% plot_richness(x = "fishType", measures=c("Observed","Shannon","Simpson"))
w <- ps %>% subset_samples(sampleType=="water") %>% plot_richness(x = "fishType", measures=c("Observed","Shannon","Simpson"))

## make the table of alpha diversity values
bind_rows(p$data, w$data) %>% group_by(fishType, variable) %>% 
  summarise(Mean = str_c(format(mean(value), digits=3, scientific=FALSE), " (", n(), ")")) %>% 
  pivot_wider(id_cols=fishType, names_from=variable, values_from = Mean) %>% 
  rename("Sample type" = fishType) %>% 
  write_csv(here("figures/alpha_diversity.csv"))

# p_plot <- ggplot(p$data, aes(x=site_short, y=value, color=wild_or_F2)) + geom_boxplot() + scale_color_manual(name="Fish type", values=cols.wild_or_F2)+geom_point(data=w$data, aes(x=site_short, y=value), color="black") + 
#   facet_wrap(~variable, scales="free_y") + labs(title="Alpha diversity metrics", y="", x="Site") + theme(legend.position = "none")
# 
# lw <- ggplot(w$data, aes(x=site, y=value, color=site))+geom_point()+scale_color_manual(values=c("black", "black"), labels=c("Water"), breaks=c("New Bedford Harbor"),name="")+theme(legend.position = "top")
# legend_w <- get_legend(lw + theme(legend.position = "bottom", legend.margin = margin(), legend.box.margin = margin(), legend.spacing = unit(0,"cm")))
# 
# legend_p <- get_legend(p_plot + theme(legend.position = "bottom", legend.margin = margin(), legend.box.margin = margin(), legend.spacing = unit(0,"cm")))
# 
# p_plot/legend_p/legend_w + plot_layout(heights=c(1,.05,.05))
# 
# ggsave(filename=here("figures/plot_alpha_diversity.png"))

# plot for presentation
p$data <- p$data %>% mutate(variable = recode(variable, "Observed" = "Richness", "Shannon" = "Shannon Index", "Simpson" = "Simpson Index"))
p$data$wild_or_F2 <- factor(p$data$wild_or_F2, levels = c("wild", "F2"))
#p$data$fishType <- factor(p$data$fishType, levels = c("New Bedford Harbor wild", "Scorton Creek wild", "New Bedford Harbor F2", "Scorton Creek F2"))
w$data <- w$data %>% mutate(variable = recode(variable, "Observed" = "Richness", "Shannon" = "Shannon Index", "Simpson" = "Simpson Index"))

# old version with facet labels on top and wrapped facet labels
# p_alpha <- ggplot(p$data %>% filter(sampleType=="gut"), aes(x=wild_or_F2, y=value)) + geom_boxplot(aes(fill=fishType)) + scale_fill_manual(name="Fish type", values=cols.fishType[1:4], labels=labels.short)+geom_point(data=w$data, aes(x="wild", y=value, shape=sampleType), color="blue") + 
#   facet_wrap(~variable, scales="free", labeller = label_wrap_gen(width=8)) + labs(y="", x="") + theme_light() + 
#   theme(legend.position = "right", legend.text = element_text(size=10), strip.text = element_text(size=10, color="black"), strip.background = element_rect(fill="white"), axis.text.x = element_text(size=10), plot.title = element_text(size=12)) + scale_shape_manual(values=17)+guides(fill=guide_legend(order=1),shape=guide_legend(title="", order=2))

# new version with facet labels substituting for y axis labels
p_alpha <- ggplot(p$data %>% filter(sampleType=="gut"), aes(x=wild_or_F2, y=value)) + geom_boxplot(aes(fill=fishType)) + scale_fill_manual(name="Fish type", values=cols.fishType[1:4], labels=labels.short)+geom_point(data=w$data, aes(x="wild", y=value, shape=sampleType), color="blue") + 
  facet_wrap(~variable, scales="free", strip.position = "left", labeller = labeller(variable=c("Richness"="Observed ASVs", "Shannon Index" = "Shannon diversity index", "Simpson Index" = "Simpson diversity index"))) + labs(y="", x="") + theme_light() + 
  theme(legend.position = "right", legend.text = element_text(size=10), 
        strip.text = element_text(size=10, color="black"), strip.background = element_rect(fill="white"), strip.placement = "outside",
        axis.text.x = element_text(size=10), plot.title = element_text(size=12)) + scale_shape_manual(values=17)+guides(fill=guide_legend(order=1),shape=guide_legend(title="", order=2))

# adding significance bars
p_alpha <- p_alpha + 
  geom_line(data=tibble(x=c(.75,2.25), y=950, variable="Richness"), aes(x=x, y=y)) +
  geom_text(data=tibble(x=1.5, y=955, variable="Richness"), aes(x=x, y=y, label="*")) +
  geom_line(data=tibble(x=c(.75,2.25), y=5, variable="Shannon Index"), aes(x=x, y=y)) +
  geom_text(data=tibble(x=1.5, y=5.1, variable="Shannon Index"), aes(x=x, y=y, label="*")) +
  geom_line(data=tibble(x=c(.75, 2.25), y=.99, variable="Simpson Index"), aes(x=x, y=y)) +
  geom_text(data=tibble(x=1.5, y=1, variable="Simpson Index"), aes(x=x, y=y, label="*")) +
  geom_line(data=tibble(x=c(1.25, 1.75), y=.9, variable="Simpson Index"), aes(x=x, y=y)) +
  geom_text(data=tibble(x=1.5, y=.91, variable="Simpson Index"), aes(x=x, y=y, label="*"))

p_legend <- get_legend(p_alpha + 
                         theme(legend.position = "bottom", legend.margin = margin(), legend.box.margin = margin(), legend.spacing = unit(0,"cm")) + 
                         guides(fill=guide_legend(title=""), shape=guide_legend(nrow=1, order=4, title="")))

ggsave(p_alpha + 
         theme(legend.position = "bottom", legend.margin = margin(), legend.box.margin = margin(), legend.spacing = unit(0,"cm")) + 
         guides(fill=guide_legend(title=""), shape=guide_legend(nrow=1, order=4, title="")),
       filename=here("figures/plot_alpha_diversity.png"))

p_alpha <- p_alpha + theme(legend.position="none")



# perform significance tests
rich_est <- estimate_richness(ps, measures=c("Observed", "Shannon","Simpson"))
# make into a tibble with relevant columns from metadata
rich_est <- rich_est %>% rownames_to_column(var="sample") %>% left_join(as(sample_data(ps), "data.frame"), by="sample") %>% filter(sampleType=="gut") %>% select(sample, Observed, Shannon, Simpson, site, fishType) %>% pivot_longer(cols = Observed:Simpson, names_to = "measurement", values_to = "values")
# do pairwise t tests on each alpha diversity measurement type
rich_est_test <- rich_est %>% group_by(measurement) %>% nest() %>% mutate(test = map(.x=data, ~pairwise.t.test(.x$values, .x$fishType, p.adjust.method = "holm") %>% tidy())) %>% select(measurement, test) %>% unnest(cols=c(test)) %>% mutate(significant = p.value < 0.05) %>% ungroup()

## Beta diversity 

## Make zcomp
f <- zCompositions::cmultRepl(t(otu_table(ps)), method="CZM", label=0, output="p-counts") %>% t()
ps.zcomp <- ps
otu_table(ps.zcomp) <- otu_table(f, taxa_are_rows = TRUE)
ps.zcomp <- microbiome::transform(ps.zcomp, "clr")

## get dispersion of aitchinson distances
dat <- ps.zcomp %>% subset_samples(sampleType=="gut") %>% sample_data() %>% as("data.frame")
disper <- ps.zcomp %>% subset_samples(sampleType=="gut") %>% prune_taxa(taxa_sums(.)>0,.) %>% phyloseq::distance("euclidean") %>% betadisper(d=., dat$fishType)
df <- data.frame(disper=disper$distances,group=disper$group)
df$group <- factor(df$group, levels = factor.fishType)
## plot beta diversity
p_beta <- ggplot(df, aes(x=group, y=disper)) + geom_boxplot(aes(fill=group)) +
  scale_fill_manual(values=cols.fishType[1:4]) +
  scale_x_discrete(labels=labels.short) +
  labs(x="", y="Distance to spatial median") +
  theme_light()+
  theme(legend.position = "none", title=element_text(size=12))

p_beta <- p_beta + 
  geom_line(data=tibble(x=c(1,4), y=85), aes(x=x, y=y)) + 
  geom_text(data=tibble(x=2.5, y=86), aes(x=x, y=y, label="*")) +
  geom_line(data=tibble(x=c(2,3), y=75), aes(x=x, y=y)) +
  geom_text(data=tibble(x=2.5, y=76), aes(x=x, y=y, label="*"))

ggsave(p_beta, filename=here("figures/plot_beta_diversity.png"))

p_full <- plot_grid(p_alpha, p_beta, labels= c("A", "B"))/p_legend + plot_layout(heights=c(1,.1))
ggsave(p_full, filename=here("figures/plot_diversity.png"), width=7, height=5, unit="in")


