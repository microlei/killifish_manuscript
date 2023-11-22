# code for plotting bar graph of important families

library(here)
source(here("code/helpful_functions.R"))

ps <- readRDS(file=here("data/processed/ps.rds"))

set.seed(1)
metadata <- sample_data(ps) %>% as("data.frame")
da_wild_f2 <- readRDS(here("data/processed/da_wild_f2.rds"))
da_wild_f2_clean <- cleanDA(da_wild_f2, "wild_or_F2")

# I used the below to randomly sample 5 fish from each fish type, but for reproducibility, I manually enter the sample numbers
# ps_sub <- metadata %>% filter(sampleType=="gut") %>% group_by(fishType) %>% sample_n(5) %>% use_series(sample)

ps_sub <- c("KLF_11", "KLF_79", "KLF_102", "KLF_243", "KLF_211", "KLF_106", "KLF_42", "KLF_107", "KLF_108", "KLF_207", 
            "KLF_228", "KLF_208", "KLF_36", "KLF_101", "KLF_13", "KLF_4", "KLF_60", "KLF_6", "KLF_69", "KLF_45" )
ps_sub <- ps %>% subset_samples(sample %in% ps_sub) %>% tax_glom(taxrank="Family") %>% transform_sample_counts(~ ./sum(.)) %>% psmelt()

# get top five most significant families from differential abundance analysis
sig_family <- da_wild_f2_clean$Family %>% unique() %>% .[1:5]

# make Family_f, a factor for plotting/labelling
ps_sub <- ps_sub %>% mutate(wild_or_F2 = factor(wild_or_F2)) %>% mutate(Family_f = factor(Family, levels=c(sig_family, ps_sub$Family %>% unique() %>% .[!. %in% sig_family])))

# labels for families top five and then "others"
labels_family <- c(levels(ps_sub$Family_f)[1:5], rep("Others", times=length(unique(ps_sub$Family_f)) - 5))

# label every family thats not in top 5 "others"
ps_sub <- ps_sub %>%
  mutate(Family_f = factor(Family_f, labels = labels_family))

# assign colors
cols_family <- c(hcl.colors(5, palette = "viridis"), rep("gray", length(ps_sub$Family_f %>% unique)-5))

# make plot
p_bar <- ggplot(ps_sub, aes(x=sample, y = Abundance, fill=Family_f)) + 
  geom_bar(stat="identity") + 
  facet_wrap(~ wild_or_F2, scales="free_x") + 
  scale_fill_manual("Family", 
                    breaks = c(levels(ps_sub$Family_f)[1:5], "Others"), # this is where the labels get made!
                    values = cols_family) + 
  ylab("Relative Abundance") + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + 
  theme_light() + 
  theme(axis.text.x = element_blank(), strip.text = element_text(size=10, color="black"), strip.background = element_rect(fill="white"))

ggsave(plot=p_bar, filename=here("figures/plot_bar_f2_wild.png"), width=4.4, height=2.7, unit="in")

## Below is the code to plot the full bar plot for every sample
ps_full <- ps %>% tax_glom(taxrank="Family") %>% transform_sample_counts(~ ./sum(.)) %>% psmelt()
ps_full <- ps_full %>% mutate(wild_or_F2 = factor(wild_or_F2)) %>% mutate(Family_f = factor(Family, levels=c(sig_family, ps_full$Family %>% unique() %>% .[!. %in% sig_family])))
labels_family_full <- c(levels(ps_full$Family_f)[1:5], rep("Others", times=length(unique(ps_full$Family_f)) - 5))
ps_full <- ps_full %>% mutate(Family_f = factor(Family_f, labels = labels_family_full))
cols_family_full <- c(hcl.colors(5, palette = "viridis"), rep("gray", length(ps_full$Family_f %>% unique)-5))
p_bar_full <- ggplot(ps_full %>% filter(wild_or_F2 !="water"), aes(x=sample, y = Abundance, fill=Family_f)) + 
  geom_bar(stat="identity") + 
  facet_wrap(~ wild_or_F2, scales="free_x") + 
  scale_fill_manual("Family", 
                    breaks = c(levels(ps_full$Family_f)[1:5], "Others"), # this is where the labels get made!
                    values = cols_family_full) + 
  ylab("Relative abundance") + 
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) + 
  theme_light() + 
  theme(axis.text.x = element_blank(), strip.text = element_text(size=10, color="black"), strip.background = element_rect(fill="white"))
p_bar_full
ggsave(plot=p_bar_full, filename=here("figures/plot_bar_full.png"), width=7, height=4, unit="in")
saveRDS(p_bar_full, file="data/processed/p_bar_full.rds")
