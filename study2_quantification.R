library(tidyverse)
library(gtools)
source("geom_hpline.R")

# Import data -------------------------------------------------------------
px2um_scale <- 1.3221  # px/micron
results_path <- "data/network2_noprune_results_csv/"
additive <- "_GM6001_"

branch_files <- list.files(results_path,
                           pattern = "branch-info",
                           full.names = TRUE) %>% mixedsort()
network_files <- list.files(results_path,
                            pattern = "network-info",
                            full.names = TRUE) %>% mixedsort()
nucleus_file <- list.files(results_path,
                           pattern = "nucleus_count",
                           full.names = TRUE)

dat_raw <-
  read_csv(nucleus_file) %>%
  select(Slice, Count) %>% 
  rename_with(tolower, everything()) %>% 
  rename(cell_count = count) %>% 
  mutate(slice = str_remove(slice, "C1-MAX_")) %>% 
  mutate(slice = str_replace(slice, "_\\-_", "_control_")) %>% 
  mutate(slice = str_replace(slice, "_\\+_", additive)) %>% 
  separate("slice", c("img_id", "day", "additive", "sample", "grid"), "_") %>% 
  mutate(grid = str_remove(grid, ".tif")) %>% 
  mutate(img_id = parse_number(img_id),
    day = parse_number(day) %>% as_factor()) %>% 
  arrange(img_id) %>% 
  mutate(network = network_files,
         branch = branch_files) %>% 
  mutate(network = map(network, read_csv),
         branch = map(branch, read_csv))
dat_raw


# Quantify network length -------------------------------------------------

get_network_length <- function(.df, .return_longest = FALSE) {
  # Calculate the total network length; alternatively, only
  # return the longest identified network
  # Plus, conversion of pixels to micron: scale = 1.3221 px/micron
  lengths <- .df %>% 
    select(junctions_voxels = `# Junction voxels`, 
           slab_voxels = `# Slab voxels`) %>% 
    transmute(length = (junctions_voxels + slab_voxels)/px2um_scale)
  
  if (.return_longest) max(lengths) else sum(lengths)
}

get_network_count <- function(.df) {
  # Count the number of networks
  nrow(.df)
}

get_median_network <- function(.df) {
  
  lengths <- .df %>% 
    select(junctions_voxels = `# Junction voxels`, 
           slab_voxels = `# Slab voxels`) %>% 
    transmute(length = (junctions_voxels + slab_voxels)/px2um_scale)
  
  median(lengths$length)
}

## network length quantified through branch files
dat_plot_branch <- dat_raw %>% select(-network) %>% unnest(branch) %>% 
  ungroup() %>% 
  group_by(day, additive, sample, grid, cell_count) %>% 
  summarise(n_skeletons = n(),
            total_length = sum(`Branch length`)/px2um_scale) %>% 
  mutate(length_per_cell = total_length/cell_count)
dat_plot_branch

write_csv(dat_plot_branch, "results/study2_new/network2_results_by_sample.csv")

dat_plot_branch %>% 
  group_by(additive, day, grid) %>% 
  summarise(n = n()) %>% 
  pluck("n") %>% 
  sum() 


# Plots ----

## total network ----
ggplot(dat_plot_branch, 
       aes(day, total_length, 
           group = interaction(day, grid), color = grid)) +
  stat_summary(geom = "hpline", width = 0.4,
               position = position_dodge(width = 0.9), alpha = 0.8) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, 
                                              dodge.width = 0.9),
              pch = 1, show.legend = FALSE) +
  facet_wrap(~additive) +
  theme_bw() +
  labs(x = "Day",
       y = "Total network length [μm]",
       color = element_blank()) +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        strip.background=element_rect(fill="white")) +
  scale_color_manual(values = c("black", "dodgerblue"))

ggsave("results/study2_new/network2_total_length.svg",
       width = 300, height = 330, units = "px", dpi = 72)

## network length per cell ----
ggplot(dat_plot_branch, 
       aes(day, length_per_cell, 
           group = interaction(day, grid), color = grid)) +
  stat_summary(geom = "hpline", width = 0.4,
               position = position_dodge(width = 0.9), alpha = 0.8) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, 
                                              dodge.width = 0.9),
              pch = 1, show.legend = FALSE) +
  facet_wrap(~additive) +
  theme_bw() +
  labs(x = "Day",
       y = "Network length per cell [μm]",
       color = element_blank()) +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        strip.background=element_rect(fill="white")) +
  scale_color_manual(values = c("black", "dodgerblue"))

ggsave("results/study2_new/network2_length_per_cell.svg",
       width = 300, height = 330, units = "px", dpi = 72)


## branch length histogram ----
dat_plot_branch_2 <- dat_raw %>% select(-network) %>% unnest(branch)

### all branches ----
ggplot(dat_plot_branch_2, 
       aes(`Branch length`/px2um_scale, 
           group = interaction(day, additive, grid), 
           color = grid)) +
  # geom_density() +
  geom_freqpoly(bins = 50) +
  facet_wrap(day~additive, nrow = 3) +
  theme_bw() +
  labs(x = "Branch length [μm]",
       y = "Count [n]",
       color = element_blank()) +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        strip.background=element_rect(fill="white")) +
  scale_color_manual(values = c("black", "dodgerblue"))

ggsave("results/study2_new/network2_branch_lengths_histogram.svg",
       width = 300, height = 330, units = "px", dpi = 72)


### short branches ----
ggplot(dat_plot_branch_2 %>% filter(`Branch length` < 50), 
       aes(`Branch length`/px2um_scale, 
           group = interaction(day, additive, grid), 
           color = grid)) +
  # geom_density() +
  geom_freqpoly(bins = 30) +
  facet_wrap(day~additive, nrow = 3) +
  theme_bw() +
  labs(x = "Branch length [μm]",
       y = "Count [n]",
       color = element_blank(),
       title = "Short branch distribution (<50μm)") +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        strip.background=element_rect(fill="white")) +
  scale_color_manual(values = c("black", "dodgerblue"))
  
ggsave("results/study2_new/network2_short_branch_lengths_histogram.svg",
       width = 300, height = 330, units = "px", dpi = 72)


### long branches ----
ggplot(dat_plot_branch_2 %>% filter(`Branch length` >= 50), 
       aes(`Branch length`/px2um_scale, 
           group = interaction(day, additive, grid), 
           color = grid)) +
  # geom_density() +
  geom_freqpoly(bins = 30) +
  facet_wrap(day~additive, nrow = 3) +
  theme_bw() +
  labs(x = "Branch length [μm]",
       y = "Count [n]",
       color = element_blank(),
       title = "Long branch distribution (>50μm)") +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        strip.background=element_rect(fill="white")) +
  scale_color_manual(values = c("black", "dodgerblue"))

ggsave("results/study2_new/network2_long_branch_lengths_histogram.svg",
       width = 300, height = 330, units = "px", dpi = 72)


## average branch length comparison - is shit ----
dat_plot_avg_branch <- dat_raw %>% select(-branch) %>% unnest(network) %>% 
  ungroup() %>% 
  group_by(day, additive, grid, sample) %>% 
  summarise(n_skeletons = n(),
            avg_branch_length = mean(`Average Branch Length`),
            median_branch_length = median(`Average Branch Length`))

ggplot(dat_plot_avg_branch, 
       aes(day, avg_branch_length, 
           group = interaction(day, grid), color = grid)) +
  stat_summary(geom = "hpline", width = 0.4,
               position = position_dodge(width = 0.9), alpha = 0.8) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2,
                                              dodge.width = 0.9),
              pch = 1, show.legend = FALSE) +
  facet_wrap(~additive) +
  theme_bw() +
  labs(x = "Day",
       y = "Avg. branch length [μm]",
       color = element_blank()) +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        strip.background=element_rect(fill="white")) +
  scale_color_manual(values = c("black", "dodgerblue"))

