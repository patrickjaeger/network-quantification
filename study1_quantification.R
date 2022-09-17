library(tidyverse)
library(gtools)
source("geom_hpline.R")

# Import data -------------------------------------------------------------
px2um_scale <- 1.3221  # px/micron
results_path <- "data/network1_results_csv/"
additive <- "_P2CK_"

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

dat_plot <- dat_raw %>% 
  mutate(total_network = map_dbl(network,
                                 get_network_length),
         longest_network = map_dbl(network,
                                   get_network_length,
                                   .return_longest = TRUE)) %>% 
  mutate(network_count = map_int(network, get_network_count),
         mean_network_length = total_network/network_count) %>% 
  select(-network, -branch) %>% 
  mutate(norm_total_network = total_network/count,
         norm_longest_network = longest_network/count) %>% 
  mutate(day = as_factor(day)) %>% 
  # mutate(grid = ifelse(grid == "grid", "grid", "control")) %>% 
  mutate(grid = as_factor(grid) %>% fct_relevel("gel", "grid"))
dat_plot


# Plot --------------------------------------------------------------------

## Total
ggplot(dat_plot, 
       aes(day, total_network, 
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

# ggsave("results/network1_total_length.svg",
       # width = 300, height = 330, units = "px", dpi = 72)

# Per cell
ggplot(dat_plot, 
       aes(day, norm_total_network, 
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

# ggsave("results/network1_norm_length.svg",
       # width = 300, height = 330, units = "px", dpi = 72)

# Connectivity
ggplot(dat_plot, 
       aes(day, mean_network_length, 
           group = interaction(day, grid), color = grid)) +
  stat_summary(geom = "hpline", width = 0.4,
               position = position_dodge(width = 0.9), alpha = 0.8) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, 
                                              dodge.width = 0.9),
              pch = 1, show.legend = FALSE) +
  facet_wrap(~additive) +
  theme_bw() +
  labs(x = "Day",
       y = "Average network length [μm]",
       color = element_blank()) +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        strip.background=element_rect(fill="white")) +
  scale_color_manual(values = c("black", "dodgerblue"))

# ggsave("results/network1_avg_network_length.svg",
       # width = 300, height = 330, units = "px", dpi = 72)


# Numeric results ---------------------------------------------------------

dat_sample_results <- dat_plot
dat_sample_results %>% print(n = nrow(.))

# write_csv(dat_sample_results, "results/network1_results_by_sample.csv")

dat_summary <- dat_plot %>% 
  group_by(additive, day, grid) %>% 
  summarise(mean_cell_count = mean(count),
            sd_cell_count = sd(count),
            mean_total_network = mean(total_network),
            sd_total_network = sd(total_network),
            mean_norm_network = mean(norm_total_network),
            sd_norm_network = sd(norm_total_network),
            avg_network_length = mean(mean_network_length),
            sd_network_length = sd(mean_network_length))
dat_sample_results %>% print(n = nrow(.))

# write_csv(dat_summary, "results/network1_results_summary.csv")









