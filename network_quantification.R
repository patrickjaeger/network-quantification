library(tidyverse)
source("geom_hpline.R")

# Import data -------------------------------------------------------------

branch_files <- list.files("data/network1_results_csv/",
                           pattern = "branch-info",
                           full.names = TRUE)
network_files <- list.files("data/network1_results_csv/",
                            pattern = "network-info",
                            full.names = TRUE)
nucleus_file <- list.files("data/network1_results_csv/",
                           pattern = "nucleus_count",
                           full.names = TRUE)

dat_raw <-
  read_csv(nucleus_file) %>%
  select(Slice, Count) %>% 
  rename_with(tolower, everything()) %>% 
  mutate(slice = str_remove(slice, "C1-MAX_")) %>% 
  mutate(img_id = str_sub(slice, 1, 3) %>% parse_integer()) %>% 
  separate("slice", c("trash", "label"), " - ") %>% 
  select(img_id, label, count) %>% 
  mutate(label = str_remove(label, ".tif")) %>% 
  mutate(label = case_when(
    label == "d14-_1_d14-_grid" ~ "d14-_1_d14-_1 grid",
    label == "d14-_1_d14-_gel" ~ "d14-_1_d14-_1 gel",
    label == "d14-_2_d14-_ 2 gel" ~ "d14-_2_d14-_2 gel",
    TRUE ~ label)) %>% 
  mutate(label = case_when(
    img_id %in% 1:10 ~ str_sub(label, 8, str_length(label)),
    img_id %in% 11:16 ~ str_sub(label, 14, str_length(label)),
    TRUE ~ "error")) %>% 
  separate("label", c("lab1", "lab2"), "_") %>% 
  mutate(day = parse_number(lab1) %>% as.integer(),
         P2CK = ifelse(str_detect(lab1, "\\+"), "+P2CK", "-P2CK")) %>% 
  separate("lab2", c("sample", "grid"), " ") %>% 
  select(img_id, day, sample, P2CK, grid, count) %>% 
  mutate(network = map(network_files, read_csv),
         branch = map(branch_files, read_csv))
dat_raw


# Quantify network length -------------------------------------------------

get_network_length <- function(.df, .return_longest = FALSE) {
  # Calculate the total network length; alternatively, only
  # return the longest identified network
  lengths <- .df %>% 
    select(junctions_voxels = `# Junction voxels`, 
           slab_voxels = `# Slab voxels`) %>% 
    transmute(length = junctions_voxels + slab_voxels)
  
  if (.return_longest) max(lengths) else sum(lengths)
}

dat_plot <- dat_raw %>% 
  mutate(total_network = map_dbl(network,
                                 get_network_length),
         longest_network = map_dbl(network,
                                   get_network_length,
                                   .return_longest = TRUE)) %>% 
  select(-network, -branch) %>% 
  mutate(norm_total_network = total_network/count,
         norm_longest_network = longest_network/count) %>% 
  mutate(day = as_factor(day)) %>% 
  mutate(grid = ifelse(grid == "grid", "grid", "control")) %>% 
  mutate(grid = as_factor(grid) %>% fct_relevel("control", "grid"))
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
  facet_wrap(~P2CK) +
  theme_bw() +
  labs(x = "Day",
       y = "Total network length [a.u.]",
       color = element_blank()) +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        strip.background=element_rect(fill="white")) +
  scale_color_manual(values = c("black", "dodgerblue"))

# ggsave("results/network1_total_length.svg",
#        width = 300, height = 330, units = "px", dpi = 72)

# Per cell
ggplot(dat_plot, 
       aes(day, norm_total_network, 
           group = interaction(day, grid), color = grid)) +
  stat_summary(geom = "hpline", width = 0.4,
               position = position_dodge(width = 0.9), alpha = 0.8) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, 
                                              dodge.width = 0.9),
              pch = 1, show.legend = FALSE) +
  facet_wrap(~P2CK) +
  theme_bw() +
  labs(x = "Day",
       y = "Network length per cell [a.u.]",
       color = element_blank()) +
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        strip.background=element_rect(fill="white")) +
  scale_color_manual(values = c("black", "dodgerblue")) 

# ggsave("results/network1_norm_length.svg",
#        width = 300, height = 330, units = "px", dpi = 72)


# Numeric results ---------------------------------------------------------

dat_summary <- dat_plot %>% 
  group_by(P2CK, day, grid) %>% 
  summarise(mean_cell_count = mean(count),
            sd_cell_count = sd(count),
            mean_total_network = mean(total_network),
            sd_total_network = sd(total_network),
            mean_norm_network = mean(norm_total_network),
            sd_norm_network = sd(norm_total_network))
dat_summary

# write_csv(dat_summary, "results/network1_results.csv")
