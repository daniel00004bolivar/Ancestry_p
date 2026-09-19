# ============================
# iAdmix comparisons (Q1-ready)
# Pools: 1–50, 51–100, Merged 1–100, Hospital
# Panels: N10 vs N11 (by population), N5 (by continent)
# Figures:
#  Fig 1: N10 vs N11 by population, faceted by pool
#  Fig 2: N5 by continent, pools side-by-side
#  Fig 3: Delta (Pool - reference) by population, for each panel (stability)
#  Optional: Heatmap of proportions
# ============================

library(dplyr)
library(patchwork)
library(ggplot2)
library(stringr)
library(tidyr)

# ----------------------------
# 1) Input data (UPDATED)
# ----------------------------
df_pop <- tibble::tribble(
  ~pop, ~`N10-1-50`, ~`N11-1-50`, ~`N10-51-100`, ~`N11-51-100`, ~`N10-Merged-1-100`, ~`N11-Merged-1-100`, ~`N10-Hospital`, ~`N11-Hospital`,
  "IBS", 0.0371, 0.0361, 0.0372, 0.0359, 0.0362, 0.0357, 0.0361, 0.0313,
  "CEU", 0.0486, 0.0528, 0.0488, 0.0530, 0.0487, 0.0535, 0.0503, 0.0518,
  "YRI", 0.1323, 0.1222, 0.1322, 0.1222, 0.1322, 0.1219, 0.1301, 0.1222,
  "ESN", 0.1353, 0.1238, 0.1350, 0.1244, 0.1359, 0.1242, 0.1342, 0.1256,
  "GWD", 0.2397, 0.2220, 0.2394, 0.2225, 0.2407, 0.2225, 0.2406, 0.2273,
  "LWK", 0.2950, 0.2733, 0.2947, 0.2739, 0.2956, 0.2732, 0.2944, 0.2764,
  "PEL", 0.0181, 0.0223, 0.0181, 0.0218, 0.0177, 0.0220, 0.0176, 0.0217,
  "MXL", 0.0278, 0.0324, 0.0282, 0.0322, 0.0277, 0.0325, 0.0290, 0.0330,
  "PUR", 0.0422, 0.0452, 0.0424, 0.0452, 0.0418, 0.0451, 0.0435, 0.0443,
  "FIN", 0.0239, 0.0289, 0.0239, 0.0286, 0.0236, 0.0285, 0.0242, 0.0279,
  "CLM", NA,     0.0410, NA,     0.0402, NA,     0.0410, NA,     0.0385
)

df_cont <- tibble::tribble(
  ~Continent, ~`N5-1-50`, ~`N5-51-100`, ~`N5-Merged-1-100`, ~`N5-Hospital`,
  "AFR", 0.7265, 0.7251, 0.7296, 0.7263,
  "EUR", 0.0434, 0.0436, 0.0432, 0.0450,
  "AMR", 0.0434, 0.0407, 0.0395, 0.0415,
  "EAS", 0.0747, 0.0751, 0.0738, 0.0731,
  "SAS", 0.1150, 0.1154, 0.1139, 0.1141
)

# ----------------------------
# 2) Helpers: parsing column names
# ----------------------------
parse_panel_pool <- function(x) {
  # expected patterns:
  # N10-1-50, N11-51-100, N10-Merged-1-100, N11-Hospital
  panel <- ifelse(str_detect(x, "^N10"), "N10", "N11")
  pool  <- case_when(
    str_detect(x, "1-50$") ~ "Pool 1–50",
    str_detect(x, "51-100$") ~ "Pool 51–100",
    str_detect(x, "Merged-1-100$") ~ "Merged 1–100",
    str_detect(x, "Hospital$") ~ "Hospital",
    TRUE ~ NA_character_
  )
  tibble(panel = panel, pool = pool)
}

pool_levels <- c("Pool 1–50", "Pool 51–100", "Merged 1–100", "Hospital")

# ----------------------------
# 3) Long format: N10/N11 by population across pools
# ----------------------------

df_pop_long <- df_pop |>
  pivot_longer(
    cols = -pop,
    names_to = "panel_pool",
    values_to = "prop"
  ) |>
  mutate(
    panel = if_else(str_detect(panel_pool, "^N10"), "N10", "N11"),
    pool = case_when(
      str_detect(panel_pool, "1-50$") ~ "Pool 1–50",
      str_detect(panel_pool, "51-100$") ~ "Pool 51–100",
      str_detect(panel_pool, "Merged-1-100$") ~ "Merged 1–100",
      str_detect(panel_pool, "Hospital$") ~ "Hospital",
      TRUE ~ NA_character_
    ),
    panel = factor(panel, levels = c("N10", "N11")),
    pool  = factor(pool, levels = pool_levels)
  ) |>
  select(pop, panel, pool, prop)


# Order pops by overall mean (descending)
pop_order <- df_pop_long |>
  group_by(pop) |>
  summarise(mean_prop = mean(prop, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(mean_prop)) |>
  pull(pop)

df_pop_long <- df_pop_long |>
  mutate(pop = factor(pop, levels = pop_order))

# ----------------------------
# 4) Long format: N5 by continent across pools
# ----------------------------
df_cont_long <- df_cont |>
  pivot_longer(cols = starts_with("N5-"), names_to = "pool", values_to = "prop") |>
  mutate(
    pool = recode(pool,
                  "N5-1-50" = "Pool 1–50",
                  "N5-51-100" = "Pool 51–100",
                  "N5-Merged-1-100" = "Merged 1–100",
                  "N5-Hospital" = "Hospital"),
    pool = factor(pool, levels = pool_levels),
    Continent = factor(Continent, levels = c("AFR", "EUR", "AMR", "EAS", "SAS"))
  )

# ----------------------------
# 5) Palette (Okabe–Ito)
# ----------------------------
okabe_ito_panel <- c("N10" = "#0072B2", "N11" = "#E69F00") # blue, orange
okabe_ito_pool  <- c(
  "Pool 1–50"    = "#0072B2", # blue
  "Pool 51–100"  = "#009E73", # bluish green
  "Merged 1–100" = "#CC79A7", # reddish purple
  "Hospital"     = "#56B4E9"  # sky blue
)

# ----------------------------
# 6) FIG 1: N10 vs N11 by population, faceted by pool
# ----------------------------
p1 <- ggplot(df_pop_long, aes(x = pop, y = prop, fill = panel)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62, na.rm = TRUE) +
  facet_wrap(~ pool, ncol = 2) +
  scale_fill_manual(values = okabe_ito_panel, name = "Reference Panel") +
  scale_y_continuous(
    limits = c(0, 0.32),
    breaks = seq(0, 0.32, 0.04),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Ancestry proportions by reference population across pools",
    x = "Reference population",
    y = "Estimated ancestry proportion"
  ) +
  theme_classic(base_family = "sans", base_size = 12) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10),
    legend.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 11, face = "bold"),
    panel.grid = element_blank()
  )

# ============================
# FIG: Merged vs Hospital (Barplots)
# Panels: N10, N11
# ============================
df_mh <- df_pop_long |>
  filter(pool %in% c("Merged 1–100", "Hospital")) |>
  mutate(
    pool = factor(pool, levels = c("Merged 1–100", "Hospital"))
  )

# Reordenar poblaciones por el valor promedio (Merged + Hospital)
pop_order_mh <- df_mh |>
  group_by(pop) |>
  summarise(mean_prop = mean(prop, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(mean_prop)) |>
  pull(pop)

df_mh <- df_mh |>
  mutate(pop = factor(pop, levels = pop_order_mh))

# Paleta Okabe–Ito
cols_pool <- c(
  "Merged 1–100" = "#0072B2",  # blue
  "Hospital"     = "#009E73"   # bluish green
)

p_mh <- ggplot(df_mh, aes(x = pop, y = prop, fill = pool)) +
  geom_col(
    position = position_dodge(width = 0.72),
    width = 0.62,
    na.rm = TRUE
  ) +
  facet_wrap(~ panel, ncol = 1) +
  scale_fill_manual(values = cols_pool, name = "Dataset") +
  scale_y_continuous(
    limits = c(0, 0.32),
    breaks = seq(0, 0.32, 0.04),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Comparison of ancestry proportions: Uninorte vs Hospital cohort",
    x = "Reference population",
    y = "Estimated ancestry proportion"
  ) +
  theme_classic(base_family = "sans", base_size = 12) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10),
    legend.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 11, face = "bold"),
    panel.grid = element_blank()
  )

p_mh


# ----------------------------
# 7) FIG 2: N5 by continent across pools (pool as color)
# ----------------------------
p2 <- ggplot(df_cont_long, aes(x = Continent, y = prop, fill = pool)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  scale_fill_manual(values = okabe_ito_pool, name = "Pool") +
  scale_y_continuous(
    limits = c(0, 0.85),
    breaks = seq(0, 0.85, 0.1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Aancestry proportions across pools",
    x = "Reference continent",
    y = "Estimated ancestry proportion"
  ) +
  theme_classic(base_family = "sans", base_size = 12) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 10),
    legend.position = "top",
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10),
    legend.background = element_blank(),
    panel.grid = element_blank()
  )


# ============================
# N5: Merged 1–100 vs Hospital
# ============================

df_mh_cont <- df_cont_long |>
  filter(pool %in% c("Merged 1–100", "Hospital")) |>
  mutate(pool = factor(pool, levels = c("Merged 1–100", "Hospital")))

# Okabe–Ito colors (consistent, high-contrast)
cols_mh <- c(
  "Merged 1–100" = "#CC79A7",  # reddish purple
  "Hospital"     = "#56B4E9"   # sky blue
)

p_mh_cont <- ggplot(df_mh_cont, aes(x = Continent, y = prop, fill = pool)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62) +
  scale_fill_manual(values = cols_mh, name = "Pool") +
  scale_y_continuous(
    limits = c(0, 0.85),
    breaks = seq(0, 0.85, 0.1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    title = "Ancestry proportions: Uninorte vs Hospital",
    x = "Reference continent",
    y = "Estimated ancestry proportion"
  ) +
  theme_classic(base_family = "sans", base_size = 12) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 10),
    legend.position = "top",
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10),
    legend.background = element_blank(),
    panel.grid = element_blank()
  )

p_mh_cont



# ----------------------------
# 10) Assemble main figure (recommended)
# ----------------------------
# Main (double-column): Fig 1 + Fig 2
p_main <- p_mh + p_mh_cont 
p_main


# ----------------------------
# 12) “Parameters you can change”
# ----------------------------
# - pool_levels (order of pools)
# - ref_pool (reference for deltas)
# - okabe_ito_panel / okabe_ito_pool (palettes)
# - y-axis limits in p1/p2
# - facet layout (ncol) in p1
# - export sizes in ggsave


