
library(readxl)
library(dplyr)
library(stringr)
library(janitor)
library(ggplot2)
library(tidyr)

file <- "/home/laboratorio/Documentos/Ancestria/Generaltables.xls"

raw <- read_excel(file, col_names = FALSE)

header_rows <- which(str_detect(raw[[1]], "^Population \\(N[0-9]+\\)"))
header_rows

extract_block <- function(raw, start_row, end_row = NULL) {
  if (is.null(end_row)) {
    block <- raw[(start_row + 1):nrow(raw), ]
  } else {
    block <- raw[(start_row + 1):(end_row - 1), ]
  }
  
  block |>
    filter(!if_all(everything(), is.na)) |>
    janitor::remove_empty("rows")
}

# Extraer bloques
df_n5  <- extract_block(raw, header_rows[1], header_rows[2])
df_n10 <- extract_block(raw, header_rows[2], header_rows[3])
df_n11 <- extract_block(raw, header_rows[3])

# ---- CLAVE: setear nombres limpiando NA ----
set_colnames <- function(df, header_row) {
  nm <- raw[header_row, ] |> unlist(use.names = FALSE) |> as.character()
  
  # detectar nombres inválidos
  bad <- is.na(nm) | nm == ""
  
  # quitar esas columnas
  df <- df[, !bad, drop = FALSE]
  nm <- nm[!bad]
  
  colnames(df) <- nm
  df
}

df_n5  <- set_colnames(df_n5,  header_rows[1])
df_n10 <- set_colnames(df_n10, header_rows[2])
df_n11 <- set_colnames(df_n11, header_rows[3])

names(df_n5)[1]  <- "Continent"
names(df_n10)[1] <- "Population"
names(df_n11)[1] <- "Population"


# Okabe–Ito (colorblind-safe)
okabe_ito <- c(
  "Pool 1–50"    = "#0072B2",
  "Pool 51–100"  = "#009E73",
  "Merged 1–100" = "#CC79A7",
  "Hospital"     = "#56B4E9"
)

# Convertir df_n5 a formato largo
df_n5_long <- df_n5 %>%
  pivot_longer(
    cols = -Continent,
    names_to = "condition",
    values_to = "prop"
  ) %>%
  mutate(
    # Extraer cohort (Pool 1–50 / 51–100 / Hospital / Merged)
    cohort = case_when(
      str_detect(condition, "^1-50") ~ "Pool 1–50",
      str_detect(condition, "^51-100") ~ "Pool 51–100",
      str_detect(condition, "^Hospital") ~ "Hospital",
      str_detect(condition, "^Merged") ~ "Merged 1–100",
      TRUE ~ NA_character_
    ),
    # Extraer fracción (10/30/50/100) si existe; Merged no tiene % en tu tabla
    frac = case_when(
      str_detect(condition, "\\(10%\\)") ~ 0.10,
      str_detect(condition, "\\(30%\\)") ~ 0.30,
      str_detect(condition, "\\(50%\\)") ~ 0.50,
      str_detect(condition, "\\(100%\\)") ~ 1.00,
      TRUE ~ NA_real_
    ),
    cohort = factor(cohort, levels = c("Pool 1–50", "Pool 51–100", "Merged 1–100", "Hospital")),
    Continent = factor(Continent, levels = c("AFR","EUR","AMR","EAS","SAS"))
  )

p1 <- df_n5_long %>%
  filter(!is.na(frac)) %>%   # excluye Merged si no tiene (10/30/50/100)
  ggplot(aes(x = frac, y = prop, color = cohort, group = cohort)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ Continent, ncol = 3, scales = "free_y") +
  scale_color_manual(values = okabe_ito, name = "Cohort") +
  scale_x_continuous(
    breaks = c(0.10, 0.30, 0.50, 1.00),
    labels = c("10%", "30%", "50%", "100%"),
    limits = c(0.10, 1.00)
  ) +
  labs(
    title = "Stability of continental ancestry estimates under downsampling",
    x = "Read fraction (downsampling)",
    y = "Estimated ancestry proportion"
  ) +
  theme_classic(base_family = "sans", base_size = 12) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 10),
    legend.position = "top",
    legend.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 11, face = "bold"),
    panel.grid = element_blank()
  )

p1



---------------------------------
df_delta <- df_n5_long %>%
  filter(!is.na(frac), !is.na(cohort), !is.na(prop),
         cohort %in% c("Pool 1–50", "Pool 51–100", "Hospital")) %>%
  group_by(Continent, cohort) %>%
  mutate(ref_100 = prop[frac == 1.00][1]) %>%
  ungroup() %>%
  filter(!is.na(ref_100)) %>%
  mutate(delta = prop - ref_100)

p2 <- ggplot(df_delta, aes(x = frac, y = delta, color = cohort, group = cohort)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  facet_wrap(~ Continent, ncol = 3, scales = "fixed") +
  scale_color_manual(values = okabe_ito, name = "Cohort") +
  scale_x_continuous(
    breaks = c(0.10, 0.30, 0.50, 1.00),
    labels = c("10%", "30%", "50%", "100%"),
    limits = c(0.10, 1.00)
  ) +
  scale_y_continuous(
    limits = c(-ylim_abs, ylim_abs),
    labels = scales::label_number(accuracy = 0.001),
    breaks = scales::pretty_breaks(n = 5)
  ) +
  labs(
    title = "Deviation from full data estimates",
    x = "Read fraction (downsampling)",
    y = expression(Delta~"ancestry proportion")
  ) +
  theme_classic(base_family = "sans", base_size = 12) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 10),
    legend.position = "top",
    legend.background = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(size = 11, face = "bold"),
    panel.grid = element_blank()
  )

p2



##############################################################

to_long_panel <- function(df) {
  if ("pop" %in% names(df) && !"Population" %in% names(df)) {
    df <- df %>% rename(Population = pop)
  }
  
  df %>%
    pivot_longer(
      cols = -Population,
      names_to = "condition",
      values_to = "prop"
    ) %>%
    mutate(
      # ⚠️ CONVERSIÓN CRÍTICA
      prop = as.numeric(prop),
      
      cohort = case_when(
        str_detect(condition, "^1-50") ~ "Pool 1–50",
        str_detect(condition, "^51-100") ~ "Pool 51–100",
        str_detect(condition, "^Hospital") ~ "Hospital",
        TRUE ~ NA_character_
      ),
      frac = case_when(
        str_detect(condition, "\\(10%\\)") ~ 0.10,
        str_detect(condition, "\\(30%\\)") ~ 0.30,
        str_detect(condition, "\\(50%\\)") ~ 0.50,
        str_detect(condition, "\\(100%\\)") ~ 1.00,
        TRUE ~ NA_real_
      ),
      cohort = factor(cohort, levels = c("Pool 1–50", "Pool 51–100", "Hospital"))
    ) %>%
    filter(!is.na(frac), !is.na(cohort), !is.na(prop))
}


###############
## N10
###############

df10_long <- to_long_panel(df_n10)

pop_order10 <- df10_long %>%
  group_by(Population) %>%
  summarise(mean_prop = mean(prop), .groups = "drop") %>%
  arrange(desc(mean_prop)) %>%
  pull(Population)

df10_long <- df10_long %>%
  mutate(Population = factor(Population, levels = pop_order10))


# ---- p1_N10: curvas ----
p1_N10 <- ggplot(df10_long, aes(x = frac, y = prop, color = cohort, group = cohort)) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  facet_wrap(~ Population, ncol = 5, scales = "free_y") +
  scale_color_manual(values = okabe_ito, name = "Cohort") +
  scale_x_continuous(
    breaks = c(0.10, 0.30, 0.50, 1.00),
    labels = c("10%", "30%", "50%", "100%"),
    limits = c(0.10, 1.00)
  ) +
  labs(
    title = "Downsampling stability (N10)",
    subtitle = "Y-axis scales vary by population to emphasize within-population stability",
    x = "Read fraction (downsampling)",
    y = "Estimated ancestry proportion"
  ) +
  theme_classic(base_family = "sans", base_size = 12) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 9),
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(size = 9, face = "bold"),
    panel.grid = element_blank()
  )

p1_N10

# ---- p2_N10: Δ vs 100% ----
df10_delta <- df10_long %>%
  group_by(Population, cohort) %>%
  mutate(ref_100 = prop[frac == 1.00][1]) %>%
  ungroup() %>%
  mutate(delta = prop - ref_100)

p2_N10 <- ggplot(df10_delta, aes(x = frac, y = delta, color = cohort, group = cohort)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  facet_wrap(~ Population, ncol = 5, scales = "free_y") +
  scale_color_manual(values = okabe_ito, name = "Cohort") +
  scale_x_continuous(
    breaks = c(0.10, 0.30, 0.50, 1.00),
    labels = c("10%", "30%", "50%", "100%"),
    limits = c(0.10, 1.00)
  ) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.001)) +
  labs(
    title = "Deviation from full data estimates",
    x = "Read fraction (downsampling)",
    y = expression(Delta~"ancestry proportion")
  ) +
  theme_classic(base_family = "sans", base_size = 12) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 9),
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(size = 9, face = "bold"),
    panel.grid = element_blank()
  )

p2_N10



##N11

df11_long <- to_long_panel(df_n11)

pop_order11 <- df11_long %>%
  group_by(Population) %>%
  summarise(mean_prop = mean(prop, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_prop)) %>%
  pull(Population)

df11_long <- df11_long %>%
  mutate(Population = factor(Population, levels = pop_order11))

# ---- p1_N11: curvas ----
p1_N11 <- ggplot(df11_long, aes(x = frac, y = prop, color = cohort, group = cohort)) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  facet_wrap(~ Population, ncol = 5, scales = "free_y") +
  scale_color_manual(values = okabe_ito, name = "Cohort") +
  scale_x_continuous(
    breaks = c(0.10, 0.30, 0.50, 1.00),
    labels = c("10%", "30%", "50%", "100%"),
    limits = c(0.10, 1.00)
  ) +
  labs(
    title = "Downsampling stability (N11)",
    x = "Read fraction (downsampling)",
    y = "Estimated ancestry proportion"
  ) +
  theme_classic(base_family = "sans", base_size = 12) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 9),
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(size = 9, face = "bold"),
    panel.grid = element_blank()
  )

p1_N11

# ---- p2_N11: Δ vs 100% ----
df11_delta <- df11_long %>%
  group_by(Population, cohort) %>%
  mutate(ref_100 = prop[frac == 1.00][1]) %>%
  ungroup() %>%
  mutate(delta = prop - ref_100)

p2_N11 <- ggplot(df11_delta, aes(x = frac, y = delta, color = cohort, group = cohort)) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_line(linewidth = 0.8, na.rm = TRUE) +
  geom_point(size = 2, na.rm = TRUE) +
  facet_wrap(~ Population, ncol = 5, scales = "free_y") +
  scale_color_manual(values = okabe_ito, name = "Cohort") +
  scale_x_continuous(
    breaks = c(0.10, 0.30, 0.50, 1.00),
    labels = c("10%", "30%", "50%", "100%"),
    limits = c(0.10, 1.00)
  ) +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.001)) +
  labs(
    title = "Deviation from full data estimates",
    x = "Read fraction (downsampling)",
    y = expression(Delta~"ancestry proportion")
  ) +
  theme_classic(base_family = "sans", base_size = 12) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 9),
    legend.position = "top",
    strip.background = element_blank(),
    strip.text = element_text(size = 9, face = "bold"),
    panel.grid = element_blank()
  )

p2_N11

