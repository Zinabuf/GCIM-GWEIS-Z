### Load libraries ###
library(ggplot2)
library(dplyr)
library(data.table)

#---- Read data ----#
gcim_adj  <- fread("trait1_out_q_gcim_gweis-z.txt")
gcim_raw  <- fread("trait1_out_q_gcim_gweis.txt")
gweis_adj <- fread("trait1_out_q_gweis-z.txt")
gweis_raw <- fread("trait1_out_q_gweis.txt")

#---- Keep only interaction test from raw PLINK outputs ----#
gcim_raw  <- gcim_raw[TEST == "ADDxPRS", ]
gweis_raw <- gweis_raw[TEST == "ADDxtrait2", ]

#---- Standardize columns across files ----#
to_uniform <- function(df, model_label, pcol) {

  # Find chromosome column name across file types
  chr_col <- dplyr::case_when(
    "#CHROM" %in% names(df) ~ "#CHROM",
    "X.CHROM" %in% names(df) ~ "X.CHROM",
    "CHROM" %in% names(df) ~ "CHROM",
    TRUE ~ NA_character_
  )

  if (is.na(chr_col)) {
    stop("No chromosome column found. Expected one of: #CHROM, X.CHROM, CHROM")
  }

  if (!("POS" %in% names(df))) stop("POS column not found.")
  if (!(pcol %in% names(df))) stop(paste0("P column '", pcol, "' not found."))

  df %>%
    mutate(
      CHR = as.character(.data[[chr_col]]),
      BP  = as.numeric(POS),
      P   = as.numeric(.data[[pcol]]),
      MODEL = model_label
    ) %>%
    select(CHR, BP, P, MODEL) %>%
    filter(!is.na(CHR), !is.na(BP), !is.na(P), is.finite(P), P > 0, P <= 1)
}

gcim_adj_u  <- to_uniform(gcim_adj,  "GCIM-GWEIS-Z", "p_value_int_adj")
gcim_raw_u  <- to_uniform(gcim_raw,  "GCIM-GWEIS",   "P")
gweis_adj_u <- to_uniform(gweis_adj, "GWEIS-Z",      "p_value_int_adj")
gweis_raw_u <- to_uniform(gweis_raw, "GWEIS",        "P")

#---- Combine all ----#
all_data <- bind_rows(gcim_adj_u, gcim_raw_u, gweis_adj_u, gweis_raw_u)

# Convert CHR to ordered factor (1:22 + X/Y if present)
chr_order <- as.character(c(1:22, "X", "Y"))
all_data$CHR <- factor(all_data$CHR, levels = chr_order)

# Drop chromosomes not in chr_order (e.g., "0", "MT", etc. if they appear)
all_data <- all_data %>% filter(!is.na(CHR))

#---- Compute cumulative BP for plotting ----#
chr_offsets <- all_data %>%
  group_by(CHR) %>%
  summarize(chr_len = max(BP, na.rm = TRUE), .groups = "drop") %>%
  arrange(CHR) %>%
  mutate(chr_start = lag(cumsum(chr_len), default = 0)) %>%
  select(CHR, chr_start)

all_data <- all_data %>%
  left_join(chr_offsets, by = "CHR") %>%
  mutate(BPcum = BP + chr_start)

axis_df <- all_data %>%
  group_by(CHR) %>%
  summarize(center = mean(BPcum, na.rm = TRUE), .groups = "drop")

#---- Plot ----#
p <- ggplot(all_data, aes(x = BPcum, y = -log10(P), color = CHR)) +
  geom_point(alpha = 0.7, size = 0.8) +
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", linewidth = 0.8) +
  scale_x_continuous(labels = axis_df$CHR, breaks = axis_df$center) +
  facet_wrap(~MODEL, ncol = 1, scales = "free_y") +
  labs(
    x = "Chromosome",
    y = "-log10(p)",
    title = "Manhattan Plots for GCIM-GWEIS and GWEIS Models",
    subtitle = "Dashed line = genome-wide significance threshold (5 x 10^-8)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    strip.background = element_rect(fill = "lightgray", color = NA),
    strip.text = element_text(size = 12, face = "bold")
  )

ggsave("Combined_Manhattan_trait1.png", p, width = 10, height = 12, dpi = 600)