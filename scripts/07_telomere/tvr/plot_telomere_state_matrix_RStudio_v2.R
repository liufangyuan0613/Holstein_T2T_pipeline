# ==========================================
# plot_telomere_state_matrix_RStudio_v2.R
# ==========================================

pkgs <- c("readr", "dplyr", "tidyr", "ggplot2", "stringr", "ggtext")
to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
if (length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

state_matrix_file <- "assembly_telomere_v2.state_matrix.tsv"
output_prefix <- "assembly_telomere_v2_TVR_map"

tile_height <- 0.9

state_colors <- c(
  "_" = "#D9D9D9",
  "A" = "#D8A028",
  "C" = "#2F36FF",
  "D" = "#E64B35",
  "F" = "#69BE45",
  "H" = "#D39BEA",
  "L" = "#E5D84C",
  "N" = "#B5B6FF",
  "S" = "#9FD3C7",
  "T" = "#8E6BD8",
  "V" = "#8AC0E8",
  "O" = "#7F7F7F"
)

end_colors <- c("5'" = "#228B22", "3'" = "#F28C00")

chr_order_key <- function(x) {
  x <- as.character(x)
  core <- sub("^chr", "", x)
  num <- suppressWarnings(as.numeric(core))
  ifelse(!is.na(num), num, ifelse(core == "X", 1000, ifelse(core == "Y", 1001, 2000)))
}

df <- readr::read_tsv(state_matrix_file, show_col_types = FALSE)
bin_cols <- colnames(df)[stringr::str_detect(colnames(df), "^bin_[0-9]+$")]
if (length(bin_cols) == 0) stop("No bin_* columns found in state_matrix.tsv")

row_df <- df %>%
  mutate(chr_order = chr_order_key(chrom), end_order = ifelse(end_label == "5'", 1, 2)) %>%
  arrange(chr_order, end_order)
row_levels <- row_df$chr_end

long_df <- row_df %>%
  select(chrom, end_label, chr_end, telomere_length_bp, all_of(bin_cols)) %>%
  pivot_longer(cols = all_of(bin_cols), names_to = "bin", values_to = "state") %>%
  mutate(bp = as.numeric(stringr::str_remove(bin, "^bin_")), chr_end = factor(chr_end, levels = rev(row_levels)))

bin_starts <- sort(as.numeric(stringr::str_remove(bin_cols, "^bin_")))
bin_bp <- if (length(bin_starts) >= 2) min(diff(bin_starts)) else 20
tile_width <- bin_bp

label_df <- row_df %>%
  mutate(label_text = case_when(
    end_label == "5'" ~ paste0(chrom, " ", "<span style='color:", end_colors["5'"], ";font-weight:600;'>5'</span>"),
    end_label == "3'" ~ paste0(chrom, " ", "<span style='color:", end_colors["3'"], ";font-weight:600;'>3'</span>")
  ))
label_map <- label_df$label_text
names(label_map) <- label_df$chr_end
label_map <- label_map[rev(names(label_map))]

p <- ggplot(long_df, aes(x = bp, y = chr_end, fill = state)) +
  geom_tile(width = tile_width, height = tile_height) +
  scale_fill_manual(values = state_colors, drop = FALSE) +
  scale_y_discrete(labels = label_map) +
  scale_x_continuous(expand = c(0, 0), breaks = pretty(long_df$bp, n = 6)) +
  labs(title = "Assembly-based chromosome-end TVR map",
       x = "Distance from subtelomere boundary (bp)",
       y = NULL, fill = "state") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        axis.text.y = ggtext::element_markdown(size = 9),
        axis.text.x = element_text(size = 10),
        panel.grid = element_blank(),
        legend.position = "right")

print(p)
ggsave(paste0(output_prefix, ".png"), p, width = 10, height = 13, dpi = 300)
ggsave(paste0(output_prefix, ".pdf"), p, width = 10, height = 13)
