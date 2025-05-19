# ─── Setup ───────────────────────────────────────────────────────
options(repos = c(CRAN = "https://cloud.r-project.org"))
packages <- c("ff", "ggplot2", "viridis", "patchwork")
installed <- rownames(installed.packages())
to_install <- setdiff(packages, installed)
if (length(to_install)) install.packages(to_install)

library(ff)
library(ggplot2)
library(viridis)
library(patchwork)

# ─── Load Matrix ─────────────────────────────────────────────────
binfile <- "../results/global_dp_matrix.bin"

# Read matrix dimensions
con <- file(binfile, "rb")
rows <- readBin(con, integer(), n = 1, size = 4, endian = "little")
cols <- readBin(con, integer(), n = 1, size = 4, endian = "little")
int_data <- readBin(con, integer(), n = rows * cols, size = 4, endian = "little")
close(con)

mat_ff <- ff(initdata = int_data, dim = c(rows, cols), vmode = "integer")

# ─── Define tiling parameters ────────────────────────────────────
tile_size <- 1000L
row_tiles <- ceiling(rows / tile_size)
col_tiles <- ceiling(cols / tile_size)

# ─── Generate plots tile-by-tile ─────────────────────────────────
tile_plots <- list()

for (i in 1:row_tiles) {
  for (j in 1:col_tiles) {
    r_start <- (i - 1) * tile_size + 1
    r_end <- min(i * tile_size, rows)
    c_start <- (j - 1) * tile_size + 1
    c_end <- min(j * tile_size, cols)

    sub_m <- mat_ff[r_start:r_end, c_start:c_end]

    df <- expand.grid(
      x = seq_len(c_end - c_start + 1),
      y = seq_len(r_end - r_start + 1)
    )
    df$z <- as.vector(sub_m)

    p <- ggplot(df, aes(x = x, y = y, fill = z)) +
      geom_raster(interpolate = TRUE) +
      scale_fill_viridis(option = "C", guide = "none") +
      coord_fixed(expand = FALSE) +
      theme_void() +
      theme(plot.margin = margin(0, 0, 0, 0))

    tile_plots[[length(tile_plots) + 1]] <- p
  }
}

# ─── Combine plots (patchwork grid) ─────────────────────────────
plot_grid <- wrap_plots(tile_plots, ncol = col_tiles)

# ─── Save output ────────────────────────────────────────────────
tile_width_in  <- 3  # inches per tile
tile_height_in <- 3

# Cap total size to stay within 50 inches
max_width  <- min(col_tiles * tile_width_in, 50)
max_height <- min(row_tiles * tile_height_in, 50)

ggsave("joined_heatmap.png", plot = plot_grid,
       width = max_width,
       height = max_height,
       dpi = 300, limitsize = FALSE)
