###
## Load libraries
###
suppressPackageStartupMessages({
  library(tidyverse)
  library(rmetalog)
  library(MASS)       # mvrnorm
  library(ggplot2)
  library(R.utils)
})

###
## Setup directory structure
###
home_dir  <- "TrueReplica/"
data_dir  <- "TrueReplica/data/"
plots_dir <- "TrueReplica/plots/"

if (!dir.exists(home_dir))  dir.create(home_dir)
if (!dir.exists(data_dir))  dir.create(data_dir)
if (!dir.exists(plots_dir)) dir.create(plots_dir)

theme_set(theme_minimal(base_size = 16))
primary <- "#16a34a"

###
## Helper functions
###

# Safe clamp
clamp <- function(x, lo, hi) {
  pmin(pmax(x, lo), hi)
}

# Pseudo-observations: rank -> (0,1)
to_pobs <- function(x) {
  r <- rank(x, ties.method = "average", na.last = "keep")
  r / (sum(!is.na(r)) + 1)
}

# Force correlation matrix to be positive definite if needed
make_pd <- function(S, eps = 1e-6) {
  eig <- eigen(S, symmetric = TRUE)
  vals <- eig$values
  vals[vals < eps] <- eps
  S_pd <- eig$vectors %*% diag(vals) %*% t(eig$vectors)
  D <- diag(1 / sqrt(diag(S_pd)))
  D %*% S_pd %*% D
}

# Build empirical category mapping info from a factor
build_factor_info <- function(x) {
  tab <- prop.table(table(x))
  cum_p <- cumsum(tab)
  lower <- c(0, head(cum_p, -1))
  
  list(
    levels = names(tab),
    probs  = as.numeric(tab),
    cum_p  = as.numeric(cum_p),
    lower  = as.numeric(lower)
  )
}

# Map observed factor values to a representative uniform value in that category bin
factor_to_uniform_midpoint <- function(x, info) {
  idx <- match(as.character(x), info$levels)
  (info$lower[idx] + info$cum_p[idx]) / 2
}

# Map synthetic uniforms back to factor levels by empirical CDF bins
uniform_to_factor <- function(u, info, original_levels = NULL) {
  idx <- findInterval(u, vec = info$cum_p) + 1
  idx <- clamp(idx, 1, length(info$levels))
  vals <- info$levels[idx]
  
  if (is.null(original_levels)) {
    factor(vals)
  } else {
    factor(vals, levels = original_levels)
  }
}

# Try inverse metalog safely
inverse_metalog <- function(m, u, term) {
  u <- clamp(u, 1e-6, 1 - 1e-6)
  out <- qmetalog(m, y = u, term = term)
  as.numeric(out)
}

# Numeric plotting helper
plot_numeric_compare <- function(original, synthetic, var_name, out_file) {
  df <- bind_rows(
    data.frame(value = original,  Type = "Original"),
    data.frame(value = synthetic, Type = "Synthetic")
  )
  
  plt <- ggplot(df, aes(x = value, color = Type, linetype = Type)) +
    geom_density(linewidth = 1) +
    labs(
      x = var_name,
      y = "Density",
      title = paste0(var_name, ": Original vs Synthetic")
    )
  
  ggsave(out_file, plot = plt, width = 11, height = 7.5, units = "in", dpi = 300)
}

# Factor plotting helper
plot_factor_compare <- function(original, synthetic, var_name, out_file) {
  df <- bind_rows(
    data.frame(value = original,  Type = "Original"),
    data.frame(value = synthetic, Type = "Synthetic")
  )
  
  plt <- ggplot(df, aes(x = value, fill = Type)) +
    geom_bar(position = "dodge") +
    labs(
      x = var_name,
      y = "Count",
      title = paste0(var_name, ": Original vs Synthetic")
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(out_file, plot = plt, width = 11, height = 7.5, units = "in", dpi = 300)
}

###
## Read source file
###
source_file_full <- file.choose()
# source_file_full <- "TrueReplica/data/adult.csv"

source_file_name <- sub("(.*)\\..*$", "\\1", basename(source_file_full))

plots_sub_dir <- file.path(plots_dir, source_file_name)
if (!dir.exists(plots_sub_dir)) dir.create(plots_sub_dir)

source_file_data <- read.csv(source_file_full, stringsAsFactors = TRUE)

# Optional: convert character columns to factor
source_file_data <- source_file_data %>%
  mutate(across(where(is.character), as.factor))

header <- names(source_file_data)
n <- nrow(source_file_data)

cat("\nStatus: Dataset loaded with ", n, " rows and ", ncol(source_file_data), " columns.", sep = "")

###
## Identify column types
###
is_factor_col  <- sapply(source_file_data, is.factor)
is_numeric_col <- sapply(source_file_data, is.numeric)

if (!all(is_factor_col | is_numeric_col)) {
  stop("All columns must be either numeric or factor after preprocessing.")
}

###
## Fit marginals and create U matrix
###
metalog_max_terms <- 9

metalog_models <- list()
factor_infos    <- list()
integer_cols    <- logical(length(header))
names(integer_cols) <- header

U <- matrix(NA_real_, nrow = n, ncol = length(header))
colnames(U) <- header

cat("\nStatus: Fitting marginals and building pseudo-observations...")

for (j in seq_along(header)) {
  idx <- header[j]
  x <- source_file_data[[idx]]
  
  if (is.factor(x)) {
    info <- build_factor_info(x)
    factor_infos[[idx]] <- info
    U[, j] <- factor_to_uniform_midpoint(x, info)
    
    cat("\nStatus: Processed column ", j, " : ", idx, " (factor)", sep = "")
  } else {
    # numeric column
    integer_cols[idx] <- all(!is.na(x) & abs(x - round(x)) < 1e-8)
    
    m <- metalog(
      x = as.numeric(x),
      term_limit = metalog_max_terms
    )
    
    metalog_models[[idx]] <- m
    U[, j] <- to_pobs(x)
    
    cat("\nStatus: Fitted metalog for column ", j, " : ", idx, " (numeric)", sep = "")
  }
}

###
## Fit Gaussian copula using transformed normal scores
###
cat("\nStatus: Fitting Gaussian copula...")

keep <- complete.cases(U)
U_fit <- U[keep, , drop = FALSE]

if (nrow(U_fit) < 5) {
  stop("Not enough complete cases to fit copula.")
}

Z <- qnorm(clamp(U_fit, 1e-6, 1 - 1e-6))
R <- cor(Z)
R <- make_pd(R)

cat("\nStatus: Gaussian copula correlation matrix estimated.")

###
## Generate correlated uniform samples
###
cat("\nStatus: Generating correlated uniforms...")

Z_syn <- MASS::mvrnorm(
  n = n,
  mu = rep(0, ncol(U_fit)),
  Sigma = R
)

U_syn <- pnorm(Z_syn)
colnames(U_syn) <- colnames(U_fit)

###
## Invert marginals to synthetic data
###
cat("\nStatus: Inverting marginals to synthetic data...")

synthetic_df <- data.frame(matrix(nrow = n, ncol = 0))

for (j in seq_along(header)) {
  idx <- header[j]
  x <- source_file_data[[idx]]
  u <- U_syn[, idx]
  
  if (is.factor(x)) {
    info <- factor_infos[[idx]]
    synthetic_df[[idx]] <- uniform_to_factor(
      u = u,
      info = info,
      original_levels = levels(x)
    )
    
    cat("\nStatus: Generated synthetic column ", j, " : ", idx, " (factor)", sep = "")
  } else {
    m <- metalog_models[[idx]]
    syn_num <- inverse_metalog(m, u, metalog_max_terms)
    
    # If original numeric column looked integer-valued, round synthetic output
    if (integer_cols[idx]) {
      syn_num <- round(syn_num)
    }
    
    synthetic_df[[idx]] <- syn_num
    
    cat("\nStatus: Generated synthetic column ", j, " : ", idx, " (numeric)", sep = "")
  }
}

###
## Save synthetic CSV
###
cat("\nStatus: Writing synthetic CSV file...")

output_csv_file_name <- file.path(data_dir, paste0(source_file_name, ".synthetic_copula.csv"))
write.csv(synthetic_df, output_csv_file_name, row.names = FALSE)

###
## Generate comparison plots
###
cat("\nStatus: Generating comparison charts...")

for (j in seq_along(header)) {
  idx <- header[j]
  x_real <- source_file_data[[idx]]
  x_syn  <- synthetic_df[[idx]]
  
  output_png_file_name <- file.path(
    plots_sub_dir,
    paste0("column_", j, ".", idx, "_original_v_synthetic.png")
  )
  
  if (is.factor(x_real)) {
    plot_factor_compare(
      original = x_real,
      synthetic = x_syn,
      var_name = capitalize(idx),
      out_file = output_png_file_name
    )
  } else {
    plot_numeric_compare(
      original = x_real,
      synthetic = x_syn,
      var_name = capitalize(idx),
      out_file = output_png_file_name
    )
  }
}

cat("\nStatus: Done generating synthetic data with Metalog + Gaussian Copula!")
cat("\nStatus: Synthetic file saved to: ", output_csv_file_name, "\n", sep = "")

