################################################################################
# Data analysis of COVID-19 published at: (article submitted for publication)
# date of creation: 07/03/2021 (date in US format)
# R version: 4.0.5
# script name: script.R
# aim: data analysis
# input: files from the folder 'data'
# output: files saved in the subdirectories of folder 'outputs'
# external sources: packages.R, outputs.R located in the folder 'codes'
################################################################################

# USER PARAMETERS --------------------------------------------------------------
# figures
width_factor <- 4 # factor for determining figure width

height_factor <- 2 # factor for determining figure height

log_fc_cutpoint <- 1 # log fold change cutpoint

# SEED -------------------------------------------------------------------------
set.seed(2021) # fix seed

# EXTERNAL SOURCES -------------------------------------------------------------
source("./codes/packages.R") # load packages

# DATA -------------------------------------------------------------------------
# load dataset
wide_data <- read.csv("./data/my_data2.csv", header = TRUE, check.names = FALSE)

original_colnames <- colnames(wide_data) # store the original column names

k <- which(original_colnames == "Group") # column number of treatment variable

long_data <- reshape2::melt(wide_data, id.vars = "Group") # dataset in long form

# DESCRIPTIVE ANALYSIS ---------------------------------------------------------
summ_data <- summary(long_data$value) # basic statistics for all values

tab_desc <- MVN:::descriptives(wide_data[, -k]) # descriptive statistics

# INFERENTIAL ANALYSIS ---------------------------------------------------------
# non-parametric manova
colnames(wide_data)[-k] <- paste0("C", 1:(k - 1))
form <- paste0(paste(colnames(wide_data)[-k], collapse = "|"), " ~ Group")
suppressWarnings(
  npar_manova <- nonpartest(formula(form),
    data = wide_data, plot = FALSE, permtest = FALSE
  )
)

# non-parametric two sample tests
npar_t_tests <- lapply(colnames(wide_data[, -k]), function(.x) {
  test_df <- wide_data[, c(.x, "Group")]
  colnames(test_df)[1] <- "var"
  npar.t.test(var ~ Group, data = test_df, method = "t.app", info = FALSE)
})
npar_t_tests_res <- do.call(rbind, lapply(npar_t_tests, "[[", 2))
rownames(npar_t_tests_res) <- original_colnames[-k]

tab_npar_t_tests <- subset(npar_t_tests_res, p.Value < 0.05) # subset by p-value
tab_npar_t_tests$vars <- rownames(tab_npar_t_tests)
tab_npar_t_tests <- arrange(tab_npar_t_tests, Estimator)
tab_npar_t_tests$vars <- factor(tab_npar_t_tests$vars,
  levels = tab_npar_t_tests$vars, ordered = TRUE
)

# plot variables with significantly different relative effect on the treatment
fig_pval_relef <- ggplot(tab_npar_t_tests, aes(x = vars, y = Estimator)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", alpha = 0.7) +
  geom_point() +
  geom_ribbon(
    alpha = 0.4, fill = "grey",
    aes(x = seq_len(nrow(tab_npar_t_tests)), ymin = Lower, ymax = Upper)
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  theme_light(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, size = 6),
    panel.grid = element_blank()
  ) +
  labs(y = "Relative effect", x = "")

# log2 fold change (using untransformed data)
untransformed_data <- 2^wide_data[, -k]
untransformed_data$Group <- wide_data[, k]
group_means <- untransformed_data %>%
  group_by(Group) %>%
  summarise_all(mean)
log_fc <- apply(log2(group_means[, -1]), 2, diff)
names(log_fc) <- original_colnames[-k]

# volcano plot for log fold change x t-test p-value
vars <- match(rownames(npar_t_tests_res), names(log_fc))
volcano_data <- data.frame(fc = log_fc, p = npar_t_tests_res[vars, "p.Value"])
cond <- abs(volcano_data$fc) > log_fc_cutpoint & volcano_data$p < 0.05
volcano_data$sign <- ifelse(cond,
  ifelse(volcano_data$fc < 0, "neg", "pos"), "none"
)
volcano_data$labs <- names(log_fc)
volcano_data$labs[!cond] <- NA

fig_volcano_log_fc <- ggplot(volcano_data, aes(
  x = fc, y = p, label = labs, colour = sign
)) +
  geom_hline(yintercept = 0.05, linetype = "dashed", alpha = 0.3) +
  geom_vline(
    xintercept = c(-log_fc_cutpoint, log_fc_cutpoint),
    linetype = "dashed", alpha = 0.3
  ) +
  geom_point(size = 2, alpha = 0.4) +
  scale_colour_manual("Legend: ",
    values = c("neg" = "red", "pos" = "blue", "none" = "grey"),
    labels = c(
      "neg" = "ICU downregulated", "pos" = "ICU upreguladed",
      "none" = "Nonsignificant ICU/Non-ICU differences"
    )
  ) +
  geom_text_repel(
    show.legend = FALSE, max.overlaps = 10, direction = "both",
    angle = 0, segment.color = grey(0.9)
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  theme_light(base_size = 12) +
  theme(
    legend.position = "top",
    legend.box.background = ggplot2::element_rect(colour = "black", fill = NA)
  ) +
  labs(x = "Log fold change", y = "P-value for the t-test")

# OUTPUTS ----------------------------------------------------------------------
# save tables and figures
suppressMessages(suppressWarnings(source("./codes/outputs.R")))