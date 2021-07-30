suppressMessages({
  load_lib <- c(
    # data manipulation
    "dplyr",
    "tidyverse",
    "tidymodels",
    # random forest
    "randomForest",
    "randomForestExplainer",
    "ROCR",
    # graphics
    "ggplot2",
    "ggraph",
    "igraph",
    "pheatmap",
    # data description
    "DescTools",
    "MVN",
    "psych"
  )
  install_lib <- load_lib[!(load_lib %in% installed.packages())] # check package
  if (length(install_lib)) for (i in install_lib) install.packages(i) # install

  cat("Loaded Packages:\n")
  print(sapply(load_lib, require, character = TRUE)) # load
  cat("\n\n")
})