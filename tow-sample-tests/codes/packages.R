suppressMessages({
  load_lib <- c(
    # data manipulation
    "reshape2",
    "dplyr",
    # graphics
    "ggplot2",
    "ggrepel",
    "scales",
    # data description
    "MVN",
    # statistical inference
    "npmv",
    "nparcomp"
  )
  install_lib <- load_lib[!(load_lib %in% installed.packages())] # check package
  if (length(install_lib)) for (i in install_lib) install.packages(i) # install

  cat("Loaded Packages:\n")
  print(sapply(load_lib, require, character = TRUE)) # load
  cat("\n\n")
})