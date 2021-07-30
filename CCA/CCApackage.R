loadPackagesCCA <- function(force_install_all_packages = FALSE) {
  cat("Installing and loading required packages ...", "\n")

  suppressMessages({
    load_lib <- c(
    #file management
    "filesstrings",
    #data manipulation:
      "readxl",
      "purrr",
      "DT",
      "rlist",
      "dplyr",
      "reshape2",
    #data description:
      "MVN",
    #cca:
      "CCA",
      "whitening",
    #clustering:
      "textshape",
    #graphics:
      "pheatmap",
      "ggplot2",
      "ggrepel",
      "GGally",
      "sna",
      "network",
      "ggnetwork"
  )

    if (force_install_all_packages) {
      install_lib <- load_lib
    } else {
      install_lib <- load_lib[!(load_lib %in% installed.packages())]
    }
    
    if(length(install_lib))
      for (lib in install_lib)
        install.packages(lib, dependencies = TRUE)

    loading <- sapply(load_lib, require, character = TRUE)
  })

  if (!any(loading)) {
    message("Packages not working: ",  load_lib[!loading])
    stop("The above packages were not properly installed.")
  } else {
    message("Packages successfully loaded.")
  }

  return(load_lib)
}

packagesReferences <- function() {
  ref <- paste("References: ",
    paste0("R package: whitening; ", "DOI: 10.1186/s12859-018-2572-9; ", "ISBN: ---"),
    paste0("R package: ggplot2; ", "DOI: ---; ", "ISBN: 331924275X"),
    sep = "\n\n")

  return(ref)
}

dataLoad <- function() {
  cat("A FILE CHOOSER WINDOW WAS OPENED:", "\n")
  cat("CHOOSE AN Excel FILE CONTAINING SPREADSHEETS WITH DATA ...", "\n")

  file <- choose.files(multi = FALSE)
  file_format <- excel_format(file)

  if (!is.na(file_format)) {
    dfs <- file %>% excel_sheets() %>% set_names() %>% map(read_excel, path = file)
    message("File successfully chosen.")
    cat("File format: ", file_format, "\nNumber of spreadsheets: ", length(dfs), "\n")
  } else {
    stop("File is not an Excel file.")
  }

  output <- list(dfs = dfs, n_sheets = length(dfs))

  return(output)
}

dataLogTransform <- function(dfs, base = 2) {
  if (base <= 0) stop("Logarithm base should be positive")

  log_dfs <- lapply(dfs, function(df) {
    if (any(df < 0)) {
      stop("Data with negative entries are not allowed.")
    } else if (any(df == 0)) {
      message("Zero values identified in the dataset.")
      df <- log(df + 1, base)
      message("Data was transformed using log", base,"(value + 1)")
    } else {
      df <- log(df, base)
      message("Data was transformed using log", base,"()")
    }
    return(df)
  })
  names(log_dfs) <- names(dfs)

  return(log_dfs)
}

dataSqrtTransform <- function(dfs) {
  sqrt_dfs <- lapply(dfs, function(df) {
    if (any(df < 0)) {
      stop("Data with negative entries are not allowed.")
    } else {
      df <- sqrt(df)
      message("Data was transformed using sqrt()")
    }
    return(df)
  })
  names(sqrt_dfs) <- names(dfs)
  
  return(sqrt_dfs)
}

dataScale <- function(dfs) {
  scale_dfs <- lapply(dfs, scale)
  names(scale_dfs) <- names(dfs)
  message("Data was scaled using scale().")

  return(scale_dfs)
}

statCCA <- function(X, Y, regularization = FALSE, n_points = 5) {
  t1 <- Sys.time()
  if (!regularization) {
    regul <- NULL
    cca <- cc(X, Y)
    cca_cor <- whitening::cca(as.matrix(X), as.matrix(Y))$lambda
  } else {
    cat("(this proccess can take a while) Searching for the optimum regularization parameter ...", "\n")
    grid <- list(lambda1 = seq(0.001, 1, length = n_points), lambda2 = seq(0.001, 1, length = n_points))
    regul <- estim.regul(X, Y, grid$lambda1, grid$lambda2, plt = FALSE)
    lambda <- c(regul$lambda1, regul$lambda2)
    cca <- rcc(X, Y, lambda[1], lambda[2])
    cca_cor <- cca$cor
    
    if (lambda[[1]] == max(grid$lambda1) | lambda[[2]] == max(grid$lambda2)) {
      message("Regularization parameter on the bounderies:", "(", lambda[[1]], ", ", lambda[[2]], ")")
      cat("Hint: try to increase 'rCCA_n_points'.", "\n")
    }
  }
  t2 <- Sys.time() - t1
  cat("Time elapsed running ", ifelse(regularization, "regularized CCA", "CCA"), " = ", t2, "\n")
  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  x <- ncol(X)
  y <- ncol(Y)

  S <- cov(cbind(X, Y))
  Sx <- S[1:x, 1:x]
  Sy <- S[(x + 1):(x + y), (x + 1):(x + y)]
  Sxy <- S[1:x, (x + 1):(x + y)]

  Rux <- cca$scores$corr.X.xscores
  Rvy <- cca$scores$corr.X.yscores

  Ru <- sum(diag(Sx))
  Rv <- sum(diag(Sy))

  if (Ru == x & Rv == y & !any(abs(signif(S, 5)) > 1)) {
    Azi2 <- lapply(as.data.frame(Rux), function(i) {
      m <- matrix(i, ncol = 1)
      m2 <- m %*% t(m)
    })
    Bzi2 <- lapply(as.data.frame(Rvy), function(i) {
      m <- matrix(i, ncol = 1)
      m2 <- m %*% t(m)
    })

    var_expl <- do.call(rbind, lapply(1:min(x, y), function(i) {
      Ru2 <- sum(diag(Reduce("+", Azi2[1:i]))) / Ru
      Rv2 <- sum(diag(Reduce("+", Bzi2[1:i]))) / Rv
      output <- c(U = Ru2, V = Rv2)
    }))
    colnames(var_expl) <- paste0("Cumulative_Variance_of_", c("U", "V"))

    message("Do not analyse the variance explained by CCA if your data was not scaled, cosider only the canonical correlations.")
  } else {
    value <- rep(0, min(x, y))
    var_expl <- data.frame(U = value, V = value)
    colnames(var_expl) <- paste0("Cumulative_Variance_of_", c("U", "V"))
  }

  output <- list(
    summary = c(n_var_X = x, n_var_Y = y),
    results = list(cca_cor = cca_cor, var_expl = var_expl),
    coefs = list(A = cca$xcoef, B = cca$ycoef),
    scores = list(U = cca$scores$xscores, V = cca$scores$yscores),
    correlations = list(Rux = cca$scores$corr.X.xscores, Rvy = cca$scores$corr.Y.yscores, 
                        Rvx = cca$scores$corr.Y.xscores, Ruy = cca$scores$corr.X.yscores),
    cca_data = list(X = X, Y = Y),
    pearson_cor = list(corX = Sx, corY = Sy, corXY = Sxy),
    regul = regul
  )

  if (!regularization) {
    message("CCA performed successfully.")
    cat("CCA performed successfully.", "\n")
  } else {
    message("Regularized CCA performed successfully.")
    cat("Regularized CCA performed successfully.", "\n")
  }

  return(output)
}

circleFun <- function(center = c(0, 0), diameter = 1, npoints = 100) {
  r <- diameter
  tt <- seq(0, 2 * pi, length.out = npoints)
  xx <- center[1] + r * cos(tt)
  yy <- center[2] + r * sin(tt)
  return(data.frame(x = xx, y = yy))
}

newDir <- function(dir_name, at_root = TRUE, root) {
  if (at_root) root <- getwd()

  path <- paste(root, dir_name, sep = "/")

  if (!dir.exists(path)) {
    result <- dir.create(path)
    if (result) message("A folder named '", dir_name, "' was created in: ", path)
  }

  return(path)
}

saveTable <- function(publish = FALSE, p, title, ...) {
  path <- newDir(...) 

  if (dir.exists(path)) {
    if(!publish) {
      path_csv <- newDir("csv_format", FALSE, path)
    } else {
      path_csv <- path
    }
    write.csv(p, file = paste(path_csv, paste0(title, ".csv"), sep = "/"), quote = FALSE, row.names = FALSE)
    
    if(!publish) {
      path_html <- newDir("html_format", FALSE, path)
      DT::saveWidget(datatable(round(p, 5)), file = paste(path_html, paste0(title, ".html"), sep = "/"), selfcontained = FALSE)
    } else {
      path_html <- "no html returned"
    }
    
    message("table '",  title, "' saved in: ", path_csv, ", and ", path_html, ".")
  } else {
    stop("It was not possible to create a folder to store tables.")
  }
}

savePlot <- function(p = NULL, title, width, height, extension = "pdf", ...) {
  path <- newDir(...) 

  if (dir.exists(path) & !is.null(p)) {
    suppressWarnings({
      if(extension == "pdf") {
        pdf_name <- paste0(title, ".pdf")
        pdf(paste(path, pdf_name, sep = "/"), width = width, height = height)
        print(p)
        dev.off()
      } else {
        png_name <- paste0(title, ".png")
        png(paste(path, png_name, sep = "/"), width = width, height = height, units = "in", res = 300)
        print(p)
        dev.off()
      }
    }) 

    message("plot '", title , "' saved in: ", path, ".")
  } else if (is.null(p)) {
    return(paste(path, paste(title, extension, sep = "."), sep = "/"))
  } else {
    stop("It was not possible to create a folder to store figures.")
  }
}

plotScores <- function(cca, canonical_variate = 1:2, id = 1, g = 1, g_col = "auto", cv_col = "gold", x_title = "X", 
                       y_title = "Y", short_x_title = "x", short_y_title = "y") {
  max_cvs <- min(ncol(cca$cca_data$X), ncol(cca$cca_data$Y))
  
  if (length(canonical_variate) > 2) stop("plots_canonical_variate should have exactly two entries.")
  if (any(canonical_variate < 0) | any(canonical_variate > (max_cvs - 1))) stop("Entries of plots_canonical_variate should be at least equal to 1 and most equal to ", max_cvs - 1, ".")
  
  if (g_col == "auto") {
    groups <- as.character(unique(g))
    g_col <- 1:length(groups)
    names(g_col) <- groups
  } else if (!((g_col == "auto") | is.vector(g_col))) {
    stop("'plots_group_color' should be a vector or 'auto'.")
  }
  
  groups <- as.character(unique(g))
  comparison <- setdiff(groups, names(g_col))
  if (length(comparison)) stop("The folowing groups do not have corresponding color in 'plots_group_color': ", paste(comparison, collapse = ", "))
  
  cca_x <- cca$scores$U[, canonical_variate]
  cca_y <- cca$scores$V[, canonical_variate]
  df_cca_x <- data.frame(cca_x, set = "X", id = id, group = g)
  df_cca_y <- data.frame(cca_y, set = "Y", id = id, group = g)
  colnames(df_cca_y)[1:2] <- colnames(df_cca_x)[1:2] <- paste0("CV", 1:2)
  
  df_cca <- rbind(df_cca_x, df_cca_y)
  df_cca$group <- as.character(df_cca$group)
  df_cca$id <- as.character(df_cca$id)
  
  facet_label <- c(X = x_title, Y = y_title)
  
  df_cca$label_x <- df_cca$label_y <- NA
  df_cca$label_x[grep("X", df_cca$set)[1]] <- paste0(short_x_title, "-CV", canonical_variate[1])
  df_cca$label_y[grep("X", df_cca$set)[1]] <- paste0(short_x_title, "-CV", canonical_variate[2])
  df_cca$label_x[grep("Y", df_cca$set)[1]] <- paste0(short_y_title, "-CV", canonical_variate[1])
  df_cca$label_y[grep("Y", df_cca$set)[1]] <- paste0(short_y_title, "-CV", canonical_variate[2])
  
  minX_x <- min(cca_x[, 1])
  minX_y <- min(cca_x[, 2])
  minY_x <- min(cca_y[, 1])
  minY_y <- min(cca_y[, 2])
  
  min_x <- max(minX_x, minY_x)
  min_y <- max(minX_y, minY_y)
  
  p <- ggplot(data = df_cca, aes(x = CV1, y = CV2, group = set, colour = group)) +
    geom_point(alpha = 0.6, size = 2, shape = 3) +
    geom_hline(yintercept = 0, alpha = 0.6, color = cv_col, size = 1.5) +
    geom_vline(xintercept = 0, alpha = 0.6, color = cv_col, size = 1.5) +
    geom_text(aes(label = label_x), x = min_x, y = 0.3, color = cv_col, size = 4.2, fontface = "bold") +
    geom_text(aes(label = label_y), x = 0.2, y = min_y, color = cv_col, angle = 90, size = 4.2, fontface = "bold") +
    facet_wrap(~ set, labeller = as_labeller(facet_label)) +
    scale_colour_manual("Group: ", values = g_col) +
    scale_y_continuous(breaks = seq(-100, 100, 1)) +
    scale_x_continuous(breaks = seq(-100, 100, 1)) +
    theme_classic() +
    theme(
      legend.position = "top",
      legend.box.background = element_rect(colour = "black"),
      panel.border = element_rect(colour = "black", fill = NA),
      panel.grid.major = element_line(),
      legend.text = element_text(size = 11),
      axis.line = element_blank(),
      axis.title = element_text(size = 15),
      strip.text = element_text(size = 13)) +
    labs(y = paste("Canonical variate", canonical_variate[2]),
         x = paste("Canonical variate", canonical_variate[1]))
  
  output <- list(p = p, df = df_cca)
  
  return(output)
}

plotCoef <- function(cca, canonical_variate = 1:2, x_title = "X", y_title = "Y", pos_col = "green", neg_col = "red") {
  max_cvs <- min(ncol(cca$cca_data$X), ncol(cca$cca_data$Y))
  
  if (length(canonical_variate) > 2) stop("plots_canonical_variate should have exactly two entries.")
  if (any(canonical_variate < 0) | any(canonical_variate > (max_cvs - 1))) stop("Entries of plots_canonical_variate should be at least equal to 1 and most equal to ", max_cvs - 1, ".")
  
  CoefRaw_X <- data.frame(cca$coefs$A[, canonical_variate], "X")
  CoefRaw_Y <- data.frame(cca$coefs$B[, canonical_variate], "Y")
  colnames(CoefRaw_Y) <- colnames(CoefRaw_X) <- c(paste0("CV", canonical_variate), "Set")
  
  CoefRaw_XY <- rbind(CoefRaw_X, CoefRaw_Y)
  CoefRaw_XY$Vars <- rownames(CoefRaw_XY)
  
  melt_CoefRaw_XY <- melt(CoefRaw_XY, id.vars = c("Set", "Vars"))
  melt_CoefRaw_XY$lables <- paste0(melt_CoefRaw_XY$Set, " (", melt_CoefRaw_XY$variable, ")")
  melt_CoefRaw_XY <- arrange(melt_CoefRaw_XY, variable, Set, -abs(value))
  melt_CoefRaw_XY$lables <- factor(melt_CoefRaw_XY$lables, levels = unique(melt_CoefRaw_XY$lables), ordered = TRUE)
  melt_CoefRaw_XY$Vars <- factor(melt_CoefRaw_XY$Vars, levels = unique(melt_CoefRaw_XY$Vars), ordered = TRUE)
  
  p <- ggplot(melt_CoefRaw_XY, aes(x = Vars, xend = Vars, yend = abs(value), group = Set, colour = ifelse(value < 0, "Negative", "Positive"))) +
    geom_segment(aes(y = 0)) +
    geom_point(aes(y = abs(value)), show.legend = FALSE) +
    facet_wrap(~ lables, scales = "free_x", ncol = 1) +
    scale_color_manual("Coefficient sign: ", values = c("Negative" = neg_col, "Positive" = pos_col)) +
    theme_classic() +
    theme(
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 13),
      axis.line = element_blank(),
      axis.title = element_text(size = 15),
      strip.text = element_text(size = 13),
      axis.text.x = element_text(angle = 90, size = 7),
      panel.grid.major.y = element_line(),
      panel.grid.minor.y = element_line(),
      legend.position = "top",
      legend.box.background = element_rect(colour = "black")) +
    labs(y = "Absolute value of the raw estimated coefficient for canonical variate", x = "")
  
  output <- list(p = p, df = melt_CoefRaw_XY)
  
  return(output)
}

plotStruct <- function(cca, canonical_variate = 1:2, x_title = "X", y_title = "Y", pos_col = "blue", neg_col = "red", k_min = 0.7) {
  max_cvs <- min(ncol(cca$cca_data$X), ncol(cca$cca_data$Y))
  
  if (length(canonical_variate) > 2) stop("plots_canonical_variate should have exactly two entries.")
  if (any(canonical_variate < 0) | any(canonical_variate > (max_cvs - 1))) stop("Entries of plots_canonical_variate should be at least equal to 1 and most equal to ", max_cvs - 1, ".")
  
  CoefStruct_X <- data.frame(cca$correlations$Rux[, canonical_variate], "X")
  CoefStruct_Y <- data.frame(cca$correlations$Rvy[, canonical_variate], "Y")

  colnames(CoefStruct_Y) <- colnames(CoefStruct_X) <- c(paste0("CV", canonical_variate), "Set")
  
  CoefStruct_XY <- rbind(CoefStruct_X, CoefStruct_Y)
  CoefStruct_XY$Vars <- rownames(CoefStruct_XY)
  
  melt_CoefStruct_XY <- melt(CoefStruct_XY, id.vars = c("Set", "Vars"))
  melt_CoefStruct_XY$lables <- paste0(melt_CoefStruct_XY$Set, " (", melt_CoefStruct_XY$variable, ")")
  melt_CoefStruct_XY <- arrange(melt_CoefStruct_XY, variable, Set, -abs(value))
  melt_CoefStruct_XY$lables <- factor(melt_CoefStruct_XY$lables, levels = unique(melt_CoefStruct_XY$lables), ordered = TRUE)
  melt_CoefStruct_XY$Vars <- factor(melt_CoefStruct_XY$Vars, levels = unique(melt_CoefStruct_XY$Vars), ordered = TRUE)
  
  p <- ggplot(melt_CoefStruct_XY, aes(x = Vars, xend = Vars, yend = abs(value), colour = ifelse(value < 0, "Negative", "Positive"))) +
    geom_segment(aes(y = 0)) +
    geom_point(aes(y = abs(value)), show.legend = FALSE) +
    geom_hline(yintercept = k_min, colour = "black", alpha = 0.7, linetype = "dashed", size = 1) +
    facet_wrap(~ lables, scales = "free_x", ncol = 1) +
    scale_color_manual("Correlation sign: ", values = c("Negative" = neg_col, "Positive" = pos_col)) +
    theme_classic() +
    theme(
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 13),
      axis.line = element_blank(),
      axis.title = element_text(size = 15),
      strip.text = element_text(size = 13),
      axis.text.x = element_text(angle = 90, size = 7),
      panel.grid.major.y = element_line(),
      panel.grid.minor.y = element_line(),
      legend.position = "top",
      legend.box.background = element_rect(colour = "black")) +
    labs(y = "Absolute value of the estimated correlation with canonical variate", x = "")
  
  output <- list(p = p, df = melt_CoefStruct_XY)
  
  return(output)
}

plotCorrel <- function(cca, pos_col = "blue", neg_col = "red", regularization = FALSE) {
  correl <- cca$results$cca_cor
  n <- length(correl)
  df_correl <- data.frame(cca_correl = correl, cca_pair = paste0(1:n, ""))
  df_correl <- arrange(df_correl, -abs(correl))
  df_correl$cca_pair <- factor(df_correl$cca_pair, levels = unique(df_correl$cca_pair), ordered = TRUE)
  
  if (regularization) {
    p <- ggplot(df_correl, aes(x = cca_pair, xend = cca_pair, yend = abs(cca_correl)))
  } else {
    p <- ggplot(df_correl, aes(x = cca_pair, xend = cca_pair, yend = abs(cca_correl), colour = ifelse(cca_correl < 0, "Negative", "Positive"))) +
      scale_color_manual("Correlation sign: ", values = c("Negative" = neg_col, "Positive" = pos_col))
  }
  
  p <- p + geom_segment(aes(y = 0)) + geom_point(aes(y = abs(cca_correl)), show.legend = FALSE) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    theme_classic() +
    theme(
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 13),
      axis.line = element_blank(),
      axis.title = element_text(size = 15),
      strip.text = element_text(size = 13),
      axis.text.x = element_text(angle = 90, size = 7),
      panel.grid.major.y = element_line(),
      panel.grid.minor.y = element_line(),
      legend.position = "top",
      legend.box.background = element_rect(colour = "black")) +
    labs(y = "Absolute value of the estimated canonical correlation", x = "Canonical pair")
  
  output <- list(p = p, df = df_correl)
  
  return(output)
}

plotCumVarExpl <- function(cca, x_title = "X", y_title = "Y", x_col = "red", y_col = "blue") {
  var_expl <- cca$results$var_expl
  n <- nrow(var_expl)
  df_var_expl <- data.frame(var_expl, cca_pair = 1:n)
  df_var_expl$cca_pair <- factor(df_var_expl$cca_pair, levels = unique(df_var_expl$cca_pair), ordered = TRUE)
  
  melt_df_var_expl <- melt(df_var_expl, id.vars = "cca_pair")
  
  p <- ggplot(melt_df_var_expl, aes(x = cca_pair, y = value, group = variable, colour = variable)) +
    geom_line() +
    geom_point() +
    scale_y_continuous(breaks = seq(0, 1, 0.1), labels = scales::percent_format(accuracy = 1)) +
    scale_x_discrete(labels = paste0(melt_df_var_expl$cca_pair, "")) +
    scale_colour_manual("Dataset: ", 
                        values = c("Cumulative_Variance_of_U" = x_col, "Cumulative_Variance_of_V" = y_col),
                        labels = c("Cumulative_Variance_of_U" = x_title, "Cumulative_Variance_of_V" = y_title)) +
    theme_classic() +
    theme(
      legend.text = element_text(size = 13),
      legend.title = element_text(size = 13),
      axis.line = element_blank(),
      axis.title = element_text(size = 15),
      strip.text = element_text(size = 13),
      axis.text.x = element_text(angle = 90, size = 7),
      panel.grid.major = element_line(),
      panel.grid.minor = element_line(),
      legend.position = "top",
      legend.box.background = element_rect(colour = "black")) +
    labs(y = "Cumulative variance explained", x = "Canonical pair")
  
  output <- list(p = p, df = melt_df_var_expl)
  
  return(output)
}

plotCCA <- function(cca, canonical_variate = 1:2, k_min = 0.6, x_title = "X", y_title = "Y", short_x_title = "X",
  short_y_title = "Y", x_col = "red", y_col = "blue", cv_col = "gold", max_overlaps = 10^6) {
  max_cvs <- min(ncol(cca$cca_data$X), ncol(cca$cca_data$Y))

  if (length(canonical_variate) > 2) stop("plots_canonical_variate should have exactly two entries.")
  if (any(canonical_variate < 0) | any(canonical_variate > (max_cvs - 1))) stop("Entries of plots_canonical_variate should be at least equal to 1 and most equal to ", max_cvs - 1, ".")
  if (k_min > 1) stop("plotCCA_corr_min should be at between 0 and 1, inclusive.")

  cca_x <- cca$correlations$Rux[, canonical_variate]
  df_cca_x <- data.frame(cca_x, Vars = colnames(cca$cca_data$X), group = "X",
   col = sapply(seq(1, 2, 2), function(i) {
    df <- cca_x[, i:(i + 1)]
    apply(df, 1, function(j) factor((sum(abs(j) >= k_min) > 0) * 1))
   }), check.names = FALSE)
  colnames(df_cca_x)[1:2] <- paste0("CV", 1:2)

  cca_y <- cca$correlations$Rvy[, canonical_variate]
  df_cca_y <- data.frame(cca_y, Vars = colnames(cca$cca_data$Y), group = "Y",
   col = sapply(seq(1, 2, 2), function(i) {
    df <- cca_y[, i:(i + 1)]
    apply(df, 1, function(j) factor((sum(abs(j) >= k_min) > 0) * 2))
   }), check.names = FALSE)
  colnames(df_cca_y)[1:2] <- paste0("CV", 1:2)

  df_cca <- rbind(df_cca_x, df_cca_y)

  facet_label <- c(X = x_title, Y = y_title)

  df_cca$label_x <- df_cca$label_y <- NA
  df_cca$label_x[grep("X", df_cca$group)[1]] <- paste0(short_x_title, "-CV", canonical_variate[1])
  df_cca$label_y[grep("X", df_cca$group)[1]] <- paste0(short_x_title, "-CV", canonical_variate[2])
  df_cca$label_x[grep("Y", df_cca$group)[1]] <- paste0(short_y_title, "-CV", canonical_variate[1])
  df_cca$label_y[grep("Y", df_cca$group)[1]] <- paste0(short_y_title, "-CV", canonical_variate[2])

  p <- ggplot(data = df_cca, aes(x = CV1, y = CV2, colour = col, group = group)) +
   geom_point(alpha = 1, size = 2, show.legend = FALSE) +
   geom_rect(xmin = -1, xmax = 1, ymin = -1, ymax = 1, colour = "black", fill = "transparent",
    linetype = "dotted") +
   geom_rect(xmin = -k_min, xmax = k_min, ymin = -k_min, ymax = k_min, colour = "black",
    fill = "transparent", linetype = "dotted") +
   geom_hline(yintercept = 0, alpha = 0.6, color = cv_col, size = 1.5) +
   geom_vline(xintercept = 0, alpha = 0.6, color = cv_col, size = 1.5) +
   geom_text(aes(label = label_x), x = 0.2, y = 0.05, color = cv_col, size = 4.2, fontface = "bold") +
   geom_text(aes(label = label_y), x = -0.05, y = -0.2, color = cv_col, angle = 90, size = 4.2, fontface = "bold") +
   facet_wrap(~ group, labeller = as_labeller(facet_label)) +
   coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
   scale_colour_manual(values = c("0" = "grey", "1" = x_col, "2" = y_col)) +
   geom_text_repel(aes(label = ifelse(col != "0", Vars, NA)), size = 3.3,
   segment.alpha = 0.4, segment.size = 0.3, min.segment.length = 0, fontface = "bold",
   max.overlaps = max_overlaps, show.legend = FALSE) +
   theme_classic() +
   theme(
    legend.text = element_text(size = 11),
    axis.line = element_blank(),
    axis.title = element_text(size = 15),
    strip.text = element_text(size = 13)) +
   labs(y = paste("Correlation with canonical variate", canonical_variate[2]),
    x = paste("Correlation with canonical variate", canonical_variate[1]))

  output <- list(p = p, df = df_cca)

  return(output)
}

plotCor <- function(corr_mat, title, clust_method = "ward.D", color_ramp, border_color, fontsize) {
  cor_clust <- cluster_matrix(corr_mat, dim = "both", method = clust_method)
  
  p <- pheatmap(cor_clust, cluster_rows = FALSE, cluster_cols = FALSE, 
    color = color_ramp, border_color = border_color, fontsize = fontsize)
  
  output <- list(p = p, df = cor_clust)
  
  return(output)
}

plotCorNetwork <- function(corr_mat, what = "X", cca = NULL, k = 0.9, seed = 100, node_size = 5,
  label_size = 2.3, x_col = "red", y_col = "blue", cv_col = "gold", neg_corr_col = "red", pos_corr_col = "blue",
  net_method = "kamadakawai") {
  diag(corr_mat) <- 0
  
  message("Just a check: correlation matrix provided to plotCorNetwork is symmetric? ", isSymmetric(corr_mat))
  
  if (what %in% c("all", "XY")) {
    if (is.null(cca)) stop("Parameter 'cca' is NULL in plotCorNetwork.")
    
    X <- cca$cca_data$X
    Y <- cca$cca_data$Y
    
    if (what == "XY") corr_mat <- cor(data.frame(X, Y))
    
    type_Vars <- ifelse(
      colnames(corr_mat) %in% colnames(X), "X",
      ifelse(colnames(corr_mat) %in% colnames(Y), "Y", "CV")
    )
  } else if (what %in% c("X", "Y")) {
    type_Vars <- rep(what, ncol(corr_mat))
  } else {
    stop("Parameter 'what' should be one of 'all', 'XY', 'X' or 'Y'.")
  }
  
  cor_Vars_cv_k <- corr_mat
  cor_Vars_cv_k[abs(cor_Vars_cv_k) < k] <- 0
  cor_Vars_cv_k[abs(cor_Vars_cv_k) >= k] <- 1
  
  name_Vars <- colnames(corr_mat)
  name_Vars[colSums(cor_Vars_cv_k) == 0] <- ""
  
  net <- network(cor_Vars_cv_k, directed = FALSE)
  net %v% "type" <- type_Vars
  net %v% "name" <- name_Vars
  
  set.seed(seed)
  df_network <- ggnetwork(net, layout = net_method)
  
  pairs <- match(
    paste0(df_network[, "xend"], "/", df_network[, "yend"]),
    paste0(df_network[, "x"], "/", df_network[, "y"]))
  df_network$connect_to <- df_network$vertex.names[pairs]
  
  correlations <- melt(corr_mat)
  pairs <- match(
    paste0(df_network$vertex.names, "/", df_network$connect_to),
    paste0(correlations$Var1, "/", correlations$Var2))
  df_network$corr <- correlations$value[pairs]
  
  df_network$sign <- ifelse(df_network$corr >= 0, "pos", "neg")
  
  p <- ggplot(data = df_network, aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_edges(alpha = 0.4, size = 0.5, aes(color = sign), show.legend = FALSE) +
    geom_nodes(aes(color = type), alpha = 0.5, size = node_size, show.legend = FALSE) +
    geom_nodetext(aes(label = name), fontface = "bold", size = label_size, alpha = 0.7) +
    scale_colour_manual("", values = c("X" = x_col, "Y" = y_col, "CV" = cv_col, "pos" = pos_corr_col, "neg" = neg_corr_col)) +
    scale_size_identity() +
    theme_blank()
  
  output <- list(p = p, df = df_network)
  
  return(output)
}

corCCA <- function(cca, n = 2, canonical_variate = NULL, short_x_title = 'X', short_y_title = 'Y') {
  max_cvs <- min(ncol(cca$cca_data$X), ncol(cca$cca_data$Y))
  
  if (!is.null(canonical_variate)) {
    if (n <= length(canonical_variate)) {
      if(any(canonical_variate < 0)) stop("There cannot be a negative number in plots_canonical_variate.")
      n <- length(canonical_variate)
    } else {
      canonical_variate <- NULL
    }
  }
  if (n < 0 | n > max_cvs) stop("the number of canonical variates should be at least equal to 1 and at most equal to ", max_cvs, ".")

  X <- cca$cca_data$X
  Y <- cca$cca_data$Y
  signed_cor <- cca$results$cca_cor

  cor_Vars_cca <- cor(cbind(X, Y))

  if (!is.null(canonical_variate)) {
    n_cor <- canonical_variate
  } else {
    n_cor <- seq_len(n)
  }
  
  cor_cv <- signed_cor[n_cor]

  u_cv <- do.call(cbind, lapply(n_cor, function(i) {
    cor_ux_cv <- cca$correlations$Rux[, i]
    cor_uy_cv <- cca$correlations$Rvx[, i]
    u_cv <- c(cor_ux_cv, cor_uy_cv)
  }))

  v_cv <- do.call(cbind, lapply(n_cor, function(i) {
    cor_vx_cv <- cca$correlations$Ruy[, i]
    cor_vy_cv <- cca$correlations$Rvy[, i]
    v_cv <- c(cor_vx_cv, cor_vy_cv)
  }))

  diag_cor <- diag(cor_cv)
  diag_1 <- diag(2 * n)
  diag_1[1:n, (n + 1):ncol(diag_1)] <- diag_cor
  diag_1[(n + 1):nrow(diag_1), 1:n] <- diag_cor

  cv_cols <- cbind(u_cv, v_cv)
  cv_rows <- cbind(t(cv_cols), diag_1)
  bind_cols <- cbind(cor_Vars_cca, cv_cols)
  bind_rows <- rbind(bind_cols, cv_rows)

  colnames(bind_rows) <- rownames(bind_rows) <- c(
   colnames(cor_Vars_cca), paste0(short_x_title, "-CV", n_cor), paste0(short_y_title, "-CV", n_cor)
  )

  return(bind_rows)
}

formatCCATables <- function(cca) {
  results <- cca$results
  cancor_tables <- c(cca$coefs, cca$scores, cca$correlations[c("Rux", "Rvy")])
  pearson_tables <- cca$pearson_cor

  if (!is.null(results$var_expl)) {
    df_results <- round(data.frame(results$var_expl, Canonical_Correlation = results$cca_cor), 5)
  } else {
    df_results <- round(data.frame(-1, -1, Canonical_Correlation = results$cca_cor), 5)
  }
  colnames(df_results) <- c("Cumulative Variance of CV (X)", "Cumulative Variance of CV (Y)", "Canonical Correlation of CV Pair")
  formated_results <- list(df_results)
  names(formated_results) <- "CumulativeVarExpl_CanonicalCorr"
  
  formated_cancor <- lapply(cancor_tables, function(.x) {
    colnames(.x) <- paste0("CVpair", 1:ncol(.x))
    .x <- round(.x, 5)
    return(.x)
  })
  names(formated_cancor) <- names(cancor_tables)
  
  formated_pearson <- lapply(pearson_tables, function(.x) {
    .x <- round(.x, 5)
    return(.x)
  })
  names(formated_pearson) <- names(pearson_tables)
  
  output <- c(formated_results, formated_cancor, formated_pearson)
  output <- output[c("corX", "corY", "corXY", "Rux", "Rvy", "A", "B", "U", "V", "CumulativeVarExpl_CanonicalCorr")]
  names(output) <- c("[p1 and p4]", "[p2 and p5]", "[p3 and p6]", "[p7 and p10 for X]", "[p7 and p10 for Y]", "[p9 for X]", "[p9 for Y]", "[p13 for X]", "[p13 for Y]", "[p11 and p12]")
  
  return(output)
}

tableLegends <- function(cca_tables, x_title, y_title) {
  intro1 <- paste0("This legends make references to figures and tables.",
    " Use the identifications between brackets [] to correspond figures and tables.")
  
  leg_CorX <- paste0("Matrix of Pearson correlations for variables in ", x_title, ".")
  
  leg_CorY <- paste0("Matrix of Pearson correlations for variables in ", y_title, ".")
  
  leg_CorXY <- paste0("Matrix of Pearson correlations for variables of ", x_title, " against variables of ", y_title, ".")

  leg_CoefStructX <- paste0("Estimated correlations of variables from ", x_title, " with their corresponding canonical variates.")
  
  leg_CoefStructY <- paste0("Estimated correlations of variables from ", y_title, " with their corresponding canonical variates.")
  
  leg_CoefEstX <- paste0("Raw estimated coefficients for the canonical variates for ", x_title, ".")
  
  leg_CoefEstY <- paste0("Raw estimated coefficients for the canonical variates for ", y_title, ".")
  
  leg_CumVarCC <- paste0("(leftmost column) Cumulative variance of ", x_title, " explained by its corresponding canonical variates.",
    "(middle column) Cumulative variance of ", y_title, " explained by its corresponding canonical variates.",
    "(rightmost column) Canonical correlation for the n-th pair of canonical variates.")
  
  leg_ScoresX <- paste0("Scores of the individuals on each canonical variate obtained for ", x_title, ".",
    " Individuals are sorted by ID.")
  
  leg_ScoresY <- paste0("Scores of the individuals on each canonical variate obtained for ", y_title, ".",
    " Individuals are sorted by ID.")
  
  output <- paste0(intro1, "\n\n", paste(paste(names(cca_tables), c(leg_CorX, leg_CorY, leg_CorXY, leg_CoefStructX, leg_CoefStructY, leg_CoefEstX, leg_CoefEstY, leg_ScoresX, leg_ScoresY, leg_CumVarCC), sep = ": "), sep = "\n\n"))
  
  return(output)
}

descMeasures <- function(df1, df2, print_out = TRUE) {
  mvn_X <- mvn(df1)
  mvn_Y <- mvn(df2)
  
  desc_X <- arrange(round(mvn_X$Descriptives, 5), Skew, Kurtosis)
  desc_Y <- arrange(round(mvn_Y$Descriptives, 5), Skew, Kurtosis)
  n_var <- floor(min(nrow(df1), nrow(df2)) / 20)
  
  mvn_X$Descriptives <- desc_X
  mvn_Y$Descriptives <- desc_Y
  
  if (print_out) {
    cat("Dataset X -------------------------------------------------------------------------", "\n\n")
    print(mvn_X)
    cat("Dataset Y -------------------------------------------------------------------------", "\n\n")
    print(mvn_Y)
    cat("Suggested Number of Variables -----------------------------------------------------", "\n\n")
    cat("It is suggested to use at most ", n_var, " variables in a non-regularized CCA.", "\n")
    cat("User has chosen to keep ", nrow(desc_X) + nrow(desc_Y), " variables.", "\n\n")
    cat("-----------------------------------------------------------------------------------", "\n")
  }
  
  output <- list(desc_X = mvn_X, desc_Y = mvn_Y)
  
  return(output)
}

figLegends <- function(x_title, y_title, x_col, y_col, cv_col, pos_col, neg_col, graph_corr_min, cca_corr_min, canonical_variate, n_canonical_variates, regularization) {
  if(length(canonical_variate) < n_canonical_variates) {
    text1 <- paste0(canonical_variate[1], " and ", canonical_variate[2])
  } else {
    text1 <- paste0("first ", n_canonical_variates)
  }
  
  intro1 <- paste0("This legends make references to figures and tables.",
    " Use the identifications between brackets [] to correspond figures and tables.")
  
  leg_p1 <- paste0(
    "[p1](heatmap_matrix-of-pearson-correlations-for-variables-from-set-X): ",
    "Heatmap of Pearson correlations between variables of ", x_title, ".")
  
  leg_p2 <- paste0(
    "[p2](heatmap_matrix-of-pearson-correlations-for-variables-from-set-Y): ",
    "Heatmap of Pearson correlations between variables of ", y_title, ".")
  
  leg_p3 <- paste0(
    "[p3](heatmap_matrix-of-pearson-correlations-for-variables-from-set-XY): ",
    "Heatmap of Pearson correlations of ", x_title, " variables (on the vertical position) against ", y_title,
    " variables (on the horizontal position).")
  
  leg_p4 <- paste0(
    "[p4](graph_matrix-of-pearson-correlations-for-variables-from-set-X): ",
    "Pearson correlations for variables of ", x_title,
    " represented as a graph. Correlated variables are represented by connected nodes.",
    " Nodes with correlation in absolute value below ", graph_corr_min, " are not connected and have their names ommited.",
    " Positive and negative correlations are represented by edges colored ", pos_col, " and ", neg_col,", respectively.")
  
  leg_p5 <- paste0(
    "[p5](graph_matrix-of-pearson-correlations-for-variables-from-set-Y): ",
    "Pearson correlations for variables of ", y_title,
    " represented as a graph. Correlated variables are represented by connected nodes.",
    " Nodes with correlation in absolute value below ", graph_corr_min, " are not connected and have their names ommited.",
    " Positive and negative correlations are represented by edges colored ", pos_col, " and ", neg_col,", respectively.")
  
  leg_p6 <- paste0(
    "[p6](graph_matrix-of-pearson-correlations-for-variables-from-set-XY): ",
    "Pearson correlations (represented as a graph) for variables of ", x_title,
    " against variables of ", y_title, ".",
    " Variables from ", x_title, " are colored ", x_col,
    " while variables from ", y_title, " are colored ", y_col,
    ". Correlated variables are represented by connected nodes.",
    " Nodes with correlation in absolute value below ", graph_corr_min, " are not connected and have their names ommited.",
    " Positive and negative correlations are represented by edges colored ", pos_col, " and ", neg_col,", respectively.")
  
  leg_p7 <- paste0(
    "[p7](2Dplane_correlations-for-variables-with-canonical-variates): ",
    "(on the left) Estimated correlations for variables of ", x_title, " against its corresponding ", text1, " canonical variates.",
    " (on the right) Estimated correlations for variables of ", y_title, " against its corresponding ", text1, " canonical variates.",
    " Variables colored grey (and with names ommited) are those with correaltion in absolute value below ", cca_corr_min, " with any of its two corresponding canonical variates.",
    " Inner dotted lines mark the limits for correlation between ", -cca_corr_min, " and ", cca_corr_min, ", while outer dotted lines mark the limits for correlations between -1 and 1.")

  leg_p8 <- paste0(
    "[p8](graph_matrix-of-pearson-correlations-for-variables-from-set-XandYandCanonicalVariates): ",
    "Pearson correlations (represented as a graph) for all variable of ", x_title, " and ", y_title,
    ", also including the ", text1, " pairs of canonical variates.",
    " Variables from ", x_title, " are colored ", x_col,
    " while variables from ", y_title, " are colored ", y_col,
    ". Canonical variates are colored ", cv_col, ".",
    " Nodes with correlation in absolute value below ", graph_corr_min, " are not connected and have their names ommited.",
    " Positive and negative correlations are represented by edges colored ", pos_col, " and ", neg_col,", respectively.")
  
  leg_p9 <- paste0(
    "[p9](rw-estimated-coefficients-of-canonical-variables): ",
    " Raw estimated coefficients of the ", text1, " pairs of canonical variates",
    " for ", x_title, " and ", y_title, ".",
    " Positive and negative coefficients are represented by lines colored ", pos_col, " and ", neg_col,", respectively.")
  
  leg_p10 <- paste0(
    "[p10](1Dlines_correlations-for-variables-with-canonical-variates): ",
    "Estimated correlations of variables from ", x_title, " and ", y_title,
    " with their corresponding canonical variates. Only the ", text1, " canonical variates are shown.",
    " Positive and negative correlations are represented by lines colored ", pos_col, " and ", neg_col,", respectively.",
    " Dashed black horizontal line marks the correlation level (in absolute value) of ", cca_corr_min, ".")
  
  leg_p11 <- paste0(
    "[p11](cumulative-variance-explained-by-canonical-variates): ",
    "Cumulative variance explained by the canonical variates. Cumulative variance explained for the sets ", x_title, " and ", x_title,
    " are shown by lines colored ", x_col, " and ", y_col, ", respectively.")
  
  leg_p12 <- paste0(
    "[p12](canonical-correlation-between-pairs-of-canonical-variates): ",
    "Estimated Canonical correlations for pairs of canonical variates.",
    ifelse(regularization,
      " The estimated canonical correlation are unsigned.",
      paste0(" Positive and negative correlations are represented by lines colored ", pos_col, " and ", neg_col,", respectively.")
    )
  )
  
  leg_p13 <- paste0(
    "[p13](scores-of-individuals-on-canonical-variates): ",
    "(left) Scores for individuals on canonical variates for ", x_title, ".",
    " (right) Scores for individuals on canonical variates for ", y_title, ".",
    " Groups of individuals are colored according to the plot legend."
  )
  
  output <- paste(intro1, leg_p1, leg_p2, leg_p3, leg_p4, leg_p5, leg_p6, leg_p7, leg_p8, leg_p9, leg_p10, leg_p11, leg_p12, leg_p13, sep = "\n\n")
  
  return(output)
}

checkParameters <- function(list_param) {
  parameters <- list(
    empty = FALSE,
    force_install_all_packages = FALSE,
    print_descriptive_statistics = TRUE,
    rCCA_regularization = FALSE,
    rCCA_n_points = 10,
    heatmap_clust_method = "ward.D2",
    heatmap_color_ramp = RColorBrewer::brewer.pal(n = 7, name = "RdYlBu"),
    heatmap_border_color = grey(0.4),
    heatmap_fontsize = 5,
    heatmap_width = 13,
    heatmap_height = 7,
    plotCorCCA_n_canonical_variates = 2,
    plotCorCCA_corr_min = 0.7,
    plotCorCCA_seed = 10,
    plotCorCCA_node_size = 5,
    plotCorCCA_label_size = 3,
    plotCorCCA_net_method = "kamadakawai",
    plotCorCCA_width = 13,
    plotCorCCA_height = 7,
    plotCCA_corr_min = 0.7,
    plotCCA_width = 13,
    plotCCA_height = 7,
    plotCCA_max_overlaps = 1000,
    plotEst_width = 20,
    plotEst_height = 15,
    plotCumVar_width = 20,
    plotCumVar_height = 10,
    plotCC_width = 20,
    plotCC_height = 10,
    plots_canonical_variate = c(1, 2),
    plots_group_color = "auto",
    plots_x_title = "Dataset X",
    plots_y_title = "Dataset Y",
    plots_short_x_title = "x",
    plots_short_y_title = "y",
    plots_x_color = "green",
    plots_y_color = "pink",
    plots_cv_color = "gold",
    plots_neg_color = "red",
    plots_pos_color = "blue",
    plots_file_extension = "pdf",
    data_log_trans = FALSE,
    data_log_base = 2,
    data_sqrt_trans = FALSE,
    data_scale_trans = TRUE,
    data_X = NA,
    data_Y = NA,
    data_sheet1 = 1,
    data_sheet2 = 2
  )

  if (!is.null(list_param)) {
    user_param <- names(list_param)
    cond <- any(is.na(match(user_param, names(parameters))))
    if (cond) {
      stop("there are some undefined parameters on 'list_param'.")
    } else {
      for (i in user_param) parameters[[i]] <- list_param[[i]]
    }
  }

  return(parameters)
}

routineCCA <- function(list_param = list(empty = TRUE), publish = FALSE) {
  old_path <- new_path <- newDir("output_CCA")
  
  message_log <- file(paste(new_path, "messages.log", sep = "/"), open = "a")
  sink(message_log, type = "message", append = FALSE)
  
  outputs_log <- file(paste(new_path, "outputs.log", sep = "/"), open = "a")
  sink(outputs_log, type = "output", split = TRUE, append = FALSE)
  
  message("current work directory: ", getwd())
  
  parameters <- checkParameters(list_param)
  for (i in names(parameters)) assign(i, parameters[[i]])

  loaded_packages <- loadPackagesCCA(force_install_all_packages)

  if(!publish) {
    current_output <- gsub(":| ", "-", format(Sys.time(), "%Y-%m-%d_%X"))
    new_path <- newDir(current_output, at_root = FALSE, root = new_path)

    list.save(parameters, file = paste(new_path, "parameters.json", sep = "/"))
    cond0 <- "parameters.json" %in% list.files(new_path, pattern = "^(parameters.json)$")
    if(!cond0) stop("Could not save the file 'parameters.json' in ", new_path, ".")
  }

  if (publish) new_path <- newDir("R_CCA", at_root = FALSE, root = new_path)

  if (!is.data.frame(data_X) | !is.data.frame(data_Y)){
    message("One of the datasets, X or Y, are not data.frame or is NA.")
    cat("Load a dataset ...", "\n")

    dfs <- dataLoad()
    n_sheets <- dfs$n_sheets

    if (n_sheets < 2) {
      stop("Excel file should have at least 2 spreadsheets.")
    } else if ((data_sheet1 > 0 & data_sheet1 <= n_sheets) & (data_sheet2 > 0 & data_sheet2 <= n_sheets)) {
      cat("Spreadsheets ", data_sheet1, " and ", data_sheet2, " were selected for ", ifelse(rCCA_regularization, "regularized CCA", "CCA"), ".", "\n")
    } else {
      data_sheet1 = 1
      data_sheet2 = 2
      message("Reset variables to: data_sheet1 = 1 and data_sheet2 = 2.")
      cat("Only the first two spreadsheets were used for ", ifelse(rCCA_regularization, "regularized CCA", "CCA"), ".", "\n")
    }

    dfs$dfs <- dfs$dfs[c(data_sheet1, data_sheet2)]
  } else {
    dfs <- list(dfs = list(data_X, data_Y))
    n_sheets <- NULL
  }

  colnames_check <- sum(sapply(dfs$dfs, function(.x) sum(colnames(.x)[1:2] != c("GROUP", "ID"))))
  if (colnames_check) stop("In both datasets, first and second columns should be named GROUP and ID, respectively.")
  
  ids <- lapply(dfs$dfs, function(.x) .x$ID)
  id_check_1 <- any(ids[[1]] != ids[[2]])
  id_check_2 <- length(setdiff(ids[[1]], ids[[2]]))
  id_check_3 <- any(table(ids[[1]]) > 1)
  if (id_check_1 & id_check_2 & id_check_3) stop("There are unmached or repeated ids in the datasets.")
  
  dfs_names <- names(dfs$dfs)
  dfs_sorted <- lapply(1:2, function(.x) dfs$dfs[[.x]][order(ids[[.x]]), ])
  names(dfs_sorted) <- dfs_names
  ID <- dfs_sorted[[1]]$ID
  GROUP <- dfs_sorted[[1]]$GROUP
  dfs$dfs <- lapply(dfs$dfs, function(.x) .x[, -(1:2)])
  
  dimXY <- sapply(dfs$dfs, dim)
  if (diff(dimXY[1, ])) stop("Datasets does not have the same number of rows.")
  if (any(dimXY[2, ] == 0)) stop("Each dataset must have at least one variable.")
  if (!rCCA_regularization & (sum(dimXY[2, ]) > dimXY[1, 1])) {
    message("Data dimensionality: number of observations = ", dimXY[1, 1], " and number of variables selected = ", sum(dimXY[2, ]))
    stop("Number of variables should be less than the number of observations. Try rCCA_regulatization = TRUE.")
  }

  if (data_log_trans & data_sqrt_trans) {
    message("Log transformation will be used over square root transformation.")
    dfs$dfs <- dataLogTransform(dfs$dfs, base = data_log_base)
  } else if (data_log_trans) {
    dfs$dfs <- dataLogTransform(dfs$dfs, base = data_log_base)
  } else if (data_sqrt_trans) {
    dfs$dfs <- dataSqrtTransform(dfs$dfs)
  }
  
  if (data_scale_trans) dfs$dfs <- dataScale(dfs$dfs)

  data_X <- dfs$dfs[[1]]
  data_Y <- dfs$dfs[[2]]

  cca <- statCCA(data_X, data_Y, rCCA_regularization, rCCA_n_points)

  cat("Generating tables ...", "\n")
  
  cca_tables <- formatCCATables(cca)
  saved_tables <- sapply(1:length(cca_tables), function(.x) {
    saveTable(FALSE, cca_tables[[.x]], names(cca_tables)[.x], dir_name = "tables", at_root = FALSE, root = new_path)
  })

  if (publish) {
    newX <- data.frame(GROUP = GROUP, ID = ID, cca$cca_data$X)
    newY <- data.frame(GROUP = GROUP, ID = ID, cca$cca_data$Y)
    saveTable(TRUE, newX, "dataset_X", dir_name = "data", at_root = FALSE, root = new_path)
    saveTable(TRUE, newY, "dataset_Y", dir_name = "data", at_root = FALSE, root = new_path)

    path_code <- newDir("code", at_root = FALSE, root = new_path)

    cond1 <- "CCApackage.R" %in% list.files(path_code, pattern = "^(CCApackage.R)$")
    if (!cond1) {
      cond2 <- length(list.files(getwd(), pattern = "^(CCApackage.R)$"))
      if (!cond2) stop("file named 'CCApackage.R' not found.")

      file_to <- paste(path_code, "CCApackage.R", sep = "/")
      file_from <- paste(getwd(), "CCApackage.R", sep = "/")

      file.create(file_to)
      filecopy <- file.copy(from = file_from, to = file_to, overwrite = TRUE)
      if (!filecopy) stop("Could not copy file 'CCApackage.R' from workspace to ", path_code, ".")
    }

    parameters[["data_log_trans"]] <- FALSE
    parameters[["data_scale_trans"]] <- FALSE
    parameters[["data_X"]] <- newX
    parameters[["data_Y"]] <- newY
    list.save(parameters, file = paste(path_code, "list_param.rdata", sep = "/"))
    cond3 <- length(list.files(path_code, pattern = "^(list_param.rdata)$"))
    if(!cond3) stop("Could not save file 'list_param.rdata' in ", path_code, ".")

    write_package <- "if (!require('rlist')) install.packages('rlist')\nlibrary(rlist)"
    write_path <- "path <- paste(getwd(), 'code', sep = '/')"
    write_source <- "source(paste(path, 'CCApackage.R', sep = '/'))"
    write_param <- "list_param <- list.load(paste(path, 'list_param.rdata', sep = '/'))"
    write_runCCA <- "cca <- runCCA(list_param, publish = FALSE)"
    write.table(c(write_package, write_path, write_source, write_param, write_runCCA), file = paste(path_code, "main.R", sep = "/"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    cond4 <- length(list.files(path_code, pattern = "^(main.R)$"))
    if(!cond4) stop("Could not save file 'main.R' in ", path_code, ".")

    path_ref <- newDir("references-and-packages", at_root = FALSE, root = new_path)
    write.table(packagesReferences(), file = paste(path_ref, "main_references.txt", sep = "/"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(data.frame(loaded_packages), file = paste(path_ref, "used_packages.csv", sep = "/"), quote = FALSE, row.names = FALSE, col.names = FALSE)

    readme <-  "
      # Instructions to run the CCA
      * Start a session in R version > 4.0.1\n
      * Difene the folder 'R_CCA' as your working directory\n
      * Open the file 'main.R' located in the folder 'code'\n
      * Excute all commands in 'main.R'\n
      The results will be saved in a folder named 'output_CCA'."
    write.table(readme, file = paste(new_path, "READ-ME.txt", sep = "/"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  }
  
  path_fig <- newDir("figures", at_root = FALSE, root = new_path)
  
  cat("Generating plots ...", "\n")
  
  p1 <- plotCor(cca$pearson_cor$corX, "X", heatmap_clust_method, heatmap_color_ramp, heatmap_border_color, heatmap_fontsize)
  savePlot(p1$p, "[p1]heatmap_matrix-of-pearson-correlations-for-variables-from-set-X", heatmap_width, heatmap_height, plots_file_extension, dir_name = "pearson-correlations", at_root = FALSE, root = path_fig)
  
  p2 <- plotCor(cca$pearson_cor$corY, "Y", heatmap_clust_method, heatmap_color_ramp, heatmap_border_color, heatmap_fontsize)
  savePlot(p2$p, "[p2]heatmap_matrix-of-pearson-correlations-for-variables-from-set-Y", heatmap_width, heatmap_height, plots_file_extension, dir_name = "pearson-correlations", at_root = FALSE, root = path_fig)
  
  p3 <- plotCor(cca$pearson_cor$corXY, "XY", heatmap_clust_method, heatmap_color_ramp, heatmap_border_color, heatmap_fontsize)
  savePlot(p3$p, "[p3]heatmap_matrix-of-pearson-correlations-for-variables-from-set-XY", heatmap_width, heatmap_height, plots_file_extension, dir_name = "pearson-correlations", at_root = FALSE, root = path_fig)
  
  p4 <- plotCorNetwork(cca$pearson_cor$corX, "X", cca = NULL, plotCorCCA_corr_min, plotCorCCA_seed, plotCorCCA_node_size, plotCorCCA_label_size, plots_x_color,
    plots_y_color, plots_cv_color, plots_neg_color, plots_pos_color, plotCorCCA_net_method)
  savePlot(p4$p, "[p4]graph_matrix-of-pearson-correlations-for-variables-from-set-X", plotCorCCA_width, plotCorCCA_height, plots_file_extension, dir_name = "pearson-correlations", at_root = FALSE, root = path_fig)
  
  p5 <- plotCorNetwork(cca$pearson_cor$corY, "Y", cca = NULL, plotCorCCA_corr_min, plotCorCCA_seed, plotCorCCA_node_size, plotCorCCA_label_size, plots_x_color,
    plots_y_color, plots_cv_color, plots_neg_color, plots_pos_color, plotCorCCA_net_method)
  savePlot(p5$p, "[p5]graph_matrix-of-pearson-correlations-for-variables-from-set-Y", plotCorCCA_width, plotCorCCA_height, plots_file_extension, dir_name = "pearson-correlations", at_root = FALSE, root = path_fig)
  
  p6 <- plotCorNetwork(cca$pearson_cor$corXY, "XY", cca, plotCorCCA_corr_min, plotCorCCA_seed, plotCorCCA_node_size, plotCorCCA_label_size, plots_x_color,
    plots_y_color, plots_cv_color, plots_neg_color, plots_pos_color, plotCorCCA_net_method)
  savePlot(p6$p, "[p6]graph_matrix-of-pearson-correlations-for-variables-from-set-XY", plotCorCCA_width, plotCorCCA_height, plots_file_extension, dir_name = "pearson-correlations", at_root = FALSE, root = path_fig)
  
  p7 <- plotCCA(cca, plots_canonical_variate, plotCCA_corr_min, plots_x_title, plots_y_title, plots_short_x_title, plots_short_y_title, plots_x_color, plots_y_color, plots_cv_color)
  savePlot(p7$p, "[p7]2Dplane_correlations-for-variables-with-canonical-variables", plotCCA_width, plotCCA_height, plots_file_extension, dir_name = "cca", at_root = FALSE, root = path_fig)

  corr_cca <- corCCA(cca, plotCorCCA_n_canonical_variates, plots_canonical_variate, plots_short_x_title, plots_short_y_title)
  p8 <- plotCorNetwork(corr_cca, "all", cca, plotCorCCA_corr_min, plotCorCCA_seed, plotCorCCA_node_size, plotCorCCA_label_size, plots_x_color,
    plots_y_color, plots_cv_color, plots_neg_color, plots_pos_color, plotCorCCA_net_method)
  savePlot(p8$p, "[p8]graph_matrix-of-pearson-correlations-for-variables-from-set-XandYandCV", plotCorCCA_width, plotCorCCA_height, plots_file_extension, dir_name = "cca", at_root = FALSE, root = path_fig)

  p9 <- plotCoef(cca, plots_canonical_variate, plots_x_title, plots_y_title, plots_pos_color, plots_neg_color) 
  savePlot(p9$p, "[p9]raw-estimated-coefficients-of-canonical-variable", plotEst_width, plotEst_height, plots_file_extension, dir_name = "cca-extra-figures", at_root = FALSE, root = path_fig)
  
  p10 <- plotStruct(cca, plots_canonical_variate, plots_x_title, plots_x_title, plots_pos_color, plots_neg_color, plotCCA_corr_min)
  savePlot(p10$p, "[p10]1Dlines_correlations-for-variables-with-canonical-variables", plotEst_width, plotEst_height, plots_file_extension, dir_name = "cca-extra-figures", at_root = FALSE, root = path_fig)
  
  p11 <- plotCumVarExpl(cca, plots_x_title, plots_y_title, plots_x_color, plots_y_color)
  savePlot(p11$p, "[p11]cumulative-variance-explained-by-canonical-variables", plotCumVar_width, plotCumVar_height, plots_file_extension, dir_name = "cca-extra-figures", at_root = FALSE, root = path_fig)
  
  p12 <- plotCorrel(cca, plots_pos_color, plots_neg_color, rCCA_regularization)
  savePlot(p12$p, "[p12]canonical-correlation-between-pairs-of-canonical-variates", plotCC_width, plotCC_height, plots_file_extension, dir_name = "cca-extra-figures", at_root = FALSE, root = path_fig)
  
  p13 <- plotScores(cca, plots_canonical_variate, ID, GROUP, plots_group_color, plots_cv_color, plots_x_title, plots_y_title, plots_short_x_title, plots_short_y_title)
  savePlot(p13$p, "[p13]scores-of-individuals-on-canonical-variates", plotCCA_width, plotCCA_height, plots_file_extension, dir_name = "cca", at_root = FALSE, root = path_fig)
  
  leg_fig <- figLegends(plots_x_title, plots_y_title, plots_x_color, plots_y_color, plots_cv_color, plots_pos_color, plots_neg_color, plotCorCCA_corr_min, plotCCA_corr_min, plots_canonical_variate, plotCorCCA_n_canonical_variates, rCCA_regularization)
  write.table(leg_fig, file = paste(new_path, "figures", "all_legends.txt", sep = "/"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  leg_tables <- tableLegends(cca_tables, plots_x_title, plots_y_title)
  write.table(leg_tables, file = paste(new_path, "tables", "all_legends.txt", sep = "/"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  cat("Calculating statistics ...", "\n")
  
  desc <- descMeasures(cca$cca_data$X, cca$cca_data$Y, print_descriptive_statistics)
  saveTable(FALSE, desc$desc_X$Descriptives, "description_for_X", dir_name = "descriptive-analysis", at_root = FALSE, root = new_path)
  saveTable(FALSE, desc$desc_Y$Descriptives, "description_for_Y", dir_name = "descriptive-analysis", at_root = FALSE, root = new_path)
  message("If you have chosen to transform your data (log, sqrt, scale), consider descriptive statistics calculated for the transformed data.")
  
  cat("CCA is finished. Results are in path: ", new_path, "\n")
  
  output <- list(
    CCA = cca,
    tables = cca_tables,
    plots = lapply(list(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13), "[[", 1),
    data_plots = lapply(list(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13), "[[", 2),
    stat_desc = desc,
    fig_legends = leg_fig,
    table_legends = leg_tables)
  
  names(output$plots) <- paste0("p", 1:13)
  names(output$data_plots) <- paste0("data_p", 1:13)
  
  return(output)
}

possiblyCCA <- purrr::possibly(routineCCA, "error")

runCCA <- function(...) {
  tryCCA <- possiblyCCA(...)
  closeAllConnections()
  
  message("Moving .log files ...")
  
  output_dir <- paste(getwd(), "output_CCA", sep = "/")
  dirs <- file.info(list.dirs(output_dir, full.names = TRUE))
  recent_dir <- rownames(dirs)[which.max(dirs$mtime)]
  
  file.move(paste(output_dir, "messages.log", sep = "/"), recent_dir)
  file.move(paste(output_dir, "outputs.log", sep = "/"), recent_dir)
  
  if (is.character(tryCCA)) {
    dir_new <- gsub(":| ", "-", format(Sys.time(), "%Y-%m-%d_%X"))
    new_path <- newDir(dir_new, at_root = FALSE, root = output_dir)
    
    file.move(paste(recent_dir, "messages.log", sep = "/"), new_path)
    file.move(paste(recent_dir, "outputs.log", sep = "/"), new_path)
    
    error <- paste0("The proccess was aborted due to errors. Look-up the files 'messages.log' and 'outputs.log' in:\n", new_path, ".")
    write.table(error, file = paste(new_path, "ERROR.log", sep = "/"), quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    unlink(recent_dir, recursive = TRUE)
    
    message("\n", "Proccess has failed.", "\n", "Look-up the files in ", new_path, ".", "\n")
  }
  
  return(tryCCA)
}