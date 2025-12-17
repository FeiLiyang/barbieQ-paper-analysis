# ---------- Negative Simulation ----------
# function for sampling once making two equal groups
Sample_Groups <- function(barbieQ, size_group = 3) {

  pool <- ncol(barbieQ) %>% seq()

  repeat{
    # First sampling without replacement
    X1_samples <- sample(pool, size = size_group, replace = FALSE)
    # Remove the first sample from the pool
    remaining_pool <- setdiff(pool, X1_samples)

    # Second sampling without replacement from the remaining pool
    X2_samples <- sample(remaining_pool, size = size_group, replace = FALSE)

    # Create a group vector
    group_vec <- vector("character", length = length(pool))
    group_vec[X1_samples] <- "X1"
    group_vec[X2_samples] <- "X2"

    barbieQ$sampleMetadata$Group <- group_vec

    # Hem data by group vector
    sample_bb <- barbieQ[, group_vec %in% c("X1", "X2")]

    # sample_bb$sampleMetadata$Time <- ifelse(
    #   sample_bb$sampleMetadata$Timepoint %in% c("TP1", "TP2"),
    #   "early", "late")

    # repeat if only one level in Time
    # if (length(unique(sample_bb$sampleMetadata$Time)) < 2L) next

    # model design matrix
    mydesign <- sample_bb$sampleMetadata %>%
      as.data.frame() %>%
      with(model.matrix(~0 + Group))
    # repeat sampling if mydesign is not full rank
    # case1: only one level in LibSize; case2: LibSize perfectly in line with grouping
    # break and use the sampling when design is full rank
    if (qr(mydesign)$rank == ncol(mydesign)) break
  }

  # apply test for Bias_Prop
  prop_asin <- testBarcodeSignif(
    barbieQ = sample_bb,
    sampleMetadata = sample_bb$sampleMetadata[,c("Group"),drop = FALSE],
    sampleGroup = "Group", designMatrix = mydesign,
    contrastFormula = "GroupX1 - GroupX2")

  prop_logit <- testBarcodeSignif(
    barbieQ = sample_bb,
    sampleMetadata = sample_bb$sampleMetadata[,c("Group"),drop = FALSE],
    sampleGroup = "Group", designMatrix = mydesign,
    contrastFormula = "GroupX1 - GroupX2", transformation = "logit")

  prop_noTrans <- testBarcodeSignif(
    barbieQ = sample_bb,
    sampleMetadata = sample_bb$sampleMetadata[,c("Group"),drop = FALSE],
    sampleGroup = "Group", designMatrix = mydesign,
    contrastFormula = "GroupX1 - GroupX2", transformation = "none")

  occ_firth <- testBarcodeSignif(
    barbieQ = sample_bb,
    sampleMetadata = sample_bb$sampleMetadata[,c("Group"),drop = FALSE],
    sampleGroup = "Group", designMatrix = mydesign,
    contrastFormula = "GroupX1 - GroupX2", method = "diffOcc" )

  # occurrence <- assays(sample_bb)$occurrence
  # results <- lapply(seq_len(nrow(occurrence)), function(i) {
  #   df <- data.frame(response = as.numeric(occurrence[i,]), mydesign)
  #   invisible(logistf::logistf(occurrence[i, ] ~ GroupX1, data = df))
  # })
  # ## extract stats
  # P.Value <- lapply(results, function(x) {
  #   x$prob["GroupX1"]
  # }) %>%
  #   unlist()
  #
  # adj.P.Value <- stats::p.adjust(P.Value, method = "BH")

  stats_prop_asin <- rowData(prop_asin)$testingBarcode %>% as.data.frame()
  stats_prop_logit <- rowData(prop_logit)$testingBarcode %>% as.data.frame()
  stats_prop_noTrans <- rowData(prop_noTrans)$testingBarcode %>% as.data.frame()
  stats_occ_firth <- rowData(occ_firth)$testingBarcode %>% as.data.frame()

  FPR_prop_asin <- sum(stats_prop_asin$P.Value < 0.05) / nrow(sample_bb)
  FPR_prop_logit <- sum(stats_prop_logit$P.Value < 0.05) / nrow(sample_bb)
  FPR_prop_noTrans <- sum(stats_prop_noTrans$P.Value < 0.05) / nrow(sample_bb)
  FPR_occ_firth <- sum(stats_occ_firth$P.Value < 0.05) / nrow(sample_bb)

  Diff_prop <- stats_prop_noTrans$meanDiff
  Amean_prop <- stats_prop_noTrans$Amean

  sample_X1 <- sample_bb[,sample_bb$sampleMetadata$Group == "X1"]
  sample_X2 <- sample_bb[,sample_bb$sampleMetadata$Group == "X2"]
  Diff_occ <- rowSums(assays(sample_X1)$occurrence) / ncol(sample_X1) - rowSums(assays(sample_X2)$occurrence) / ncol(sample_X2)
  Amean_occ <- rowSums(assays(sample_bb)$occurrence) / ncol(sample_bb)

  LFC_prop <- rowMeans(log2(assays(sample_X1)$CPM + 0.5)) - rowMeans(log2(assays(sample_X2)$CPM + 0.5))
  Gmean_prop <- rowMeans(log2(assays(sample_bb)$CPM + 0.5))

  return(
    list(
      P.Val_prop_asin = stats_prop_asin$P.Value,
      P.Val_prop_logit = stats_prop_logit$P.Value,
      P.Val_prop_noTrans = stats_prop_noTrans$P.Value,
      P.Val_occ_firth = stats_occ_firth$P.Value,

      FPR_prop_asin = FPR_prop_asin,
      FPR_prop_logit = FPR_prop_logit,
      FPR_prop_noTrans = FPR_prop_noTrans,
      FPR_occ_firth = FPR_occ_firth,

      Diff_prop = Diff_prop,
      Amean_prop = Amean_prop,
      LFC_prop = LFC_prop,
      Gmean_prop = Gmean_prop,
      Diff_occ = Diff_occ,
      Amean_occ = Amean_occ,

      group_vec = group_vec)
  )

}

## run 100 loops for certain group size
Loop_Samples <- function(barbieQ, size_group = 3, loop_times = 100) {

  FPR_prop_asin_list <- list(length = loop_times)
  FPR_prop_logit_list <- list(length = loop_times)
  FPR_prop_noTrans_list <- list(length = loop_times)
  FPR_occ_firth_list <- list(length = loop_times)

  P.Val_prop_asin_list <- list(length = loop_times)
  P.Val_prop_logit_list <- list(length = loop_times)
  P.Val_prop_noTrans_list <- list(length = loop_times)
  P.Val_occ_firth_list <- list(length = loop_times)

  Diff_prop_list <- list(length = loop_times)
  Amean_prop_list <- list(length = loop_times)
  LFC_prop_list <- list(length = loop_times)
  Gmean_prop_list <- list(length = loop_times)
  Diff_occ_list <- list(length = loop_times)
  Amean_occ_list <- list(length = loop_times)

  group_vec_list <- list(length = loop_times)

  for(i in 1:loop_times) {
    results_i <- Sample_Groups(
      barbieQ = barbieQ, size_group = size_group)

    FPR_prop_asin_list[[i]] <- results_i$FPR_prop_asin
    FPR_prop_logit_list[[i]] <- results_i$FPR_prop_logit
    FPR_prop_noTrans_list[[i]] <- results_i$FPR_prop_noTrans
    FPR_occ_firth_list[[i]] <- results_i$FPR_occ_firth

    P.Val_prop_asin_list[[i]] <- results_i$P.Val_prop_asin
    P.Val_prop_logit_list[[i]] <- results_i$P.Val_prop_logit
    P.Val_prop_noTrans_list[[i]] <- results_i$P.Val_prop_noTrans
    P.Val_occ_firth_list[[i]] <- results_i$P.Val_occ_firth

    Diff_prop_list[[i]] <- results_i$Diff_prop
    Amean_prop_list[[i]] <- results_i$Amean_prop
    LFC_prop_list[[i]] <- results_i$LFC_prop
    Gmean_prop_list[[i]] <- results_i$Gmean_prop
    Diff_occ_list[[i]] <- results_i$Diff_occ
    Amean_occ_list[[i]] <- results_i$Amean_occ

    group_vec_list[[i]] <- results_i$group_vec
  }

  FPR <- data.frame(
    Prop_asin = FPR_prop_asin_list %>% unlist(),
    Prop_logit = FPR_prop_logit_list %>% unlist(),
    Prop_noTrans = FPR_prop_noTrans_list %>% unlist(),
    Occ_firth = FPR_occ_firth_list %>% unlist(),
    N_Samples = size_group,
    Loop_N = 1:loop_times
  )
  FPR$SimulationID <- paste0("sample_", FPR$N_Samples, "_loop_", FPR$Loop_N)

  P.Val <- data.frame(
    Prop_asin = P.Val_prop_asin_list %>% unlist(),
    Prop_logit = P.Val_prop_logit_list %>% unlist(),
    Prop_noTrans = P.Val_prop_noTrans_list %>% unlist(),
    Occ_firth = P.Val_occ_firth_list %>% unlist(),
    N_Samples = size_group,
    Loop_N = rep(1:loop_times, each = nrow(barbieQ))
  )
  P.Val$SimulationID <- paste0("sample_", P.Val$N_Samples, "_loop_", P.Val$Loop_N)

  MA <- data.frame(
    Diff_prop = Diff_prop_list %>% unlist(),
    Amean_prop = Amean_prop_list %>% unlist(),
    LFC_prop = LFC_prop_list %>% unlist(),
    Gmean_prop = Gmean_prop_list %>% unlist(),
    Diff_occ = Diff_occ_list %>% unlist(),
    Amean_occ = Amean_occ_list %>% unlist(),
    N_Samples = size_group,
    Loop_N = rep(1:loop_times, each = nrow(barbieQ))
  )
  MA$SimulationID <- paste0("sample_", MA$N_Samples, "_loop_", MA$Loop_N)

  Group_Vec <- do.call(rbind, group_vec_list)
  rownames(Group_Vec) <- paste0("sample_", size_group, "_loop_", 1:loop_times)
  colnames(Group_Vec) <- colnames(barbieQ)

  return(list(
    FPR = FPR,
    P.Val = P.Val,
    MA = MA,
    Group_Vec = Group_Vec))

}

# run the loop for each size of samples
random_sampling <- function(barbieQ, loop_times = 100) {

  num_samples <- ncol(barbieQ)
  end_sampling <- floor(num_samples / 2)

  # stop running loop if end_sampling <= 3
  if(num_samples < 3) {stop("not enough samples!")}

  set.seed(2025) # for reproducibility

  # run the loops
  all_loops <- lapply(3:end_sampling, function(n) {
    Loop_Samples(
      barbieQ = barbieQ, size_group = n, loop_times = loop_times)
  })
  names(all_loops) <- paste0("N_Sample_", 3:end_sampling)

  # save it globally
  global_all_loops <<- all_loops

  ## default wd: "../../analysis/DEBRA"
  # save(all_loops, file = "../../output/DEBRA/all_loops.rda")

  # extract results from the loops
  ## extract FPR
  all_FPR <- lapply(seq_along(3:end_sampling), function(n) {
    all_loops[[n]]$FPR
  })
  df_FPR <- do.call(rbind, all_FPR)
  df_FPR$N_Samples <- as.factor(df_FPR$N_Samples)
  ## extract P.Val
  all_P.Val <- lapply(seq_along(3:end_sampling), function(n) {
    all_loops[[n]]$P.Val
  })
  df_P.Val <- do.call(rbind, all_P.Val)
  df_P.Val$N_Samples <- as.factor(df_P.Val$N_Samples)

  gFDR_prop_asin <- ggplot(as.data.frame(df_FPR),
                           aes(x = N_Samples, y = Prop_asin, fill = N_Samples)) +
    # geom_violin(width = 0.5) +  # Set width for better visibility +
    geom_boxplot(width = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.1, height = 0, size = 1) +  # Add jittered points
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    labs(title = "Diff_Prop_asin",
         y = "Prop. barcodes with p.value < 0.05",
         x = "No. samples per group") +
    theme_classic()

  gFDR_prop_logit <- ggplot(as.data.frame(df_FPR),
                            aes(x = N_Samples, y = Prop_logit, fill = N_Samples)) +
    # geom_violin(width = 0.5) +  # Set width for better visibility
    geom_boxplot(width = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.1, height = 0, size = 1) +  # Add jittered points
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    labs(title = "Diff_Prop_logit",
         y = "Prop. barcodes with p.value < 0.05",
         x = "No. samples per group") +
    theme_classic()

  gFDR_prop_noTrans <- ggplot(as.data.frame(df_FPR),
                            aes(x = N_Samples, y = Prop_noTrans, fill = N_Samples)) +
    # geom_violin(width = 0.5) +  # Set width for better visibility
    geom_boxplot(width = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.1, height = 0, size = 1) +  # Add jittered points
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    labs(title = "Diff_Prop_noTrans",
         y = "Prop. barcodes with p.value < 0.05",
         x = "No. samples per group") +
    theme_classic()

  gFDR_occ_firth <- ggplot(as.data.frame(df_FPR),
                           aes(x = N_Samples, y = Occ_firth, fill = N_Samples)) +
    # geom_violin(width = 0.5) +  # Set width for better visibility +
    geom_boxplot(width = 0.5, outlier.shape = NA) +
    geom_jitter(width = 0.1, height = 0, size = 1) +  # Add jittered points
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    labs(title = "Diff_Occ_firth",
         y = "Prop. barcodes with p.value < 0.05",
         x = "No. samples per group") +
    theme_classic()


  gP_prop_asin <- ggplot(as.data.frame(df_P.Val),
                         aes(x = N_Samples, y = Prop_asin, fill = N_Samples)) +
    geom_violin(width = 0.5) +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    labs(title = "Diff_Prop_asin",
         y = "p.value")

  gP_prop_logit <- ggplot(as.data.frame(df_P.Val),
                          aes(x = N_Samples, y = Prop_logit, fill = N_Samples)) +
    geom_violin() +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    labs(title = "Diff_Prop_logit",
         y = "p.value")

  gP_prop_noTrans <- ggplot(as.data.frame(df_P.Val),
                          aes(x = N_Samples, y = Prop_noTrans, fill = N_Samples)) +
    geom_violin() +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    labs(title = "Diff_Prop_noTrans",
         y = "p.value")

  gP_occ_firth <- ggplot(as.data.frame(df_P.Val),
                         aes(x = N_Samples, y = Occ_firth, fill = N_Samples)) +
    geom_violin() +
    geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") +
    labs(title = "Diff_Occ_firth",
         y = "p.value")

  return(
    list(gFDR_prop_asin = gFDR_prop_asin,
         gFDR_prop_logit = gFDR_prop_logit,
         gFDR_prop_noTrans = gFDR_prop_noTrans,
         gFDR_occ_firth = gFDR_occ_firth,
         gP_prop_asin = gP_prop_asin,
         gP_prop_logit = gP_prop_logit,
         gP_prop_noTrans = gP_prop_noTrans,
         gP_occ_firth = gP_occ_firth))

}

get_random_sampling_loops <- function(barbieQ, loop_times = 100) {
  num_samples <- ncol(barbieQ)
  end_sampling <- floor(num_samples / 2)
  
  # stop running loop if end_sampling <= 3
  if(num_samples < 3) {stop("not enough samples!")}
  
  set.seed(2025) # for reproducibility
  
  # run the loops
  all_loops <- lapply(3:end_sampling, function(n) {
    Loop_Samples(
      barbieQ = barbieQ, size_group = n, loop_times = loop_times)
  })
  names(all_loops) <- paste0("N_Sample_", 3:end_sampling)
  
  return(all_loops)

}

# ---------- Estimate vs Design ----------

calc_theory <- function(a, trans = "none"){
  theory_pool1 = (1/2 + a) / (1+a)
  theory_pool2 = (1/2) / (1+a)
  
  trans_fun = ifelse(trans == "none", function(x) {x},
                     ifelse(trans == "asin-sqrt", function(x) {asin(sqrt(x))}, 
                            ifelse(trans == "logit", function(x) {log((x)/(1-x))}, NA))) ## {log((x+1e-10)/(1-x+1e-10))},or {logit(x)}
  
  # the ratio of a barcode in pool1 or pool2; ~= pool1+ pool2
  X_Pool1 = pool1_vec
  X_Pool2 = pool2_vec
  # baseline sample (null sample) is equal mixture of pool1 and pool2
  X_0 = (X_Pool1 + X_Pool2) * 1/2
  # in perturbed sample, BC of pool1 and BC of pool2 show a certain fold change
  X_a = theory_pool1 * X_Pool1 + theory_pool2 * X_Pool2
  
  # difference between X_a and X_0 following transformation
  diff = trans_fun(X_a) - trans_fun(X_0)
  
  # theory_FC1 = (trans_fun(theory_pool1 * pool1_vec )) - (trans_fun(1/2 * pool1_vec ))
  # theory_FC2 = (trans_fun(theory_pool2 * pool2_vec )) - (trans_fun(1/2 * pool2_vec ))
  # df <- data.frame(diff_1 = theory_FC1, diff_2 = theory_FC2)
  
  ## FC
  # theory_FC1 = (trans(theory_pool1) + 1e-10) / (trans(1/2) + 1e-10)
  # theory_FC2 = (trans(theory_pool2) + 1e-10) / (trans(1/2) + 1e-10)
  
  df <- data.frame(diff = diff)
  
  return(df)
}

designed_diff <- function(a, trans = "none") {
  ratio_pool1 = (1/2 + a) / (1+a)
  ratio_pool2 = (1/2) / (1+a)
  
  ## an offset has been added to raw count before - no need to add offset in logit function
  trans_fun = ifelse(trans == "none", function(x) {x},
                     ifelse(trans == "asin-sqrt", function(x) {asin(sqrt(x))}, 
                            ifelse(trans == "logit", function(x) {log((x)/(1-x))}, NA)))
  
  ## baseline sample (equal mixture of barcodes)
  x0 = 1/2 * pool1_vec + 1/2 * pool2_vec
  
  
}

apply_model <- function(mixed, mixed_design, myblock, a = 0.18, trans = "none", scale = 1) {
  
  mixed@assays@data$proportion <- mixed@assays@data$proportion * scale
  
  tester <- testBarcodeSignif(
    barbieQ = mixed, sampleMetadata = mixed_target, sampleGroup = "Perturbation", 
    contrastFormula = paste0("Perturbation", a, " - Perturbation0"), transformation = trans,
    designMatrix = mixed_design, block = myblock)
  
  veri_diff <- cbind(
    tester@elementMetadata$testingBarcode,
    calc_theory(a = a, trans = trans)) %>% as.data.frame()
  
  veri_diff <- rownames_to_column(veri_diff, var = "BC") 
  
  veri_diff$Contrast = paste0(a, " - 0")
  veri_diff$Method = trans
  
  # veri_diff$uncertain_origin <- flag_small_gap
  # veri_diff$noise_in_origin <- flag_no_zero
  
  return(veri_diff)
  
}

get_residules <- function(mixed, mixed_design, myblock, a = 0.18, trans = "none") {
  
  ## get props
  props <- mixed@assays@data$proportion
  
  ## set trans fun
  trans_fun = ifelse(
    trans == "none", function(x) {x},
    ifelse(trans == "asin-sqrt", function(x) {asin(sqrt(x))}, 
           ifelse(trans == "logit", function(x) {log((x)/(1-x))}, NA))) ## {log((x+1e-10)/(1-x+1e-10))},or {logit(x)}
  
  mydata <- trans_fun(props)
  if(is.null(myblock)) {
    myfit1 <- limma::lmFit(
      object = mydata, design = mixed_design, block = myblock)
  } else {
    dup <- limma::duplicateCorrelation(object = mydata, design = mixed_design, block = myblock)
    myfit1 <- limma::lmFit(
      object = mydata, design = mixed_design, block = myblock, correlation = dup$consensus.correlation)
  }
  
  mycontrast <- limma::makeContrasts(contrasts = paste0("Perturbation", a, " - Perturbation0"), levels = colnames(mixed_design))
  
  ## fit contrast
  myfit2 <- limma::contrasts.fit(fit = myfit1, contrasts = mycontrast)
  
  ## fit eBayes, moderated p.values obtained after applying empirical Bayes
  myfit3 <- limma::eBayes(myfit2)
  
  ## extract residuals
  resi <- limma::residuals.MArrayLM(object = myfit1, y = mydata)
  
  resi_df <- as.data.frame(as.table(resi))
  
  colnames(resi_df) <- c("BC", "Sample", "Residual")
  
  resi_df$Contrast <- paste0(a, " - 0")
  
  resi_df$Method <- trans
  
  return(resi_df)
  
}

# -----------------Function to N_BC --------------
  
get_null_sample_by_filter <- function(barbieQ, filter_level = 0.95, loop_times = 3) {
  ## tag top BC from raw barbieQ
  barbieQ <- tagTopBarcodes(barbieQ, nSampleThreshold = 8, proportionThreshold = filter_level)
  ## select top BC from barbieQ_plus1
  barbieQ_top <- barbieQ[rowData(barbieQ)$isTopBarcode$isTop,]
  ## subset baseline "null" samples
  nulls <- barbieQ_top[, barbieQ_top$sampleMetadata$Subset == "null"]
  ## re-calculate proportions
  nulls@assays@data$proportion <- ((nulls@assays@data$proportion %>% t()) / colSums(nulls@assays@data$proportion)) %>% t()
  
  ## run random sampling loops
  result_loops <- suppressMessages({
    get_random_sampling_loops(barbieQ = nulls, loop_times = loop_times)
  })
  
  ## extract results from the loops
  end_sampling <- 6
  ## extract FPR
  all_FPR <- lapply(seq(3:end_sampling), function(n) {result_loops[[n]]$FPR})
  df_FPR <- do.call(rbind, all_FPR)
  ## extract P.Val
  all_P.Val <- lapply(seq(3:end_sampling), function(n) {result_loops[[n]]$P.Val})
  df_P.Val <- do.call(rbind, all_P.Val)
  ## extract MA
  all_MA <- lapply(seq(3:end_sampling), function(n) {result_loops[[n]]$MA})
  df_MA <- do.call(rbind, all_MA)
  
  all_Group_Vec <- lapply(seq(3:end_sampling), function(n) {result_loops[[n]]$Group_Vec})
  df_Group_Vec <- do.call(rbind, all_Group_Vec)
  
  stats_barcodes <- cbind(df_P.Val[,1:4], df_MA)
  
  return(stats_barcodes)
}

plot_fdr_on_quantile <- function(stats_barcodes_i, num_quantiles = 10) {
  
  quantile_breaks <- quantile(stats_barcodes_i$Amean_prop, probs = seq(0, 1, length.out = num_quantiles+1)) 
  
  ## for scientific format of x laebls for quantiles breaks
  bin_labels <- mapply(
    function(lo, hi) bquote("("*.(format_math(lo))*","*.(format_math(hi))*"]"),
    head(quantile_breaks, -1),
    tail(quantile_breaks, -1),
    SIMPLIFY = FALSE
  )
  
  ## specify label here
  stats_binned <- stats_barcodes_i %>%
    mutate(bin = cut(Amean_prop, breaks = quantile_breaks, labels = seq_len(num_quantiles), include.lowest = TRUE, right = TRUE)) %>% 
    pivot_longer(
      cols = c(Prop_asin, Prop_logit, Prop_noTrans),
      names_to = "method",
      values_to = "pvalue"
    )
  
  FPR_lines <- stats_binned %>%
    group_by(bin, method, N_Samples, Loop_N) %>%
    summarize(Frac_signif = sum(pvalue < 0.05) / n(), .groups = "drop") %>% 
    mutate(Method = recode(method, !!!method_labels))
  
  p <- ggplot(FPR_lines, aes(x = bin, y = Frac_signif, color = Method, group = paste0(bin, Method))) +
    geom_boxplot(alpha = 0.6, position = position_dodge(width = 0.8), outliers = T) +
    labs(
      x = "Quantiles of barcode mean proportion (raw)",
      y = "Fraction of P-Value < 0.05", color = "Method", shape = "N_sample", linetype = "N_sample"
      # caption = paste("Quantile mapping:\n", legend_text)
    ) +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    theme_linedraw() +
    scale_x_discrete(labels = scales::label_parse()) +
    # theme(axis.text.x = element_text(angle = 75, vjust = 0.5, size = 10, hjust = 0.5), aspect.ratio = 0.5) + 
    theme(
      legend.title = element_text(size = 12), legend.text = element_text(size = 11),
      axis.title = element_text(size = 13), axis.text = element_text(size = 12))
  
  return(p)
}

## Helper to format in "×10^" notation
format_math <- function(x) {
  exp <- floor(log10(x))
  base <- x / 10^exp
  # Return expression object (not a string!)
  bquote(.(formatC(base, digits = 2, format = "f")) %*% 10^.(exp))
}

get_mixted_by_filter <- function(barbieQ, filter_level = 0.95) {
  ## tag top BC from raw barbieQ
  barbieQ <- tagTopBarcodes(barbieQ, nSampleThreshold = 8, proportionThreshold = filter_level)
  ## select top BC from barbieQ_plus1
  barbieQ_top <- barbieQ[rowData(barbieQ)$isTopBarcode$isTop,]
  ## subset baseline "null" samples and perturbed samples
  mixed <- barbieQ_top[, barbieQ_top$sampleMetadata$Subset %in% c("null", "perturb")]
  
  ## re-calculate proportions
  mixed@assays@data$proportion <- ((mixed@assays@data$proportion %>% t()) / colSums(mixed@assays@data$proportion)) %>% t()
  
  mixed_target <- mixed$sampleMetadata %>% as.data.frame()
  mixed_target$Perturbation <- as.factor(mixed_target$Perturbation)
  # levels(mixed_target$Perturbation)
  mixed_target$BarcodeSize <- as.factor(mixed_target$BarcodeSize)
  # levels(mixed_target$BarcodeSize)
  
  mixed_design <- mixed_target %>%
    with(model.matrix(~0 + Perturbation + BarcodeSize))
  myblock <- mixed_target %>% with(paste(BarcodeSize, Perturbation))
  
  es_de_none <- lapply(c(0.18, 0.27, 0.35), function(i) test_perturb(mixed, mixed_target, mixed_design, NULL, a = i, trans = "none"))
  
  es_de_asin <- lapply(c(0.18, 0.27, 0.35), function(i) test_perturb(mixed, mixed_target, mixed_design, NULL, a = i, trans = "asin-sqrt"))
  
  es_de_logit <- lapply(c(0.18, 0.27, 0.35), function(i) test_perturb(mixed, mixed_target, mixed_design, NULL, a = i, trans = "logit"))
  
  es_de_df <- rbind(do.call(rbind, es_de_none), do.call(rbind, es_de_asin), do.call(rbind, es_de_logit))
  
  mean_prop <- rowMeans(mixed@assays@data$proportion)
  
  es_de_df <- es_de_df %>% 
    mutate(raw_Amean = rep(mean_prop, 9)) %>% 
    mutate(Signif. = ifelse(P.Value<0.05, "signif.", "n.s.")) %>% 
    mutate(Signif. = factor(Signif., levels = c("signif.", "n.s.")))
  
  es_de_df$filter_level = filter_level
  es_de_df$N_BC <- filter_to_N[filter_level %>% as.character()]
  
  return(es_de_df)
}

plot_tpr_on_quantile <- function(test_all_levels, num_quantiles = 10) {
  
  quantile_breaks <- quantile(test_all_levels$raw_Amean, probs = seq(0, 1, length.out = num_quantiles+1)) 
  
  ## for scientific format of x laebls for quantiles breaks
  bin_labels <- mapply(
    function(lo, hi) bquote("("*.(format_math(lo))*","*.(format_math(hi))*"]"),
    head(quantile_breaks, -1),
    tail(quantile_breaks, -1),
    SIMPLIFY = FALSE
  )
  
  ## specify label here
  stats_binned <- test_all_levels %>%
    mutate(bin = cut(raw_Amean, breaks = quantile_breaks, labels = seq_len(num_quantiles), include.lowest = TRUE, right = TRUE))
  
  FPR_lines <- stats_binned %>%
    group_by(bin, Method, Contrast) %>%
    summarize(Frac_signif = sum(P.Value < 0.05) / n(), .groups = "drop")
  
  p <- ggplot(FPR_lines, aes(x = bin, y = Frac_signif, color = Method, linetype = Contrast, shape = Contrast)) +
    # geom_boxplot(alpha = 0.6, position = position_dodge(width = 0.8), outliers = T) +
    geom_line(aes(group = paste0(Contrast, Method))) +
    geom_point() +
    labs(
      x = "Quantiles of barcode mean proportion (raw)",
      y = "Fraction of P-Value < 0.05", color = "Method", shape = "Contrast", linetype = "Contrast"
      # caption = paste("Quantile mapping:\n", legend_text)
    ) +
    geom_hline(yintercept = 0.05, linetype = "dashed") +
    theme_linedraw() +
    scale_x_discrete(labels = scales::label_parse()) +
    theme(axis.text.x = element_text(angle = 75, vjust = 0.5, size = 10, hjust = 0.5), aspect.ratio = 0.5) + 
    theme(
      legend.title = element_text(size = 12), legend.text = element_text(size = 11),
      axis.title = element_text(size = 13), axis.text = element_text(size = 12))
  
  return(p)
}
