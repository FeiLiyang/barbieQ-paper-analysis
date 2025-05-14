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
