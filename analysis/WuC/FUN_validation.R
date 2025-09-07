tagTop10 <- function(col, n=10) {
  ranks <- rank(-col)
  # summary(ranks)
  isTop <- ranks <= n
  # summary(isTop)
  return(isTop)
}

plot_line_validation <- function(barcode_set_list, barbieQ, get_long_data=FALSE) {
  ## Example for barcode set list
  # barcode_set_list <- list(
  #   test_30sample = rownames(top10each_30_test30)[top10each_30_test30@elementMetadata$testingBarcode$tendencyTo == "CelltypeNK_CD56n_CD16p"],
  #   test_early_mid = rownames(top10each_30_test_ealy_mid)[top10each_30_test_ealy_mid@elementMetadata$testingBarcode$tendencyTo == "CelltypeNK_CD56n_CD16p"],
  #   heuristic = union_LFC_top10each_30
  # )

  ## check barcode set list and barbieQ contents
  print(names(barcode_set_list))
  print(dim(barbieQ))

  # Create named vector with default value "none"
  selectedBy <- setNames(rep("none", nrow(barbieQ)), rownames(barbieQ))

  # Assign sources
  for (source_name in names(barcode_set_list)) {
    matching_barcodes <- barcode_set_list[[source_name]]
    selectedBy[matching_barcodes] <- ifelse(
      selectedBy[matching_barcodes] == "none",
      source_name,
      paste(selectedBy[matching_barcodes], source_name, sep = ",")
    )
  }

  ## extract props. and sample conditions
  prop. <- assays(barbieQ)$proportion
  cond. <- barbieQ$sampleMetadata %>% as.data.frame()

  cond.$Celltype <- factor(cond.$Celltype,
                           levels = c("T", "B", "Gr", "NK_CD56p_CD16n", "NK_CD56n_CD16p",
                                      "NK_NKG2Ap_CD16p", "NK_NKG2Ap_CD16p_KIR3DL01n",
                                      "NK_NKG2Ap_CD16p_KIR3DL01p"))

  ## combine slots
  prop_df <- data.frame(cond., t(prop.))
  ## keep only early and mid samples
  prop_df <- prop_df %>%
    # dplyr::filter(Phase %in% c("early", "mid")) %>%
    select(-Phase)
  ## classify BC by selected by test or LFC

  ## picot to long data
  prop_long <- prop_df %>%
    pivot_longer(
      cols = -c(Celltype, Months),     # keep Celltype and Months fixed
      names_to = "Barcode",
      values_to = "Proportion"
    )
  ## assign to barcode selectedBy status
  prop_long$selectedBy <- selectedBy[as.character(prop_long$Barcode)]

  if(get_long_data) {
    return(prop_long)
  }

  ## plot
  p <- ggplot(prop_long) +
    geom_line(aes(x = Months, y = Proportion, group = Barcode, color = selectedBy)) +
    facet_grid(selectedBy ~ Celltype, scales = "free_y") +
    theme_minimal() +
    theme(legend.position = "top")

  return(p)

}


plot_validation <- function(myObject) {

  flag <- rowData(myObject)$testingBarcode$direction == 1
  myObject <- myObject[flag,]

  NK_clone_contributionU <- assays(myObject)$proportion %>% colSums()
  dat <- data.frame(
    contr = NK_clone_contributionU,
    Celltype = factor(myObject$sampleMetadata$Celltype,
                      levels = c("T", "B", "Gr", "NK_CD56p_CD16n", "NK_CD56n_CD16p")),
    Months = myObject$sampleMetadata$Months)
  ## sort samples by celltype and months
  dat <- dat %>%
    arrange(Celltype, Months) %>%
    mutate(id = row_number())
  ## add dummy rows for spacing between celltypes
  dat <- dat %>%
    group_by(Celltype) %>%
    mutate(dummy = n() + 1) %>%    # Add a dummy count for spacing
    ungroup() %>%
    bind_rows(dat %>%             # Add dummy rows for spacing
                group_by(Celltype) %>%
                slice_tail(n = 1) %>%
                mutate(id = as.numeric(id) + 0.5, contr = NA, Months = NA)) %>%
    arrange(as.numeric(id))        # Arrange with dummy rows in place

  celltype_colors <- c(
    "T" = "yellowgreen",
    "B" = "yellowgreen",
    "Gr" = "yellowgreen",
    "NK_CD56p_CD16n" = "yellowgreen",
    "NK_CD56n_CD16p" = "slateblue1" # "orchid"
  )

  p <- ggplot(dat, aes(x = as.factor(id), y = contr, color = Celltype)) +
    geom_bar(stat = "identity", fill = NA) +
    theme_minimal() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 0.5, size = 10),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          panel.grid.major.x = element_blank()) +
    scale_x_discrete(labels = ifelse(is.na(dat$Months), "", dat$Months)) + # Relabel with `Celltype`
    scale_color_manual(values = celltype_colors) +  # <- set bar outline color
    scale_fill_manual(values = celltype_colors) +   # <- set tile fill color
    scale_y_continuous(limits = c(-0.09, 0.8)) +
    labs(x = "Samples of cell types at different months", y = "Total prop. of selected barcodes") +  # Update axis labels
    geom_tile(
      data = dat,
      aes(x = ifelse(is.na(Months), NA,order(id)), y = -0.03, fill = Celltype),
      height = 0.03, width = 1
    ) +
    geom_text(
      data = dat %>%
        mutate(id = order(id)) %>%
        group_by(Celltype) %>%
        mutate(Celltype = dplyr::recode(Celltype,"NK_CD56n_CD16p" = "NK.CD56-CD16+","NK_CD56p_CD16n" = "NK.CD56+CD16-")) %>%
        group_by(Celltype) %>%
        summarize(mid_point = mean(as.numeric(id))),
      aes(x = mid_point-0.5, y = -0.03, label = Celltype),
      size = 3, vjust = 0.5, color = "black"
    )

  return(p)
}

calc_bias_LFC <- function(barbieQ, get_FC=FALSE, get_Mean=FALSE, get_AMean=FALSE) {
  ## extract props. and sample conditions
  prop. <- assays(barbieQ)$proportion
  cond. <- barbieQ$sampleMetadata %>% as.data.frame()
  ## group samples by NKnp or others
  cond. <- cond. %>%
    mutate(isNK = ifelse(Celltype == "NK_CD56n_CD16p", "NKnp", "others"))
  cond.$Celltype %>% table()

  ## calc mean per group at each month
  group_labels <- paste(cond.$Months, cond.$isNK, sep = "__")
  group_means <- sapply(unique(group_labels), function(g) {
    idx <- which(group_labels == g)
    rowMeans(prop.[, idx, drop = FALSE])
  })

  if(get_Mean){
    return(group_means)
  }


  ## calc group LFC at each month
  colnames(group_means)
  months <- gsub("(\\1)__.*", "\\1",colnames(group_means))

  ## calc mean of all groups at each month
  group_Ameans <- sapply(unique(months), function(m) {
    idx <- which(months == m)
    rowMeans(prop.[, idx, drop = FALSE])
  })

  if(get_AMean){
    return(group_Ameans)
  }

  FC_months <- sapply(unique(months), function(m) {
    idx <- which(months == m)
    sub <- group_means[, idx, drop = FALSE]
    group <- gsub(".*__()", "\\1",colnames(sub))
    NK <- sub[, which(group == "NKnp"), drop = F]
    others <- sub[, which(group=="others"),drop = F]
    (NK+1e-6) / (others+1e-6)
  })

  if(get_FC){
    return(FC_months)
  }

  ## decide if each barcode LFC is high at each month
  is_LFC_high <- FC_months > 10
  colnames(is_LFC_high)

  ## calc if each barcode NKnp  prop > 0.01 at each month
  NK_prop <- prop.[,cond.$isNK == "NKnp", drop=F]
  is_NKnp_high <- NK_prop > 0.01
  colnames(is_NKnp_high)

  ## barcodes meet both metric at each month
  select_months <- (is_LFC_high & is_NKnp_high)

  ## select barcodes at each month
  select_barcodes_months <- lapply(colnames(select_months), function(j)
    rownames(prop.)[select_months[,j]])

  names(select_barcodes_months) <- paste0("heuristic at ", colnames(is_LFC_high), "m")

  return(select_barcodes_months)
}
