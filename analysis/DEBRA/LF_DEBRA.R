EV1 <- read.csv("/Users/feiliyang/Desktop/published_analysis_tools/DEBRA/paper_data/msb199195-sup-0004-datasetev1.csv", 
                row.names = 1)
EV1 <- as.matrix(EV1)
is.numeric(EV1)
plot(density(log2(EV1 + 1)))
dim(EV1)
colSums(EV1)
ggplot(data.frame(libSize = colSums(EV1)[1:14],
                  sample = colnames(EV1)[1:14])) + 
  geom_point(aes(x = seq_len(14), y = libSize)) +
  geom_text(aes(x = seq_len(14), y = libSize, label = sample), 
            vjust = -1,
            size = 3,
            angle = 45)

ggplot(data.frame(P1 = log2(EV1[, "P1"] +1))) + 
  geom_histogram(aes(x=P1))

ggplot(data.frame(P2 = log2(EV1[, "P2"] +1))) + 
  geom_histogram(aes(x=P2))

ggplot(data.frame(mix = log2(rowSums(EV1[, c("P1", "P2")]) +1))) + 
  geom_histogram(aes(x=mix))

samples_1 <- c("null_660.1","null_660.2","null_330.1","null_330.2","m_null_20.p35.1","m_null_20.p35.2")
EV1_s1 <- EV1[, samples_1]
ggplot(data.frame(libSize = colSums(EV1_s1),
                  sample = colnames(EV1_s1))) + 
  geom_point(aes(x = seq_len(6), y = libSize)) +
  geom_text(aes(x = seq_len(6), y = libSize, label = sample), 
            vjust = -1,
            size = 3,
            angle = 45)

install_github("Oshlack/Barbie")
library(Barbie)
EV1_s1_bb <- Barbie::createBarbie(
  object = EV1_s1,
  target = data.frame(
    perturb = c(rep("null", 4), rep("p35", 2)),
    libSize = rep(c(660, 330, 20), each=2)))
EV1_s1_bb <- Barbie::tagTopBarcodes(Barbie = EV1_s1_bb)
EV1_s1_bb$isTop$vec %>% sum()
Barbie::plotBarcodePareto(Barbie = EV1_s1_bb)
Barbie::plotBarcodeSankey(Barbie = EV1_s1_bb)
EV1_s1_bb <- Barbie::subsetBarcodes(Barbie = EV1_s1_bb, retainedRows = EV1_s1_bb$isTop$vec)
Barbie::plotBarbieHeatmap(Barbie = EV1_s1_bb)
Barbie::plotSamplePairCorrelation(Barbie = EV1_s1_bb)
EV1_s1_bb <- Barbie::testBarcodeBias(Barbie = EV1_s1_bb)
Barbie::plotBarcodeBiasHeatmap(Barbie = EV1_s1_bb)
EV1_s1_bb$testBarcodes$diffProp_perturb$methods
Barbie::plotBarcodeBiasScatterPlot(Barbie = EV1_s1_bb, pValuesAdjusted = FALSE, elementName = "diffProp_perturb")
results1 <- EV1_s1_bb$testBarcodes$diffProp_perturb$results %>% 
  dplyr::filter(p.value < 0.05)

EV1_s1_bb <- Barbie::testBarcodeBias(Barbie = EV1_s1_bb, method = "diffOcc")
Barbie::plotBarcodeBiasHeatmap(Barbie = EV1_s1_bb)
EV1_s1_bb$testBarcodes$diffOcc_perturb$methods
Barbie::plotBarcodeBiasScatterPlot(Barbie = EV1_s1_bb, pValuesAdjusted = FALSE)
results <- EV1_s1_bb$testBarcodes$diffOcc_perturb$results %>% 
  dplyr::filter(direction != "n.s.")

EV1_null <- EV1[, 1:14]
colnames(EV1_null)
EV1_null_bb <- Barbie::createBarbie(
  object = EV1_null,
  target = data.frame(
    perturb = c("P1", "P2", rep("null", 12)),
    libSize = c("P1", "P2", rep(c(660, 330, 160, 80, 40, 20), each=2))
  )
)

pca <- prcomp(t(EV1_null_bb$CPM))
v1v2 <- pca$x[, 1:2]
v1v2 <- cbind(v1v2, EV1_null_bb$metadata)
ggplot(v1v2) +
  geom_point(aes(x=PC1, y=PC2, color=libSize)) +
  theme_classic()

EV1_null_bb <- Barbie::subsetSamples(
  Barbie = EV1_null_bb, retainedColumns = c(F,F,rep(T, 12)))

pca <- prcomp(t(EV1_null_bb$CPM))
v1v2 <- pca$x[, 1:2]
v1v2 <- cbind(v1v2, EV1_null_bb$metadata)
ggplot(v1v2) +
  geom_point(aes(x=PC1, y=PC2, color=libSize)) +
  theme_classic() + 
  theme(aspect.ratio = 1)

set.seed(42) 

Sample_Groups <- function(Barbie, size_group = 3) {
  
  pool <- Barbie$metadata %>% nrow() %>% seq()
  
  # First sampling without replacement
  group1_samples <- sample(pool, size = size_group, replace = FALSE)
  # Remove the first sample from the pool
  remaining_pool <- setdiff(pool, group1_samples)
  
  # Second sampling without replacement from the remaining pool
  group2_samples <- sample(remaining_pool, size = size_group, replace = FALSE)
  
  # Create a group vector
  group_vec <- vector("character", length = length(pool))
  group_vec[group1_samples] <- "group1"
  group_vec[group2_samples] <- "group2"
  
  Barbie$metadata$group <- group_vec
  
  # Hem data by group vector
  sample_bb <- Barbie::subsetSamples(
    Barbie = Barbie, 
    retainedColumns = group_vec %in% c("group1", "group2")) 
  
  
  mytargets <- data.frame(sample_bb$metadata)
  
  # model design matrix
  mydesign <- mytargets %>%
    with(model.matrix(~0 + group))
  
  # set up block
  myblock <- NULL
  
  # apply test for Bias_Prop
  sample_bb <- Barbie::testBarcodeBias(
    Barbie = sample_bb,
    targets = mytargets,
    designMatrix = mydesign,
    sampleGroups = "group"
  )
  
  #penalized logistic madel
  sample_bb <- Barbie::testBarcodeBias(
    Barbie = sample_bb, 
    method = "diffOcc",
    regularization = "firth"
    )
  
  # extract p.values and FPR
  pvalues_occ <- sample_bb$TestBarcodes$diffProp_group$results
  pvalues_prop <- sample_bb$Bias_Prop$pvalue
  pvalues_LR_occ_classic <- sample_bb$LR_Occ_classic$p.values
  pvalues_LR_occ_penalized <- sample_bb$LR_Occ_penalized$p.values
  
  FPR_occ <- sum(sample_bb$Bias_Occ$group != "Unbiased") / length(sample_bb$Bias_Occ$group)
  
  FPR_prop <- sum(sample_bb$Bias_Prop$group != "Unbiased") / length(sample_bb$Bias_Prop$group)
  
  FPR_LR_occ_classic <- sum(sample_bb$LR_Occ_classic$group != "Unbiased") / length(sample_bb$LR_Occ_classic$group)
  
  FPR_LR_occ_penalized <- sum(sample_bb$LR_Occ_penalized$group != "Unbiased") / length(sample_bb$LR_Occ_penalized$group)
  
  return(
    list(FPR_occ = FPR_occ,
         FPR_prop = FPR_prop,
         FPR_LR_occ_classic = FPR_LR_occ_classic,
         FPR_LR_occ_penalized = FPR_LR_occ_penalized,
         pvalues_occ = pvalues_occ,
         pvalues_prop = pvalues_prop,
         pvalues_LR_occ_classic = pvalues_LR_occ_classic,
         pvalues_LR_occ_penalized = pvalues_LR_occ_penalized,
         group_vec = group_vec)
  )
  
}
