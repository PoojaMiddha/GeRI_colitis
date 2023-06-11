# LDPred2 on CD summary stats already matched with HapMap3 variants to obtain the best PRS 
# Approach used: Infinitesimal model and Auto model (Auto model will be preferred)
# Author: Pooja Middha
# Date: Oct 12, 2022

# Load library
library(data.table)
library(dplyr)
library(bigsnpr)
library(bigsparser)
library(ggplot2)
library(pROC)
library(stringr)

# Number of cores
NCORES <- nb_cores()

#info <- readRDS("/zivlab/data3/pmiddha/Data/ukbb/ldpred2/map_hm3_ldpred2.rds")
#str(info)

# Load UKB data for CD
obj.bigSNP <- snp_attach("/zivlab/data3/pmiddha/Proj_GeRI_GI-irAEs/PRS_UKBB/Data/genetic/CD/cd_ukb_allchr.rds")
str(obj.bigSNP, max.level = 2, strict.width = "cut")

# Read genotype matrix
G   <- obj.bigSNP$genotypes
big_counts(G, ind.col = 1:10)

# For some reason the UKB data despite being imputed has small proportion of NA
# Impute the missinng genptype because LDPred2 doesnt work with NA in genotype matrix
Gimp <- snp_fastImputeSimple(G, method = "mean2", ncores = NCORES)
big_counts(Gimp, ind.col = 1:10)

# Chromosome and Position 
CHR <- obj.bigSNP$map$chromosome
POS <- obj.bigSNP$map$physical.pos

# Phenotype data
pheno <- fread("/zivlab/data3/pmiddha/Proj_GeRI_GI-irAEs/PRS_UKBB/Data/pheno/gi_pheno_final_10112022.txt",
               header = T)

# ID is obtained from genetic data
ID <- obj.bigSNP$fam$family.ID

# Subset pheno data to only those IDs that are also in genetic data
subpheno <- subset(pheno, pheno$FID %in% ID)
dim(subpheno)
y <- subpheno$CD

# Import HapMap3 matched summary stats
sumstats <- readRDS("/zivlab/data3/pmiddha/Data/sumstats_delange_mrcieu/ldpred2_cd_sumstats_hapmapmatched.rds")
#colnames(sumstats)[7] <- "beta_se"

# Create validation and test dataset 
set.seed(1)
ind.val <- sample(nrow(Gimp), floor(0.70*nrow(Gimp))) # 70%
ind.test <- setdiff(rows_along(Gimp), ind.val) # 30%
print("Size of validation group")
length(ind.val)
print("Size of testing group")
length(ind.test)

# Extract main variables from target data (map section)
map <- setNames(obj.bigSNP$map[-3], c("chr", "rsid", "pos", "a1", "a0"))
head(map)

# Matching summary stats (base data) and target data
df_beta <- snp_match(sumstats, map)
dim(df_beta)

# Subset SNPs from summary stats after snp_match()
#sumstats_sub <- subset(sumstats, sumstats$rsid %in% df_beta$rsid)

# Obtain SNPs from matched summary stats which were reversed during the snp_match() part
#rev2 <- subset(sumstats_sub, sumstats_sub$a0 != df_beta$a0)

# Save reversed SNPs and the matched SNPs (df_beta)
#save(rev2, file = "/zivlab/data3/pmiddha/Data/sumstats_delange_mrcieu/cd_sumstats_reversed_ldpred2.Rdata")
#save(df_beta, file = "/zivlab/data3/pmiddha//Data/sumstats_delange_mrcieu/cd_sumstats_snpmatched_ldpred2.Rdata")

# Compute correlations between variants (Recommendation: 3cM window size)
tmp <- tempfile(tmpdir = "/zivlab/data3/pmiddha/Proj_GeRI_GI-irAEs/PRS_UKBB/Data/genetic/CD/tmp-data")
POS2 <- snp_asGeneticPos(CHR, POS, dir = "/zivlab/data3/pmiddha/Proj_GeRI_GI-irAEs/PRS_UKBB/Data/genetic/CD/tmp-data", 
                         ncores = NCORES)

for (chr in 1:22) {
  
  print(chr)
  
  ## indices in 'df_beta'
  ind.chr <- which(df_beta$chr == chr)
  ## indices in 'G'
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  
  corr0 <- snp_cor(Gimp, ind.col = ind.chr2, size = 3 / 1000,
                   infos.pos = POS2[ind.chr2], ncores = NCORES)

  
  if (chr == 1) {
    ld <- Matrix::colSums(corr0^2)
    corr <- as_SFBM(corr0, tmp, compact = TRUE)
  } else {
    ld <- c(ld, Matrix::colSums(corr0^2))
    corr$add_columns(corr0, nrow(corr))
  }
}

file.size(corr$sbk) / 1024^3

# Observed heritability
(ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                sample_size = n_eff, blocks = NULL)))
#       int        h2 
#  1.1091209 0.2878494 
h2_est <- ldsc[["h2"]]
#paste("The estimated heritability is", print(h2_est), sep = "\t")
h2_est

## LDPred2-inf
beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = h2_est)
#head(beta_inf)
pred_inf <- big_prodVec(Gimp, beta_inf, ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]]) 
#head(pred_inf)
cor(pred_inf, y[ind.test])
AUCBoot(pred_inf, y[ind.test])
#paste("AUC-inf is", AUCBoot(pred_inf, y[ind.test]), sep = "\t")

## LDpred2-auto
multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                               vec_p_init = seq_log(1e-4, 0.5, length.out = 30),
                               ncores = NCORES)
str(multi_auto, max.level = 1)
str(multi_auto[[1]], max.level = 1)

# Plot the convergence chains
auto <- multi_auto[[1]]
jpeg("/zivlab/data3/pmiddha/Proj_GeRI_GI-irAEs/PRS_UKBB/Results/CD/cd_ldpred2_chain_convergence_fplot.jpeg")
plot_grid(
  qplot(y = auto$path_p_est) +
    theme_bigstatsr() +
    geom_hline(yintercept = auto$p_est, col = "blue") +
    scale_y_log10() +
    labs(y = "p"),
  qplot(y = auto$path_h2_est) +
    theme_bigstatsr() +
    geom_hline(yintercept = auto$h2_est, col = "blue") +
    labs(y = "h2"),
  ncol = 1, align = "hv"
)
dev.off()

# Beta estimation and prediction in test set
beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
#head(beta_auto)
pred_auto <- big_prodMat(Gimp, beta_auto, ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]])
#head(pred_auto)

sc <- apply(pred_auto, 2, sd)
keep <- abs(sc - median(sc)) < 3 * mad(sc)
final_beta_auto <- rowMeans(beta_auto[, keep])
#head(final_beta_auto)
final_pred_auto <- rowMeans(pred_auto[, keep])
#head(final_pred_auto)

cor(final_pred_auto, y[ind.test])
AUCBoot(final_pred_auto, y[ind.test])

# Make sure it is on some appropriate scaling
# (the same as LDpred2-inf or a bit larger if the prediction is better)
# to safely merge predictions from multiple chromosomes
c(mad(pred_inf), mad(final_pred_auto))

# Combine SNP with beta_inf and final_beta_auto
beta_inf_snp <- cbind(df_beta, beta_inf)
final_beta_auto_snp <- cbind(df_beta, final_beta_auto)

# Store betas
write.csv(beta_inf_snp, file = "/zivlab/data3/pmiddha/Proj_GeRI_GI-irAEs/PRS_UKBB/Results/CD/cd_ldpred2_inf_model.csv", quote = F)
write.csv(final_beta_auto_snp, file = "/zivlab/data3/pmiddha/Proj_GeRI_GI-irAEs/PRS_UKBB/Results/CD/cd_ldpred2_auto_model.csv", quote = F)

# AUC plot
jpeg("/zivlab/data3/pmiddha/Proj_GeRI_GI-irAEs/PRS_UKBB/Results/CD/cd_ldpred2_auto_auc.jpeg")
rocobj <- plot.roc(y[ind.test], final_pred_auto,
                   main = "AUC-CD", 
                   percent=FALSE,
                   ci = TRUE,                  # compute AUC (of AUC by default)
                   print.auc = TRUE) 
dev.off()

# Predicted probability plot for cases and controls 
jpeg("/zivlab/data3/pmiddha/Proj_GeRI_GI-irAEs/PRS_UKBB/Results/CD/cd_ldpred2_auto_predicted_caco.jpeg")
ggplot(data.frame(
  Phenotype = factor(y[ind.test], levels = 0:1, labels = c("Control", "Case")),
  Probability = 1 / (1 + exp(-final_pred_auto)))) + 
  theme_bigstatsr() +
  geom_density(aes(Probability, fill = Phenotype), alpha = 0.3)
dev.off()

# Florian's paper section trial
#betas <- cbind(beta_inf, beta_auto = final_beta_auto)
#pred_test <- big_prodMat(Gimp, betas, ind.row = ind.test, ind.col = df_beta[["_NUM_ID_"]])
#resflo <- list(pred = setNames(as.data.frame(pred_test), colnames(betas)), auto = multi_auto[keep])
#saveRDS(resflo, file = "/zivlab/data3/pmiddha/Results/ukbb/gi_prs/ldpred2/florian_paper_style")

rm(corr)
gc()
file.remove(paste0(tmp, ".sbk"))
q("no")

