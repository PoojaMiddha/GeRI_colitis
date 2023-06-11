# Prepare summary stats from deLange paper as per LDPred2 tutorial (UC)
# Basically the summary stats are matched/restricted to HapMap3 variants 
# Author: Pooja Middha
# Date: Apr 20, 2022

# Load library
library(bigsnpr)
library(bigreadr)
library(data.table)
library(dplyr)
library(ggplot2)

# Number of cores
NCORES <- nb_cores() # or mention 30 or some number based on core capacity

# Information for the variants provided in the LD reference
map <- readRDS("/zivlab/data3/pmiddha/Proj_GeRI_GI-irAEs/PRS_UKBB/Data/genetic/map_hm3_ldpred2.rds") 

# Load summary statistics
load("/zivlab/data3/pmiddha/Data/sumstats_delange_mrcieu/uc_summary_stats.Rdata")
sumstats <- df2 %>% select(rsid, chr, pos, ref, alt, beta, se, p)
colnames(sumstats)
colnames(sumstats)[4] <- "a0"
colnames(sumstats)[5] <- "a1"
colnames(sumstats)[7] <- "beta_se"
colnames(sumstats)

# Adding number of cases and controls from the paper
sumstats$n_cases <- as.numeric(12366)
sumstats$n_controls <- as.numeric(33609)
# The numbers are from Supplementary table 1 of DeLange paper (2017)

# Calculate n_eff as per LDPred2 tutorial
sumstats$n_eff <- 4/(1 / sumstats$n_cases + 1 / sumstats$n_controls)

# Match reference data from LDPred2 and summary stats
info_snp <- tibble::as_tibble(snp_match(sumstats, map))
#9,459,549 variants to be matched.
#0 ambiguous SNPs have been removed.
#1,017,787 variants have been matched; 0 were flipped and 0 were reversed.

# Standard deviation for variants in reference and summary stats
sd_ldref <- with(info_snp, sqrt(2 * af_UKBB * (1 - af_UKBB)))
sd_ss <- with(info_snp, 2 / sqrt(n_eff * beta_se^2)) 

# Assess how many variants are not fit to take forward based on standard deviation calculations
is_bad <-
  sd_ss < (0.5 * sd_ldref) | sd_ss > (sd_ldref + 0.1) | sd_ss < 0.1 | sd_ldref < 0.05

# Plot the above standard deviations and show bad variants
jpeg("/zivlab/data3/pmiddha/Data/sumstats_delange_mrcieu/ldpred2_prepare_uc_sumstats_sd.jpeg")
qplot(sd_ldref, sd_ss, color = is_bad, alpha = I(0.5)) +
  theme_bigstatsr() +
  coord_equal() +
  scale_color_viridis_d(direction = -1) +
  geom_abline(linetype = 2, color = "red") +
  labs(x = "Standard deviations derived from allele frequencies of the LD reference",
       y = "Standard deviations derived from the summary statistics",
       color = "To remove?")
dev.off()

# Remove bad variants from summary stats
df_beta <- info_snp[!is_bad, ]
dim(df_beta)

sumstats_new <- df_beta %>% select(rsid, chr, pos, a0, a1, beta, beta_se, p, n_cases, n_controls, n_eff)

saveRDS(sumstats_new, file = "/zivlab/data3/pmiddha/Data/sumstats_delange_mrcieu/ldpred2_uc_sumstats_hapmapmatched.rds") 

q("no")







