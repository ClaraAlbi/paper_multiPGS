library(tidyverse)
library(bigsnpr)
library(glue)

G <- snp_attach("~/iPSYCH2015/HRC_Imputed/bigsnp_r_format/ipsych_hm3.rds")$genotypes
ldref <- readRDS("data/LD_with_blocks_UKB/map.rds")

files_sub <- readRDS("data/ncrr_gwascatalog_traits.rds") %>%
  filter(file.exists(res_file) & file.exists(gwas_file_qc)) %>%
  filter(!file.exists(glue("results/prs/{basename}.rds")))


library(future.batchtools)
NCORES <- 8
plan(batchtools_slurm(resources = list(
  t = "0-04:00", c = NCORES, mem = "250g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))


furrr::future_pmap(files_sub[, c("basename", "res_file", "gwas_file_qc")], function(basename, res_file, gwas_file_qc) {
  
  sum_file <- glue("results/prs/{basename}.rds")
  
  sumstats <- ldref %>% select(rsid, chr) %>%
    left_join(readRDS(gwas_file_qc)[["sumstats"]]) %>%
    mutate(idx = row_number()) %>%
    filter(!is.na(beta))
  
  # get M
  M_ldpred <- nrow(sumstats)
  
  #Filter chains
  multi_auto <- readRDS(res_file)[["auto"]]
  ldsc <- readRDS(res_file)[c("int", "int_se", "h2", "h2_se")]
  
  range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
  keep <- (range > (0.9 * quantile(range, 0.9)))
  k <- sum(keep)
  
  # Chenge h2 and p
  h2_est <- sapply(multi_auto, function(auto) auto$h2_est)[keep]
  p_est <- sapply(multi_auto, function(auto) auto$p_est)[keep]
  
  # Change pred
  beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
  
  beta_auto_m <- rowMeans(beta_auto[, keep])
  pred_auto <- big_prodVec(G, beta_auto_m, ind.col = sumstats[["_NUM_ID_"]])
  
  res <- c(list(M_ldpred = M_ldpred,
       k = k,
       h2_est = h2_est,
       p_est = p_est,
       pred_auto = pred_auto), ldsc)
  
  saveRDS(res, sum_file)
  
})

library(bigparallelr)
registerDoParallel(cl <- makeCluster(6))

f <- list.files("results/prs", full.names = T)
res <- foreach(res_file = f) %dopar% readRDS(res_file)
res_l <- do.call(rbind, res)

# Parameters

pars <- as_tibble(res_l) %>% mutate(basename = f,
                                    basename = str_match(basename, "results/prs/*(.*?).rds")[,2]) %>%
  select(-pred_auto) %>%
  mutate(across(c(h2_est, p_est), ~map_dbl(.x, median))) %>%
  unnest(c(int, int_se, h2, h2_se, k, M_ldpred))

#saveRDS(pars, "results/all_pars_v2.rds")

# PRS

pars <- readRDS("results/all_pars_v2.rds")

pars %>%
  filter(M_ldpred > 5e5) %>%
  ggplot(aes(x = h2_est, y = h2, color = p_est)) + 
  geom_errorbar(aes(ymin = h2 - h2_se*1.96, ymax = h2 + h2_se*1.96),
                color = "grey") +
  geom_point(alpha = 0.5) +
  xlab("h2_ldpred2") + 
  ylab("h2_ldsc") + 
  theme_minimal() +
  scale_color_viridis_c() + 
  geom_abline(slope = 1, linetype = 2)

pars %>%
  filter(M_ldpred > 200000) %>%
  filter(h2_est < 1 & h2 < 1) %>%
  filter(p_est < 0.1) %>%
  ggplot(aes(x = h2_est, y = h2, color = log10(p_est))) + 
  geom_errorbar(aes(ymin = h2 - h2_se*1.96, ymax = h2 + h2_se*1.96),
                color = "grey") +
  geom_point(alpha = 0.5) +
  xlab("h2_ldpred2") + 
  ylab("h2_ldsc") + 
  theme_minimal() +
  scale_color_viridis_c() + 
  geom_abline(slope = 1, linetype = 2)
  

prs <- as_tibble(res_l) %>% mutate(basename = f,
                                   basename = str_match(basename, "results/prs/*(.*?)_")[,2]) %>%
  select(basename, pred_auto) %>%
  pivot_wider(names_from = basename, values_from = pred_auto) %>%
  unnest()

#saveRDS(prs, "results/all_prs_v2.rds")

prs_std <- prs %>% mutate(across(everything(), ~scale(.x)[,1]))

#saveRDS(prs_std, "../multi_prs/data/prs_all_std_v2.rds")


