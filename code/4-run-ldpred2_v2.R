
library(tidyverse)
library(bigsnpr)
library(purrr)

ipsych_test <- snp_attach("~/iPSYCH2015/HRC_Imputed/bigsnp_r_format/ipsych_hm3.rds")
G <- ipsych_test$genotypes

ldref <- readRDS("data/LD_with_blocks_UKB/map.rds")

name_corr <- "data/LD_with_blocks_UKB/LD_with_blocks_chr"

bigassertr::assert_dir("results/ldpred2_v2")
bigassertr::assert_dir("plots/conv_v2")

files_sub <- readRDS("data/ncrr_gwascatalog_traits.rds") %>%
  filter(!file.exists(res_file) & file.exists(gwas_file_qc)) %>%
  left_join(readRDS("data/sumstats_qc.rds")) %>%
  filter(M_qc > 200000)


library(future.batchtools)
NCORES <- 32
plan(batchtools_slurm(resources = list(
  t = "0-24:00", c = NCORES, mem = "250g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_pmap(files_sub[, c("basename", "res_file", "gwas_file_qc")], function(basename, res_file, gwas_file_qc) {
  
  
  sumstats <- ldref %>% select(rsid, chr) %>%
    left_join(readRDS(gwas_file_qc)[["sumstats"]]) %>%
    mutate(idx = row_number()) %>%
    filter(!is.na(beta))
  
  tmp <- tempfile(tmpdir = "tmp-data")
  on.exit(file.remove(paste0(tmp, ".sbk")), add = TRUE)
  
  for (chr in 1:22) {
    print(chr)
    ## indices in 'sumstats'
    ind.chr <- which(sumstats$chr == chr)
    ## indices in 'G'
    ind.chr2 <- sumstats$idx[ind.chr]
    ## indices in 'corr'
    ind.chr3 <- match(ind.chr2, which(ldref$chr == chr))
    
    corr0 <- readRDS(paste0(name_corr, chr, ".rds"))
    ld_chr <- Matrix::colSums(corr0^2)
    df_beta_chr <- sumstats[ind.chr, c("beta", "beta_se", "n_eff", "_NUM_ID_")]
    
    if (chr == 1) {
      df_beta <- df_beta_chr
      ld <- ld_chr[ind.chr]
      corr <- bigsparser::as_SFBM(corr0[ind.chr3, ind.chr3], tmp, compact = TRUE)
    } else {
      df_beta <- rbind(df_beta, df_beta_chr)
      ld <- c(ld, ld_chr[ind.chr3])
      corr$add_columns(corr0[ind.chr3, ind.chr3], nrow(corr))
    }
  }

  (ldsc <- with(df_beta, snp_ldsc(ld, length(ld), chi2 = (beta / beta_se)^2,
                                  sample_size = n_eff)))
  h2_est <- ldsc[["h2"]]
  
  # LDpred2-inf
  beta_inf <- snp_ldpred2_inf(corr, df_beta, h2_est)
  
  # LDpred2-auto
  multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = h2_est,
                                 vec_p_init = seq_log(1e-4, 0.9, 30),
                                 allow_jump_sign = FALSE, shrink_corr = 0.95,
                                 ncores = 30, burn_in = 800, num_iter = 400)

  beta_auto <- sapply(multi_auto, function(auto) auto$beta_est)
  pred_auto <- big_prodMat(G, beta_auto, ind.col = df_beta[["_NUM_ID_"]])
  
  sc <- apply(pred_auto, 2, sd)
  keep <- abs(sc - median(sc)) < 3 * mad(sc)
  final_beta_auto <- rowMeans(beta_auto[, keep])
  
  # compute prs
  all_beta <- cbind(beta_inf, final_beta_auto)
  pred <- big_prodMat(G, all_beta, ind.col = df_beta[["_NUM_ID_"]])
  
  res <- c(list(pred = setNames(as.data.frame(pred), colnames(all_beta)), auto = multi_auto, keep = keep, M = nrow(df_beta)), ldsc)
  
  saveRDS(res, res_file)
  
  auto <- multi_auto[[1]]
  plot2 <- plot_grid(
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
  ggsave(paste0("plots/conv_v2/", sub('\\.rds$', '', basename), ".png"), plot = plot2, width = 12, height = 7)
  
  
})


