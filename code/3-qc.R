library(tidyverse)
library(bigsnpr)
library(future.batchtools)
library(data.table)
source("code/match_fun.R")

# Map file & sd of G
ipsych15 <- bigsnpr::snp_attach("~/iPSYCH2015/HRC_Imputed/bigsnp_r_format/ipsych_hm3.rds")
G <- ipsych15$genotypes
map <- ipsych15$map %>%
  rename(chr = CHR, rsid = SNP, pos = POS, a0 = a2)
NCORES <- 16

ind.val <- readRDS("data/ind_row_ipsych.rds")
sd <- runonce::save_run(
  big_parallelize(G, function(X, ind, ind.val) {
    bigstatsr::big_scale()(X, ind.val, ind)$scale
  }, p.combine = "c", ncores = NCORES, ind.val = ind.val),
  file = "data/sd.rds"
)
# user  system elapsed 
# 0.397   0.995 734.342

sumstats_f <- readRDS("data/ncrr_gwascatalog_traits.rds") %>%
  inner_join(readRDS("data/sumstats_f.rds")) %>%
  filter(file.exists(gwas_file_f)) %>%
  rename(trait = `Reported trait`) 

NCORES <- 1
plan(batchtools_slurm(resources = list(
  t = "0-01:00", c = NCORES, mem = "16g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))


furrr::future_pmap(sumstats_f, function(trait, basename, gwas_file_f, gwas_file_qc, qc_plot, rsid, pos, se, chr) {

  ss_f <- readRDS(gwas_file_f)[[1]]
  stopifnot(length(ss_f) != 1)
  
  if (se == "zscore") {
    ss_f$beta_se <- with(ss_f, 1/sqrt(2*af*(1 - af)*(n_eff + beta^2)))
    ss_f$beta <- with(ss_f, beta * beta_se)
    print("beta and se inferred from zscore")
  }
  
  #MATCHING
  if (!is.na(chr)) {if (sum(is.na(ss_f[,chr])) > 0.3*nrow(ss_f)) chr <- NA} 
  if (!is.na(rsid)) {if (sum(is.na(ss_f[,rsid])) > 0.3*nrow(ss_f)) rsid <- NA} 
  if (!is.na(pos)) {if (sum(is.na(ss_f[,pos])) > 0.3*nrow(ss_f)) pos <- NA} 
  
  #stopifnot(!is.na(chr))
  
  if (!is.na(rsid)) {sumstats_m_rsid <- match_sumstats(sumstats = ss_f, info_snp = map[,c(1:3, 7:8)], join_by_pos = FALSE, match.min.prop = 0) 
  } 
  if (!is.na(pos)) {sumstats_m_pos <- match_sumstats(sumstats = ss_f, info_snp = map[, c(1:3, 7:8)], join_by_pos = TRUE, match.min.prop = 0) 
  }
  
  
  if (exists("sumstats_m_rsid") & exists("sumstats_m_pos")) {
    if (nrow(sumstats_m_rsid) > nrow(sumstats_m_pos)) {sumstats_m <- sumstats_m_rsid
    } else {sumstats_m <- sumstats_m_pos}
  } else if (exists("sumstats_m_rsid") & !exists("sumstats_m_pos")) {
    sumstats_m <- sumstats_m_rsid 
    print("Joining by rsid")
  } else if (!exists("sumstats_m_rsid") & exists("sumstats_m_pos")) {
    sumstats_m <- sumstats_m_pos
    print("Joining by pos")
  }
  print(paste("# Matched SNPs:", nrow(sumstats_m)))
  
  
  # Get slope between sd_G and se_val*sqrt(Neff_val)
  
  sd_val <- sd[sumstats_m$`_NUM_ID_`]
  if (se == "beta") {
  #if (TRUE) {
    slope <- median(sd_val * sumstats_m$beta_se * sqrt(sumstats_m$n_eff))
    print("Slope from lin reg")
  } else if (se %in% c("log_or", "or"))  {
    slope <- 2
    print("Slope from log reg")
  }
  
  print(slope)
  sd_ss <- with(sumstats_m, slope / (beta_se * sqrt(n_eff)))
  
  # Remove SNPs with large diff between sd_ss-sd_val
  is_bad <- 
    sd_ss < (0.5 * sd_val) | sd_ss > (sd_val + 0.1) | sd_ss < 0.1 | sd_val < 0.05
  
  print(paste("Filtered", sum(is_bad), "SNPs"))
  
  #Plot summary of QC
  plot_qc <- qplot(sd_val, color = is_bad, sd_ss, alpha = I(0.5)) +
    coord_fixed() +
    scale_color_viridis_d(direction = -1) +
    geom_abline(linetype = 2, color = "red") +
    labs(x = "Sd val",
         y = "Sd ss",
         color = "Removed?",
         title = trait,
         subtitle = paste0("# removed snps: ", sum(is_bad, na.rm = TRUE), " out of ", nrow(sumstats_m)),
         caption = paste0("slope: ", slope, "\n"))
  ggsave(qc_plot, plot = plot_qc, width = 10, height = 5)
  
  saveRDS(list(sumstats = sumstats_m[!is_bad, ], slope_qc = slope, M_or = nrow(ss_f), M_m = nrow(sumstats_m), M_qc = nrow(sumstats_m[!is_bad, ])), gwas_file_qc)
  
})



library(bigparallelr)
registerDoParallel(cl <- makeCluster(24))

all_res_files <- list.files("data/sumstats_qc", ".rds", full.names = T)

res <- foreach(res_file = all_res_files) %dopar% readRDS(res_file)[2:5]
res_l <- do.call(bind_rows, res)
stopCluster(cl)

qc_res <- cbind(all_res_files, res_l) %>%
  mutate(basename = str_match(all_res_files, "data/sumstats_qc/*(.*?).rds")[,2]) %>%
  select(-all_res_files)

#saveRDS(qc_res, "data/sumstats_qc.rds")


