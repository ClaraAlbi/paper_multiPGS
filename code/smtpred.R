library(tidyverse)
library(bigsnpr)

sumstats <- readRDS("../gwas_prs/results/all_sumstats.rds")
phen <- readRDS("data/phenotypes_ipsych2015_2016j.rds")

res_lasso <- readRDS("results/weights_lasso_ipsych.rds") %>%
  filter(!n %in% c("sexM", "birth_year")) %>%
  filter(abs(w_m) > 0.01) %>%
  mutate(N_p1 = map_dbl(p, ~4/(1/sum(phen[[.x]] == 1, na.rm = T) + 1/sum(phen[[.x]] == 0, na.rm = T))),
         N_p2 = map_dbl(n, ~case_when(is.na(sumstats$Ncase[sumstats$id == .x]) ~ sumstats$N[sumstats$id == .x],
                                      !is.na(sumstats$Ncase[sumstats$id == .x]) ~ 4/(1/sumstats$Ncase[sumstats$id == .x] + 1/sumstats$Ncontrol[sumstats$id == .x])))) %>%
  filter(file.exists(paste0("results/ldsc/", p, "_", n, ".log"))) %>%
  mutate(filename = paste0("results/ldsc/", p, "_", n, ".log"),
         pars = map(filename, ~data.table::fread(.x, skip = "p1", fill = T, nrows = 1))) %>%
  rename(phen = p) %>%
  unnest(pars) %>%
  group_by(phen) %>%
  summarise(n = list(n),
            h2_obs = list(h2_obs),
            rg = list(rg),
            neff = list(N_p2)) %>%
  left_join(n, c("phen" = "p"))

ipsych <- snp_attach("~/iPSYCH2015/HRC_Imputed/bigsnp_r_format/ipsych_hm3.rds")

pmap(res_lasso[6,c(1:2, 6)], function(phen, n, id) {
  
  
  # Get all GWAS sumstats for each PRS and each phen
  sumstats <- list.files("../gwas_prs/data/sumstats_qc", pattern = paste0(paste0(n, "_"), collapse="|"), full.names = T)
  ss_path <- c()
  for (file in n) {
    tmp_1 <- paste0("tmp-data/", file)
    ss_path <- c(ss_path, tmp_1)
    if (!file.exists(tmp_1)) {
      sumstats_p2 <- readRDS(sumstats[str_detect(sumstats, paste0(file, "_"))])$sumstats %>%
        select(chr, pos, rsid, a1, a0, beta, p, n_eff) %>%
        rename(CHR = chr, BP = pos, SNP = rsid, A1 = a1, A2 = a0, N = n_eff)
      
      write_delim(sumstats_p2, tmp_1)
    }}
  
  
  system(glue::glue("python ~/REPOS/smtpred/ldsc_wrapper.py ",
                    " --ldscpath ~/REPOS/ldsc/",
                    " --files {paste0(ss_path, collapse = ' ')}",
                    " --snplist /project/NCRR-PRS/faststorage/clara/gwas_VDBP/data/eur_w_ld_chr/w_hm3.snplist",
                    " --ref_ld /project/NCRR-PRS/faststorage/clara/gwas_VDBP/data/eur_w_ld_chr/",
                    " --w_ld /project/NCRR-PRS/faststorage/clara/gwas_VDBP/data/eur_w_ld_chr/",
                    " --out results/smtpred/"))
  
  # GET all PRSs for each phen
  dir.create(paste0("tmp-data/prs/", phen, "/"), showWarnings = F)
  prs_paths <- c()
  for (s in n) {
    tmp_2 <- paste0("tmp-data/prs/", phen, "/", s)
    prs_paths <- c(prs_paths, tmp_2)
    if (!file.exists(tmp_2)) {
      prs_p2 <- readRDS(list.files("../gwas_prs/results/prs/", pattern = paste0(regex(s), "_"), full.names = T))$pred_auto
      
      write_delim(tibble(FID = ipsych$fam$family.ID, IID = ipsych$fam$sample.ID, PHENO = NA, CNT = NA, CNT2 = NA, SCORE =  scale(prs_p2)[,1]), tmp_2)
      }
  } 
  id_file <- prs_paths[str_detect(prs_paths, id)]
  prs_paths <- c(id_file, prs_paths[-which(prs_paths == id_file)])


  system(glue::glue("python ~/REPOS/smtpred/smtpred.py ",
                    " --h2file results/smtpred/ldsc_h2s.txt",
                    " --rgfile results/smtpred/ldsc_rgs.txt",
                    " --nfile results/smtpred/ldsc_ns.txt",
                    " --scorefiles {paste0(prs_paths, collapse = ' ')}",
                    #" --scorepath tmp-data/prs/{phen}/",
                    " --out results/smtpred/{phen}"))
  on.exit(file.remove(prs_paths), add = TRUE)
})




