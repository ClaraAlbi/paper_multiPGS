library(tidyverse)
library(future.batchtools)
library(bigsnpr)

# Supplementary functions
se_from_or <- function(effect, se, maf) {
  varSNP <- 2*maf*(1 - maf)  
  beta <- effect/(effect^2 * varSNP + (pi^2)/3)^.5
  se.beta <- (se/exp(effect))/((effect^2 * varSNP + (pi^2)/3)^.5)
  m <- median((beta/se.beta)^2)
  t <- m > qchisq(0.05, 2) & m < qchisq(0.95, 2)
  return(list(t, beta, se.beta))
}
is_chi2 <- function(effect, se) {
  m <- median((effect/se)^2)
  t <- m > qchisq(0.05, 2) & m < qchisq(0.95, 2)
  return(t)
}

# Load list of files to process (from script 1)
sumstats_l <- readRDS("data/ncrr_gwascatalog_traits.rds") 

files_sub <- sumstats_l %>%
  filter(!is.na(beta_se) & (!is.na(or) | !is.na(beta))) %>%
  rename(trait = `Reported trait`)

NCORES <- 4
plan(batchtools_slurm(resources = list(
  t = "0-02:00", c = NCORES, mem = "250g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_pmap(files_sub[, c("trait", "ext", "basename", "gwas_file_o", "gwas_file_f", "gwas_file_qc", "qc_plot", "N", "Ncase", "Ncontrol", "chr", "pos", 'rsid', "a1", "a0", "beta", "or", "beta_se", "p", "Nsnp", "Nca", "Nco", "af")], 
                   function(trait, ext, basename, gwas_file_o, gwas_file_f, gwas_file_qc, qc_plot, N, Ncase, Ncontrol, chr, pos, rsid, a1, a0, beta, or, beta_se, p, Nsnp, Nca, Nco, af) {

  print(trait)
  
  print(c(chr, pos, rsid, a1, a0, beta, or, beta_se, p, Nsnp, Nca, Nco, af))
  
  map_maf <- data.table::fread("/home/clara/REPOS/BOLT-LMM_v2.3.4/tables/LDSCORE.1000G_EUR.tab.gz")
  
  # load sumstats
  if (str_detect(basename, "25673412")) {ss <- as_tibble(data.table::fread(gwas_file_o, skip = 8))
  } else if (any(str_detect(basename, as.character(288555:288664)))) {
    ss <- readr::read_table2(gwas_file_o)
  } else {
    if (ext == "zip") ss <- as_tibble(data.table::fread(unzip(gwas_file_o), fill = TRUE))
    if (ext != "zip") ss <- as_tibble(data.table::fread(gwas_file_o)) 
  } 
  print(ss)

  
  # if beta or or is.na
  if (!is.na(beta) | !is.na(or)) {
    if (!is.na(beta)) {if (sum(is.na(ss[,beta])) > 0.3*nrow(ss)) beta <- NA} 
    if (!is.na(or)) {if (sum(is.na(ss[,or])) > 0.1*nrow(ss)) or <- NA} 
    if (!is.na(af)) {if (sum(is.na(ss[,af])) > 0.50*nrow(ss)) af <- NA} 
    
    # Only contintue if there are variants with beta/or
    if (!is.na(beta) | !is.na(or)) {

      # fix if chr or/and pos empty
      if (is.na(chr) & !is.na(rsid)) {
        m <- sum(map_maf[["SNP"]] %in% ss[[rsid]])
        if (m > 0) {
          var_rsid <- colnames(ss)[rsid]
          join_cols = c("SNP")
          names(join_cols) <- var_rsid
          ss <- ss %>% left_join(map_maf[,1:3], by = join_cols)
          chr <- which(colnames(ss) == "CHR")
          pos <- which(colnames(ss) == "BP")
          #rsid <- which(colnames(ss) == join_cols)
        }
        if (m == 0) {
          ss[["chr"]] <- as.numeric(str_extract(rsid, "[^:]+"))
          ss[["pos"]] <- as.numeric(str_match(rsid, ":(.*?)_")[,2])
          chr <- which(colnames(ss) == "chr")
          pos <- which(colnames(ss) == "pos")
          stopifnot(sum(is.na(ss[["chr"]])) < 0.2*nrow(ss))
        }
      }
      # Fix when rsid is the original name
      #if (!is.na(rsid) & colnames(ss)[rsid] == "rsid") ss <- ss %>% rename(RSID_o = rsid)
      
      print(ss)
      print(c(chr, rsid, pos))
      # rename location columns
      if (isTRUE(all(chr, pos, rsid))) {
        print("All location columns present")
        ss_f <- ss %>% rename(chr = all_of(chr), pos = all_of(pos), rsid = all_of(rsid), a1 = all_of(a1), a0 = all_of(a0),  p = all_of(p)) %>%
          mutate(a1 = toupper(a1), a0 = toupper(a0)) %>%
          mutate(chr = as.numeric(str_extract(chr, "[^chr]+")),
                 rsid = str_extract(rsid, "[^:]+"))
      }
      if (isTRUE(all(chr, pos)) & is.na(rsid)) {
        print("Missing rsid")
        ss_f <- ss %>% rename(chr = all_of(chr), pos = all_of(pos),  a1 = all_of(a1), a0 = all_of(a0),  p = all_of(p)) %>%
          mutate(a1 = toupper(a1), a0 = toupper(a0)) %>%
          mutate(chr = as.numeric(str_extract(chr, "[^chr]+")))
      }
      if (isTRUE(all(chr, rsid)) & is.na(pos)) {
        print("Missing pos")
        ss_f <- ss %>% rename(chr = all_of(chr), rsid = all_of(rsid),  a1 = all_of(a1), a0 = all_of(a0),  p = all_of(p)) %>%
          mutate(a1 = toupper(a1), a0 = toupper(a0)) %>%
          mutate(chr = as.numeric(str_extract(chr, "[^chr]+")),
                 rsid = str_extract(rsid, "[^:]+"))
      }
      # If is.na(af), replace with 1KG AF
      if (is.na(af)) {
        print("Setting AF to 1KG")
        if (isTRUE(all(chr, rsid))) {
          ss_f <- ss_f %>% inner_join(map_maf[,c(1:2, 4)], by = c("rsid" = "SNP")) %>% rename(af = MAF)
          af <- which(colnames(ss_f) == "af")
        }
        
        if (isTRUE(all(chr, pos)) & is.na(rsid)) {
          ss_f <- ss_f %>% inner_join(map_maf[,2:4], by = c("chr" = "CHR", "pos" = "BP")) %>% rename(af = MAF)
          af <- which(colnames(ss_f) == "af")
        }
      } else {ss_f <- ss_f %>% rename(af = all_of(af)) %>% mutate(af = as.numeric(af))}

      # which n_eff
      if (!is.na(Nca) & !is.na(Nco)) ss_f$n_eff <- 4/(1/ss_f[[Nca]] + 1/ss_f[[Nco]])
      if (is.na(Nca) & !is.na(Nsnp)) ss_f$n_eff <- ss_f[[Nsnp]]
      if (is.na(Nca) & is.na(Nsnp)) {
        if (is.na(Ncase)) ss_f$n_eff <- N
        if (!is.na(Ncase)) ss_f$n_eff <- 4/(1/Ncase + 1/Ncontrol)
      }

      # GET EFFECTS AND SE
      
      if (beta == beta_se) {
        ss_f <- ss_f %>% filter(n_eff != -9) %>%
          rename(beta = all_of(beta)) %>%
          mutate(beta_se = 1) %>%
          drop_na(beta, beta_se, n_eff, af) %>%
          mutate_at(vars(beta, beta_se), as.numeric)
        print("Score without se")
        se_from <- "zscore"
      }

      # case when beta and se
      if (!is.na(beta) & !is.na(beta_se) & is.na(or)) {
        ss_f <- ss_f %>% filter(n_eff != -9) %>%
          rename(beta = all_of(beta), beta_se = all_of(beta_se)) %>%
          drop_na(beta, beta_se, n_eff, af) %>%
          filter(beta_se != 0) %>%
          mutate_at(vars(beta, beta_se), as.numeric)
        print("No transformation of beta/se")
        se_from <- "beta"

      }
      # case when or and se
      if (!is.na(or) & !is.na(beta_se) & is.na(beta)) {
        ss_f <- ss_f %>% filter(n_eff != -9) %>%
          rename(beta_se = all_of(beta_se), or = all_of(or)) %>%
          mutate(beta = log(as.numeric(or))) %>%
          drop_na(or, beta, beta_se, n_eff, af) %>%
          filter(beta_se != 0) %>%
          mutate_at(vars(beta, or, beta_se), as.numeric)

        v_or <- with(ss_f, se_from_or(or, beta_se, af))

        if (with(ss_f, is_chi2(beta, beta_se))) {
          print("No transformation")
          se_from <- "log_or"

        } else if (v_or[[1]]) {
          ss_f[["beta"]] <- v_or[[2]]
          ss_f[["beta_se"]] <- v_or[[3]]
          print("SE from OR, Transform effect and se from se_or")
          se_from <- "or"
          }
      }

      # case when beta and or and se
      if (!is.na(beta) & !is.na(beta_se) & !is.na(or)) {
        ss_f <- ss_f %>% filter(n_eff != -9) %>%
          rename(beta = all_of(beta), beta_se = all_of(beta_se), or = all_of(or)) %>%
          mutate(log_or = log(as.numeric(or))) %>%
          drop_na(beta, beta_se, log_or, n_eff, or, af) %>%
          filter(beta_se != 0) %>%
          mutate_at(vars(beta, or, beta_se), as.numeric)

        v_or <- with(ss_f, se_from_or(or, beta_se, af))

        if (with(ss_f, is_chi2(log_or, beta_se))) {
          ss_f$beta <- ss_f$log_or
          print("No transformation, se from log_or")
          se_from <- "log_or"
        } else if (with(ss_f, is_chi2(beta, beta_se))) {
          print("No transformation, se from beta")
          se_from <- "beta"

        } else if (v_or[[1]]) {
          ss_f[["beta"]] <- v_or[[2]]
          ss_f[["beta_se"]] <- v_or[[3]]
          print("SE from OR, Transform effect and se from se_or")
          se_from <- "or"
        }
      }


      chi2 <- with(ss_f, median((beta/beta_se)^2))
      stopifnot(with(ss_f, is_chi2(beta, beta_se)))
      print(paste0("chi2: ", chi2))
      print(head(ss_f))
      
      
      # Run LDSC
      
      map_ldsc <- map_maf %>% rename(rsid = SNP, chr = CHR, pos = BP) %>% select(-MAF)
      
      if (isTRUE(all(chr, rsid))) {
        ss_ld <- ss_f %>% inner_join(map_ldsc, by = c("rsid")) %>% filter(!is.na(LDSCORE)) 
        
      }
      
      if (isTRUE(all(chr, pos)) & is.na(rsid)) {
        ss_ld <- ss_f %>% inner_join(map_ldsc, by = c("chr", "pos")) %>% filter(!is.na(LDSCORE)) %>% drop_na(beta, beta_se, n_eff)
       
      }
      
      print("Calculating h2 with LDSC")
      
      
      (ldsc <- with(ss_ld, snp_ldsc(LDSCORE, nrow(map_maf), chi2 = (beta / beta_se)^2,
                                       sample_size = n_eff, ncores = 16)))
      print(ldsc)
      M_or <- nrow(ss_f)
      
    } else {
      print("No not NA variants")
      se_from <- NA
      ss_f <- NA
      M_or <- NA
      chi2 <- NA
      ldsc <- rep(NA, 4)
      names(ldsc) <- c("int", "int_se", "h2", "h2_se")
    }
  } 
  
  saveRDS(c(list(sumstats = ss_f, M_or = M_or, chi2 = chi2, se = se_from), ldsc), gwas_file_f)
})


library(bigparallelr)
registerDoParallel(cl <- makeCluster(24))

all_res_files <- list.files("data/sumstats_f", ".rds", full.names = T)

res <- foreach(res_file = all_res_files) %dopar% readRDS(res_file)
res_l <- do.call(bind_rows, res)
stopCluster(cl)

qc_res <- cbind(all_res_files, res_l) %>%
  mutate(basename = str_match(all_res_files, "data/sumstats_f/*(.*?).rds")[,2]) %>%
  select(-all_res_files)

#saveRDS(qc_res, "data/sumstats_f.rds")


