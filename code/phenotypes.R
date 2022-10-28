
library(bigsnpr)
library(tidyverse)

# Define the European-ancestry unrelated sample

# ipsych <- snp_attach("~/iPSYCH2015/HRC_Imputed/bigsnp_r_format/ipsych_hm3.rds")
# PC <- readRDS("~/iPSYCH2015/HRC_Imputed/bigsnp_r_format/PC.rds")
# PC <- readRDS("../gwas_VDBP/data/PCA_homogeneous45.rds")
# hist(log_dist <- log(bigutilsr::dist_ogk(PC)))
# is_homogeneous <- log_dist < 4.5
# 
# rel <- readRDS("~/iPSYCH2015/HRC_Imputed/bigsnp_r_format/rel.rds")
# is_rel <- ipsych$fam$sample.ID %in% subset(rel, KINSHIP > 2^-3.5)$IID2
# mean(is_rel) # 10.98%
# 
# mean(KEEP <- !is_rel & is_homogeneous) # 80.2%
# ind_keep <- which(KEEP)
# 
# saveRDS(tibble(pid = ipsych$fam$family.ID, is_homogeneous = is_homogeneous, is_rel = is_rel), "data/ipsych2015_pid.rds")


ipsych_hm3 <- readRDS("data/ipsych2015_pid.rds")

phen_file <- read.csv("/project/Register/2019_06/csv/ipsych2015design_v2.csv") %>%
  left_join(read_csv("/project/Register/2021_01/phenotype2016j.csv"), by = c("pid", "fdato")) %>%
  mutate(across(c(30:146), ~ ifelse(!is.na(.x), 1, 0))) %>%
  slice(match(ipsych_hm3$pid, pid))

idx <- match(ipsych_hm3$pid, phen_file$pid)

phen <- phen_file[idx,] %>% select(1, 21, 13:20, 30:146) %>%
  mutate(across(3:127, ~case_when(.x == 1 ~ 1,
                              kontrol2015I == 1 ~ 0)))

#saveRDS(phen, "data/phenotypes_ipsych2015_2016j.rds")

table_phen <- phen %>% 
  pivot_longer(2:127) %>%
  group_by(name, value) %>%
  tally() %>%
  pivot_wider(names_from = value, names_prefix = "N_", values_from = n)

#saveRDS(table_phen, "data/counts_ipsych2015_2016j.rds")


# Phen names
library(tabulizer)
pdf <- tabulizer::extract_tables("~/Register/2021_01/phenotype2016j.pdf", pages = 2:5, method = "stream")

tables <- do.call(rbind, lapply(pdf[-5], function(x) x[4:nrow(x),])) %>% as_tibble()

col_names <- c("Diagnosis", "ICD_10", "ICD_8", "Earliest_diagnosis", "Origin", "Variable")
  
# Overwrite names and remove Row 1
table_f <- tables %>%
  set_names(col_names) %>%
  mutate(Variable = case_when(Variable == "" ~ NA_character_,
                              TRUE ~ Variable)) %>%
  fill(Variable) %>%
  group_by(Variable) %>%
  summarise(diagnosis = toString(Diagnosis),
            icd10 = toString(ICD_10),
            icd8 = toString(ICD_8),
            origin = toString(Origin)) %>%
  mutate(origin = gsub(",", "", origin, fixed = TRUE),
         diagnosis = gsub(',', "", diagnosis, fixed = TRUE),
         diagnosis = gsub('- ', "", diagnosis, fixed = TRUE),
         diagnosis = gsub('-', "", diagnosis, fixed = TRUE)) %>%
  filter(Variable != "Variable")

#saveRDS(table_f, "data/phenotype_names_ipsych2015_2016j.rds")

# Covariates
PC <- readRDS("../gwas_VDBP/data/PCA_homogeneous45.rds")
colnames(PC) <- paste0("PC", 1:20)

cov <- bind_cols(tibble(birth_year = lubridate::year(lubridate::dmy(phen_file[idx,]$fdato)),
              sex = phen_file[idx,]$gender), is_2012 = ipsych$fam$is_2012, as_tibble(PC)) %>% 
  bigstatsr::covar_from_df(.) %>% as_tibble

#saveRDS(cov, "data/covariates_ipsych2015_2016j.rds")


#  5-fold Cross-validation

pars <- tibble(p = colnames(phen)[2:127]) %>%
  mutate(idx = map(p, ~which(!is.na(phen[[.x]]) & ipsych_hm3$is_homogeneous & !ipsych_hm3$is_rel)),
         cv_sub = map(idx, ~sample(1:5, length(.x), replace = TRUE))) %>%
  expand_grid(cv = 1:5) %>%
  mutate(ind.train = pmap(list(idx, cv_sub, cv), function(a,b,c) a[b != c]),
         ind.test = pmap(list(idx, cv_sub, cv), function(a,b,c) a[b == c])) %>%
  select(-idx)

#saveRDS(pars, "data/cv_subsets_ipsych2015_2016j.rds")

phen <- readRDS("data/phenotypes_ipsych2015_2016j.rds")
pars <- readRDS("data/cv_subsets_ipsych2015_2016j.rds")

counts_j <- pars %>%
  mutate(n_cases_train = map2_dbl(p, ind.train, ~sum(phen[[.x]][.y] == 1)),
         n_control_train = map2_dbl(p, ind.train, ~sum(phen[[.x]][.y] == 0)),
         n_cases_test = map2_dbl(p, ind.test, ~sum(phen[[.x]][.y] == 1)),
         n_control_test = map2_dbl(p, ind.test, ~sum(phen[[.x]][.y] == 0))) %>%
  select(p, cv, contains("n_"))

#saveRDS(counts_j, "results/counts_ipsych2015_2016j.rds")
