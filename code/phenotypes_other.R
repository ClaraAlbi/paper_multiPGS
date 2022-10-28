library(tidyverse)

ipsych_hm3 <- readRDS("data/ipsych2015_pid.rds")

phen_file <- read.csv("/project/Register/2019_06/csv/ipsych2015design_v2.csv") %>%
  left_join(read_csv("~/Register/Supplements/MBR_2019_09/csv/iPSYCH2015_MBR.csv"), by = c("pid", "fdato")) %>%
  left_join(read_csv("~/Register/Supplements/Severity/severity2016_ipsych2015.csv"), by = c("pid")) %>%
  slice(match(ipsych_hm3$pid, pid))

idx <- match(ipsych_hm3$pid, phen_file$pid)

phen <- phen_file[idx,] #%>% select(1, apgar5, gest_age, fvagt, d2100_ptype0_contacts) %>%
  mutate(across(2:4, as.numeric)) 

phen$scz_hos <- qnorm((rank(phen$d2100_ptype0_contacts, na.last = "keep") - 0.5) / sum(!is.na(phen$d2100_ptype0_contacts)))

#saveRDS(phen, "data/phenotypes_ipsych2015_2016MBR.rds")

#  5-fold Cross-validation

pars <- tibble(p = colnames(phen)[2:6]) %>%
  mutate(idx = map(p, ~which(!is.na(phen[[.x]]) & ipsych_hm3$is_homogeneous & !ipsych_hm3$is_rel)),
         cv_sub = map(idx, ~sample(1:5, length(.x), replace = TRUE))) %>%
  expand_grid(cv = 1:5) %>%
  mutate(ind.train = pmap(list(idx, cv_sub, cv), function(a,b,c) a[b != c]),
         ind.test = pmap(list(idx, cv_sub, cv), function(a,b,c) a[b == c])) %>%
  select(-idx)

#saveRDS(pars, "data/cv_subsets_ipsych2015_2016MBR_v2.rds")

# Bipolar or depression?
sum(phen$bipol2015I & phen$affek2015I, na.rm = T)
sum(phen$bipol2015I, na.rm = T)
sum(phen$affek2015I, na.rm = T)

phen <- phen %>% 
  mutate(bip = case_when(bipol2015I == 1 ~ 1,
                         affek2015I == 1 & bipol2015I == 0 ~ 0))

#saveRDS(phen, "data/phenotypes_ipsych2015_2016bip.rds")

pars <- tibble(p = "bip") %>%
  mutate(idx = map(p, ~which(!is.na(phen[[.x]]) & ipsych_hm3$is_homogeneous & !ipsych_hm3$is_rel)),
         cv_sub = map(idx, ~sample(1:5, length(.x), replace = TRUE))) %>%
  expand_grid(cv = 1:5) %>%
  mutate(ind.train = pmap(list(idx, cv_sub, cv), function(a,b,c) a[b != c]),
         ind.test = pmap(list(idx, cv_sub, cv), function(a,b,c) a[b == c])) %>%
  select(-idx)

#saveRDS(pars, "data/cv_subsets_ipsych2015_2016bip.rds")

# ADHD OR ASD???

sum(phen$adhd2015I & phen$autism2015I, na.rm = T) #4939
sum(phen$adhd2015I, na.rm = T) #27551
sum(phen$autism2015I, na.rm = T) #23217

phen <- phen %>% 
  mutate(adhd_asd = case_when(adhd2015I == 1 & autism2015I == 0 ~ 1,
                         adhd2015I == 0 & autism2015I == 1 ~ 0))

pars <- tibble(p = "adhd_asd") %>%
  mutate(idx = map(p, ~which(!is.na(phen[[.x]]) & ipsych_hm3$is_homogeneous & !ipsych_hm3$is_rel)),
         cv_sub = map(idx, ~sample(1:5, length(.x), replace = TRUE))) %>%
  expand_grid(cv = 1:5) %>%
  mutate(ind.train = pmap(list(idx, cv_sub, cv), function(a,b,c) a[b != c]),
         ind.test = pmap(list(idx, cv_sub, cv), function(a,b,c) a[b == c])) %>%
  select(-idx)

#saveRDS(phen, "data/phenotypes_ipsych2015_2016adhd.rds")
#saveRDS(pars, "data/cv_subsets_ipsych2015_2016adhd.rds")

#COVARIATES
PC <- readRDS("../gwas_VDBP/data/PCA_homogeneous45.rds")
colnames(PC) <- paste0("PC", 1:20)

cov <- bind_cols(tibble(birth_year = lubridate::year(lubridate::dmy(phen_file[idx,]$fdato)),
                        sex = phen_file[idx,]$gender.x), is_2012 = ipsych$fam$is_2012, as_tibble(PC), 
                 is_cohort = phen_file[idx,]$kontrol2015I) %>% 
  bigstatsr::covar_from_df(.) %>% as_tibble

#saveRDS(cov, "data/covariates_ipsych2015_2016K.rds")
