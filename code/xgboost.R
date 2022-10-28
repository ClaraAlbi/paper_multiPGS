
library(tidyverse)
library(future.batchtools)
library(xgboost)
library(caret)

phen <- readRDS("data/phenotypes_ipsych2015_2016j.rds")
cov <- readRDS("data/covariates_ipsych2015_2016j.rds") %>% select(-is_2012)


pars_xgboost <- readRDS("data/cv_subsets_ipsych2015_2016j.rds") %>% 
  mutate(phen_out = map(p, ~case_when(.x == "F2100" ~ c("phen58"),
                                      .x == "F3100" ~ c("phen59"),
                                      .x == "F3000" ~ c("phen60", "phen61", "phen54"),
                                      .x == "F8101" ~ c("phen56"),
                                      .x == "F9201" ~ c("phen55"),
                                      .x == "F5100" ~ c("phen57")))) %>%
  expand_grid(includes = c("all", "out")) %>%
  mutate(phen_out = map2(includes, phen_out, ~case_when(.x == "out" ~ .y,
                                                        TRUE ~ NA_character_)))
    
NCORES <- 1
plan(batchtools_slurm(resources = list(
  t = "0-00:30", c = NCORES, mem = "16g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))

furrr::future_pmap(pars_xgboost[,c(1,3:7)], function(p, cv, ind.train, ind.test, phen_out, includes) {
  library(tidyverse)
  if (is.na(phen_out)) prs <- readRDS("data/prs_all_std_v2.rds") %>% select(-phen12, -phen67, -GCST008754, -GCST008753, -GCST008752, -GCST008751)
  if (!is.na(phen_out)) prs <- readRDS("data/prs_all_std_v2.rds") %>% select(-phen12, -phen67, -all_of(phen_out))

  x <- as.matrix(cbind(cov[ind.train,], prs[ind.train,]))
  y <- phen[[p]][ind.train]

  x.test <- as.matrix(cbind(cov[ind.test,], prs[ind.test,]))
  y.test <- phen[[p]][ind.test]
  dtrain <- xgb.DMatrix(x, label = y)
  dtest <- xgb.DMatrix(x.test, label = y.test)

  bst_prs <- xgboost(dtrain, eta = 0.01, objective = "binary:logistic", nrounds = 10)
  gb_prs <- predict(bst_prs, dtest)
  
  # Just cov
  x <- as.matrix(cov[ind.train,])
  x.test <- as.matrix(cov[ind.test,])
  dtrain <- xgb.DMatrix(x, label = y)
  dtest <- xgb.DMatrix(x.test, label = y.test)
  bst_cov <- xgboost(dtrain, eta = 0.01, objective = "binary:logistic", nrounds = 10)
  gb_cov <- predict(bst_cov, dtest)

  
  saveRDS(list(p = p,
               cv = cv,
               model_prs = list(bst_prs),
               pred_prs = list(gb_prs),
               model_cov = list(bst_cov),
               pred_cov = list(gb_cov),
               y.test = list(y.test)), paste0("results/xgboost/", p, "_", cv, "_", includes, ".rds"))
  })


# Gather result

library(bigparallelr)
all_res_files <- list.files("results/xgboost", full.names = T)

registerDoParallel(cl <- makeCluster(36))

res_l <- foreach(res_file = all_res_files) %dopar% readRDS(res_file)

res_xgboost <- res_l %>%
  pivot_longer(contains("pred"), names_to = c(".value", "set", "type"),
               names_pattern = "(.*)_(.*)_(.*)") %>%
  mutate(r2 = map2_dbl(pred, y.test, ~R2(.x, .y)),
         auc = map2_dbl(pred, y.test, ~bigstatsr::AUC(.x, .y))) %>%
  select(p, cv, set, type, r2, auc)

#saveRDS(res_xgboost, "results/pred_xgboost.rds")

#Features

pars <- res_xgboost %>%
  pivot_longer(contains("model"), names_to = c(".value", "set", "type"),
               names_pattern = "(.*)_(.*)_(.*)") %>%
  mutate(impRaw = map(model, ~xgb.importance(model = .x)),
         impClean = map(impRaw, ~.x[,`:=`(Cover=NULL, Frequency=NULL)]),
         n = map_dbl(impClean, ~nrow(.x)),
         vars_i = map(impClean, ~.x$Feature)) %>%
  group_by(p, set, type) %>%
  summarise(v = list(vars_i),
            mean_n = mean(n)) %>%
  mutate(c_vars = map(v, ~Reduce(intersect, .x))) %>%
  left_join(icd, by = c("p" = "Variable"))

library(tidytext)
sumstats <- readRDS("../gwas_prs/results/all_sumstats.rds")
icd <- readRDS("data/phenotype_names_ipsych2015_2016j.rds") %>% select(Variable, diagnosis, origin)

res_feature <- res %>%
  pivot_longer(contains("model"), names_to = c(".value", "set", "type"),
               names_pattern = "(.*)_(.*)_(.*)") %>%
  filter(set == "prs" & type == "t") %>%
  left_join(pars) %>%
  mutate(impRaw = map(model, ~xgb.importance(model = .x)),
         impClean = map(impRaw, ~.x[,`:=`(Cover=NULL, Frequency=NULL)]),
         vars = map(impClean, ~tibble(ind = .x$Feature, gain = .x$Gain)),
         common_vars = map2(vars, c_vars, ~.x %>% filter(ind %in% .y)),
         n = map_dbl(vars, ~nrow(.x)),
         n_common = map_dbl(common_vars, ~nrow(.x))) %>%
  select(p, cv, common_vars, n, n_common, mean_n, diagnosis) %>%
  group_by(cv) %>%
  unnest(common_vars) %>%
  ungroup() %>%
  filter(!ind %in% c("sexM", "birth_year", paste0("PC", 1:20))) %>%
  mutate(ind_2 = case_when(!ind %in% c("sexM", "birth_year", paste0("PC", 1:20)) ~ map_chr(ind, ~sumstats$`Reported trait`[sumstats$id == .x][1]),
                           TRUE ~ ind),
         ind_2 = str_sub(ind_2, end = 15),
         ind_2 = reorder_within(ind_2, gain, p),
         N_sumstats = case_when(!ind %in% c("sexM", "birth_year", paste0("PC", 1:20)) ~ map_dbl(ind, ~sumstats$N[sumstats$id == .x][1]),
                           TRUE ~ NA_real_),
         pheno = case_when(p == "skizo2015I" ~ "SCZ",
                           p == "skizospek2015I" ~ "SCZ_spectrum",
                           p == "bipol2015I" ~ "BD",
                           p == "affek2015I" ~ "AFF",
                           p == "autism2015I" ~ "ASD",
                           p == "adhd2015I" ~ "ADHD",
                           p == "anorek2016I" ~ "Ano",
                           str_detect(diagnosis, "substance") ~ "Substance abuse d",
                           diagnosis == "Fxxxx" ~ "Any psychiatric ",
                           p == "postpartum2015I" ~ "Post-partum d",
                           p == "kontrol2015I" ~ "Control status",
                           str_detect(diagnosis, "stress") ~ "Neurotic d",
                           TRUE ~ diagnosis))

#saveRDS(res_feature, "results/weights_xgboost.rds")
