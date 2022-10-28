# compare
library(broom.helpers)
library(tidyverse)
library(glue)

n <- tribble(~p, ~id, ~prev, ~tag, ~comp,
             "F2100", "phen58", 0.01, "SCZ", "PGC3-SCZ",
             "F3100", "phen59", 0.01, "BD", "PGC3-BD",
             "F3000", "phen61", 0.08, "AFF", "MDD-Wray2018",
             "F8101", "phen56", 0.01, "ASD", "PGC-ASD",
             "F9201", "phen55", 0.05, "ADHD", "PGC-ADHD",
             "F5100", "phen57", 0.01, "AN", "PGC-AN",
             "F1003", "phen66", 0.015, "F1-SUD", "PGC-CUD",
             "F4100", "phen30", 0.1, "F4-Neuro-stress", "Neuroticism-Okbay2016",
             "F6100", "phen61", 0.09, "F6-Personality", "MDD-Wray2018",
             "F7000", "phen56", 0.03, "F7-Mental ret", "PGC-ASD",
             "Fxxxx", "phen60", 0.2, "FX-Any", "MDD-Howard2019",
             "Fxxx1", "phen60", 0.15, "FX-Postpartum","MDD-Howard2019",
             "G4010", "GCST007349", 0.02, "G4-Epilepsy", "Partial epilepsy-ILAE",
             "G4030", "GA3603", 0.01, "G4-Migraine", "Migraine-UKB",
             "J4501", "phen51", 0.034, "J4-Asthma", "Asthma-Ferreira2017",
             "K4010", "GA3671", 0.017, "K4-Hernia", "Hernia-UKB",
             "R5010", "GCST007349", 0.04, "R5-FebSeizure", "Partial epilepsy-ILAE")

phen <- readRDS("data/phenotypes_ipsych2015_2016j.rds")
cov <- readRDS("data/covariates_ipsych2015_2016j.rds") 

prs <- readRDS("data/prs_all_std_v2.rds")

prs_metaPRS <- data.table::fread("../to_DST/metaPRS_6PD_iPSYCH15.prs")

ipsych <- bigsnpr::snp_attach("~/iPSYCH2015/HRC_Imputed/bigsnp_r_format/ipsych_hm3.rds")$fam$family.ID
idx <- match(ipsych, prs_metaPRS$pid)

res_meta <- tibble(prs_metaPRS[idx,]) %>%
  summarise(across(c(contains("bolt"), contains("meta")), list)) %>%
  pivot_longer(everything()) %>%
  mutate(tag = map_chr(name, ~tail(str_split(.x, "_")[[1]], n = 1)),
         tag = case_when(tag == "Ano" ~ "AN",
                         tag == "MDD" ~ "AFF",
                         tag == "BP" ~ "BD",
                         TRUE ~ tag),
         name = map_chr(name, ~str_split(.x, "_")[[1]][1])) %>%
  pivot_wider(names_from = name, values_from = value) %>%
  rename(pred_blup = prs, pred_metaPRS = metaPRS) %>% select(pred_blup, pred_metaPRS, tag)

res_mtblup <- tibble(f = list.files("results/smtpred", pattern = "score", full.names = T)) %>%
  mutate(p = str_match(f, "results/smtpred/*(.*?).score")[,2],
         pred_mtblup = map(f, ~read_table2(.x)[[3]])) %>% select(p, pred_mtblup)

pred_cv <- readRDS("data/cv_subsets_ipsych2015_2016j.rds") %>% 
  filter(str_detect(p, "F9201|F8101|F3100|F2100|F3000|F5100")) %>%
  select(p, cv, ind.test, ind.train) %>%
  left_join(res_lasso) %>% 
  left_join(res_xgboost) %>%
  mutate(i = map(cv.lasso, ~ which(abs(coef(.x, s = "lambda.min"))[c(-1:-24)] >= 0.01)),
         names_coef01 = map2(cv.lasso, i, ~ coef(.x, s = "lambda.min")@Dimnames[[1]][c(-1:-24)][.y]),
         weight_coef01 = map2(cv.lasso, i, ~ coef(.x, s = "lambda.min")[c(-1:-24)][.y]),
         pred_lassotr = map2(cv.lasso, ind.train, ~scale(as.matrix(prs[.y, coef(.x, s = "lambda.min")@Dimnames[[1]][c(-1:-24)]]) %*% as.vector(coef(.x, s = "lambda.min")[c(-1:-24)]))[,1]),
         pred_lassot = map2(cv.lasso, ind.test, ~scale(as.matrix(prs[.y, coef(.x, s = "lambda.min")@Dimnames[[1]][c(-1:-24)]]) %*% as.vector(coef(.x, s = "lambda.min")[c(-1:-24)]))[,1]),
         pred_lMP01 = pmap(list(ind.train, names_coef01, weight_coef01), function(a,b,c) scale(as.matrix(prs[a,b]) %*% as.vector(c))[,1]),
         pred_lMP01t = pmap(list(ind.test, names_coef01, weight_coef01), function(a,b,c) scale(as.matrix(prs[a,b]) %*% as.vector(c))[,1])
  ) %>%
  left_join(n) %>%
  left_join(res_meta) %>%
  left_join(res_mtblup) %>%
  mutate(y.train = map2(p, ind.train, ~phen[[.x]][.y]),
         data = map2(pred_lMP01, ind.train, ~bind_cols(.x, cov[.y,])), # DATA FOR SINGLE PGS
         datatr = map2(pred_lassotr, ind.train, ~bind_cols(.x, cov[.y,])),
         data_b = map2(pred_mtblup, ind.train, ~bind_cols(.x, cov)[.y,]), # DATA FOR BLUP PGS
         data_M = map2(pred_metaPRS, ind.train, ~bind_cols(.x, cov)[.y,]), #DATA META PGS
         data_test = map2(pred_lMP01t, ind.test, ~bind_cols(.x, cov[.y,])), #TEST SINGLE PGS
         data_testtr = map2(pred_lassot, ind.test, ~bind_cols(.x, cov[.y,])),
         data_test_b = map2(pred_mtblup, ind.test, ~bind_cols(.x, cov)[.y,]), #TEST BLUP PGS
         data_test_M = map2(pred_metaPRS, ind.test, ~bind_cols(.x, cov)[.y,]), #TEST BLUP PGS
         model_prs = map2(data, y.train, ~glm(.y ~ ., data = .x, family = binomial())), #model SINGLE
         model_tr = map2(datatr, y.train, ~glm(.y ~ ., data = .x, family = binomial())),
         model_b = map2(data_b, y.train, ~glm(.y ~ ., data = .x, family = binomial())), #MODEL blup
         model_M = map2(data_M, y.train, ~glm(.y ~ ., data = .x, family = binomial())), #MODEL META
         #pred_O = map2(data_test, model_prs, ~ scale(as.matrix(.x[,c(1:3, 5:24)]) %*% as.vector(.y$coefficients[c(-1, -5)]))[,1]),
         #pred_Onocov = map2(data_test, model_prs, ~ scale(as.matrix(.x[,1]) %*% as.vector(.y$coefficients[2]))[,1]),
         pred_blupnocov = map2(pred_blup, ind.test, ~scale(.x[.y])[,1]),
         pred_metaprsnocov = map2(pred_metaPRS, ind.test, ~scale(.x[.y])[,1]),
         pred_lassoC = map2(data_testtr, model_tr, ~ scale(as.matrix(.x[,c(1:3, 5:24)]) %*% as.vector(.y$coefficients[c(-1, -5)]))[,1]),
         pred_mtblup = map2(data_test_b, model_b, ~ scale(as.matrix(.x[,c(1:3, 5:24)]) %*% as.vector(.y$coefficients[c(-1, -5)]))[,1]),
         pred_metaprs = map2(data_test_M, model_M, ~ scale(as.matrix(.x[,c(1:3, 5:24)]) %*% as.vector(.y$coefficients[c(-1, -5)]))[,1]),
         pred_lMP01c = map2(data_test, model_prs, ~ scale(as.matrix(.x[,c(1:3, 5:24)]) %*% as.vector(.y$coefficients[c(-1, -5)]))[,1])
  ) %>% 
  select(-datatr, -data_testtr, -model_tr, -pred_lMP01, -y.train, -model_prs, -data, -data_test, -model_M, -data_M, -data_test_M, -ind.train, -data_test_b, -data_b, -model_b, -pred_blup, -pred_metaPRS, -pred_lassotr, -pred_lassot) %>%
  pivot_longer(contains("pred"), names_to = c(".value", "type"),
               names_pattern = "(.*)_(.*)")

#saveRDS(pred_cv, "results/prs_ipsych_com.rds")


z <- pred_cv %>%
  filter(!type %in% c("")) %>%
  mutate(r2 = map2_dbl(pred, y.test, ~R2(.x, .y)),
         auc = map2_dbl(pred, y.test, ~bigstatsr::AUC(.x, .y)),
         m_case = map_dbl(y.test, ~mean(.x)),
         coef_liab = map2_dbl(prev, m_case, ~bigsnpr::coef_to_liab(K_pop = .x, K_gwas = .y)),
         r2 = r2 * coef_liab) %>%
  pivot_wider(c(p, cv, id, prev, m_case, tag, comp), names_from = type, values_from = c(r2, auc)) %>%
  pivot_longer(c(contains("r2"), contains("auc"), -c(r2_cov, auc_cov)), names_to = c(".value", "model"), names_sep = "_") %>%
  mutate(adj_r2 = (r2 - r2_cov) / (1 - r2_cov),
         delta_auc = (auc - auc_cov) / auc_cov) %>%
  group_by(p, id, tag, prev, comp, model) %>%
  summarise(across(c(adj_r2, delta_auc, r2, auc), list(mean = mean, 
                                                       ci = function(t) map(list(t), ~quantile(sample(.x, 10e3, replace = T), probs = c(0.025, 0.975)))),
                   .names = "{.col}.{.fn}")) %>%
  unnest_wider(adj_r2.ci, names_sep = "_") %>%
  unnest_wider(delta_auc.ci, names_sep = "_") %>%
  unnest_wider(r2.ci, names_sep = "_")  %>%
  unnest_wider(auc.ci, names_sep = "_") 

x <- pred_cv %>%
  mutate(q10 = map(pred, ~as.factor(paste0("q", ntile(.x, 10)))),
         data_q10 = map2(ind.test, q10, ~tibble(cov[.x,-3], q10 = .y)),
         q5 = map(pred, ~as.factor(paste0("q", ntile(.x, 5)))),
         data_q5 = map2(ind.test, q5, ~tibble(cov[.x,-3], q5 = .y)),
         or_q1 = map2(data_q10, y.test,~glm(.y ~ ., data = .x, family = "binomial")),
         or_qm = map2(data_q10, y.test, ~glm(.y ~ ., data = .x %>% mutate(q10 = relevel(fct_recode(q10, "q5" = "q6"), ref = "q5")), family = "binomial")),
         or_q5 = map2(data_q5, y.test,~glm(.y ~ ., data = .x, family = "binomial")),
         or_qm5 = map2(data_q5, y.test, ~glm(.y ~ ., data = .x %>% mutate(q5 = relevel(q5, ref = "q3")), family = "binomial")),
         q1_t = map(or_q1, ~bind_cols(broom::tidy(.x), confint.default(.x) %>% as_tibble() %>% set_names(c("conf.low", "conf.high")))  %>% add_row(term = "q10q1", estimate = 0)),
         qm_t = map(or_qm, ~bind_cols(broom::tidy(.x), confint.default(.x) %>% as_tibble() %>% set_names(c("conf.low", "conf.high"))) %>% add_row(term = "qmid", estimate = 0)),
         q5_t = map(or_q5, ~bind_cols(broom::tidy(.x), confint.default(.x) %>% as_tibble() %>% set_names(c("conf.low", "conf.high")))  %>% add_row(term = "q5q1", estimate = 0)),
         qm5_t = map(or_qm5, ~bind_cols(broom::tidy(.x), confint.default(.x) %>% as_tibble() %>% set_names(c("conf.low", "conf.high"))) %>% add_row(term = "qmid", estimate = 0))
  ) %>%
  pivot_longer(c(q1_t, qm_t, q5_t, qm5_t), names_to = "stat") %>%
  unnest(value) %>%
  filter(str_detect(term, "q")) %>%
  select(p, cv, id, prev, tag, comp, stat, type, term, estimate, conf.low, conf.high, p.value) %>%
  group_by(p, id, prev, tag, comp, stat, type, term) %>%
  summarise(across(c(estimate), list(mean = ~mean(.x, na.rm = T), 
                                     ci = function(t) map(list(t), ~quantile(sample(.x, 10e3, replace = T), probs = c(0.025, 0.975), na.rm = T))),
                   .names = "{.col}.{.fn}")) %>%
  unnest_wider(estimate.ci, names_sep = "_")




saveRDS(z, "results/pred_cv_ipsych_n.rds")
saveRDS(x, "results/pred_q10_ipsych_n.rds")


### TOP 10 SINGLE R2

sumstats <- readRDS("../gwas_prs/results/all_sumstats.rds")

res <- readRDS("data/cv_subsets_ipsych2015_2016j.rds") %>%
  inner_join(readRDS("results/weights_lasso_ipsych.rds")) %>%
  inner_join(n) %>%
  filter(!n %in% c("sexM", "birth_year")) %>%
  group_by(p, cv, pheno) %>%
  slice_max(order_by = abs(w_m), n = 10) %>%
  ungroup() %>%
  #slice(1:5) %>%
  mutate(ind = map_chr(ind, ~str_split(.x, "\\(|- |or |\\/")[[1]][1]),
         N_p1 = map_dbl(p, ~4/(1/sum(phen[[.x]] == 1, na.rm = T) + 1/sum(phen[[.x]] == 0, na.rm = T))),
         N_p2 = map_dbl(n, ~case_when(is.na(sumstats$Ncase[sumstats$id == .x]) ~ sumstats$N[sumstats$id == .x],
                                      !is.na(sumstats$Ncase[sumstats$id == .x]) ~ 4/(1/sumstats$Ncase[sumstats$id == .x] + 1/sumstats$Ncontrol[sumstats$id == .x]))),
         N_p2r = round(N_p2/1000, 0),
         y.train = map2(p, ind.train, ~phen[[.x]][.y]),
         y.test = map2(p, ind.test, ~phen[[.x]][.y]),
         data = map2(n, ind.train, ~bind_cols(prs[[.x]], cov)[.y,]),
         data_test = map2(n, ind.test, ~bind_cols(prs[[.x]], cov)[.y,]), #TEST SINGLE PGS
         model = map2(data, y.train, ~glm(.y ~ ., data = .x, family = binomial())), #model SINGLE
         pred_prs = map2(data_test, model, ~ scale(as.matrix(.x[,c(1:3, 5:24)]) %*% as.vector(.y$coefficients[c(-1, -5)]))[,1]),
         pred_cov = map2(data_test, model, ~ scale(as.matrix(.x[,c(2:3, 5:24)]) %*% as.vector(.y$coefficients[c(-1, -2, -5)]))[,1])) %>%
  pivot_longer(contains("pred"), names_to = c(".value", "type"),
               names_pattern = "(.*)_(.*)") %>%
  mutate(
         r2 = map2_dbl(pred, y.test, ~R2(.x, .y)),
         auc = map2_dbl(pred, y.test, ~bigstatsr::AUC(.x, .y)),
         m_case = map_dbl(y.test, ~mean(.x)),
         coef_liab = map2_dbl(prev, m_case, ~bigsnpr::coef_to_liab(K_pop = .x, K_gwas = .y)),
         r2 = r2 * coef_liab
         ) %>%
  pivot_wider(c(p, cv, n, diagnosis, w_m, w_sd, ind, pheno, id, tag, N_p1, N_p2, N_p2r, r2, auc), names_from = type, values_from = c(r2, auc)) %>%
  pivot_longer(c(contains("r2"), contains("auc"), -c(r2_cov, auc_cov)), names_to = c(".value", "model"), names_sep = "_") %>%
  mutate(adj_r2 = (r2 - r2_cov) / (1 - r2_cov),
         delta_auc = (auc - auc_cov) / auc_cov) %>%
  group_by(p, n, diagnosis, w_m, w_sd, ind, pheno, id, tag, N_p1, N_p2, N_p2r) %>%
  summarise(across(c(adj_r2, delta_auc, r2, auc), list(mean = mean, 
                                                       ci = function(t) map(list(t), ~quantile(sample(.x, 10e3, replace = T), probs = c(0.025, 0.975)))),
                   .names = "{.col}.{.fn}")) %>%
  unnest_wider(adj_r2.ci, names_sep = "_") %>%
  unnest_wider(delta_auc.ci, names_sep = "_") %>%
  unnest_wider(r2.ci, names_sep = "_")  %>%
  unnest_wider(auc.ci, names_sep = "_") %>%
  filter(file.exists(paste0("results/ldsc/", p, "_", n, ".log"))) %>%
  mutate(filename = paste0("results/ldsc/", p, "_", n, ".log"),
         pars = map(filename, ~data.table::fread(.x, skip = "p1", fill = T, nrows = 1))) %>%
  rename(phen = p) %>%
  unnest(pars) 

saveRDS(res, "results/pred_cv_single.rds")


