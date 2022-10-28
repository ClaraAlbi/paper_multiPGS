
library(tidyverse)
library(glmnet)
library(caret)
library(future.batchtools)
library(glue)

### WITHOUT IPSYCH AS COVARIATE

phen <- readRDS("data/phenotypes_ipsych2015_2016j.rds")
cov <- readRDS("data/covariates_ipsych2015_2016j.rds") 

pars_lasso <- readRDS("data/cv_subsets_ipsych2015_2016j.rds") %>% 
  mutate(phen_out = map(p, ~case_when(.x == "F2100" ~ c("phen58"),
                              .x == "F3100" ~ c("phen59"),
                              .x == "F3000" ~ c("phen60", "phen61", "phen54"),
                              .x == "F8101" ~ c("phen56"),
                              .x == "F9201" ~ c("phen55"),
                              .x == "F5100" ~ c("phen57"),
                              TRUE ~ c("")))) %>%
  expand_grid(includes = c("all", "out")) %>%
  mutate(phen_out = map2(includes, phen_out, ~case_when(.x == "out" ~ .y,
                                                        TRUE ~ NA_character_))) %>%
  filter(!file.exists(glue::glue("results/lasso/{p}_{cv}.rds"))) 

vars_cov <- c("birth_year", "sexM", paste0("PC", 1:20))

NCORES <- 4
plan(batchtools_slurm(resources = list(
  t = "0-12:00", c = NCORES, mem = "16g",
  name = basename(rstudioapi::getSourceEditorContext()$path))))


furrr::future_pmap(pars_lasso[,c(1,3:7)], function(p, cv, ind.train, ind.test, phen_out, includes) {
  library(tidyverse)
  if (is.na(phen_out)) prs <- readRDS("data/prs_all_std_v2.rds") %>% select(-phen12, -phen67, -GCST008754, -GCST008753, -GCST008752, -GCST008751)
  if (!is.na(phen_out)) prs <- readRDS("data/prs_all_std_v2.rds") %>% select(-phen12, -phen67, -GCST008754, -GCST008753, -GCST008752, -GCST008751, -all_of(phen_out))
  
  
  x <- as.matrix(cbind(cov[ind.train,], prs[ind.train,]))
  y <- phen[[p]][ind.train]
  x.test <- as.matrix(cbind(cov[ind.test,], prs[ind.test,]))
  y.test <- phen[[p]][ind.test]
  
  glm_cov <- glm(y ~ ., data = cov[ind.train, vars_cov], family = binomial())
  
  p.fac <- c(rep(0, ncol(cov)), rep(1, ncol(prs)))
  cv.lasso <- cv.glmnet(x, y, alpha = 1, penalty.factor = p.fac, family = "binomial")
  
  coefs <- coef(cv.lasso, s = cv.lasso$lambda.1se)
  idx <- coefs[-1,1] != 0
  idx[3] <- FALSE
  
  pred_lasso <- scale(x.test[,idx] %*% as.vector(coefs[-1,1][idx]))[,1]
  pred_cov <- scale(x.test[,c(1:2, 4:23)] %*% as.vector(glm_cov$coefficients[-1]))[,1]
  
  saveRDS(list(p = p,
               cv = cv,
               cv.lasso = list(cv.lasso),
               glm_cov = list(glm_cov),
               pred_lasso = list(pred_lasso),
               pred_cov = list(pred_cov),
               y.test = list(y.test)), paste0("results/lasso/", p, "_", cv, "_", includes, ".rds"))
})


# Gather results
library(bigparallelr)
all_res_files <- list.files("results/lasso", full.names = T)

registerDoParallel(cl <- makeCluster(36))

res <- foreach(res_file = all_res_files) %dopar% readRDS(res_file)
res_lasso <- do.call(bind_rows, res) %>%
  mutate(name = all_res_files,
         name = str_match(name, glue("results/lasso/{p}_{cv}_*(.*?).rds"))[,2])
stopCluster(cl)


# Extract weights
library(tidytext)
sumstats <- readRDS("../gwas_prs/results/all_sumstats.rds")
icd <- readRDS("data/phenotype_names_ipsych2015_2016j.rds") %>% select(Variable, diagnosis, origin)
  
pars <- res_lasso %>%
  mutate(c = map(cv.lasso, ~coef(.x, s=.x$lambda.min)),
         vars = map(c, ~tibble(stack(.x[.x[,1]!=0,])) %>% tail(-24)),
         vars_i = map(vars, ~.x[["ind"]]),
         n = map_dbl(vars, ~nrow(.x))) %>%
  group_by(p) %>%
  summarise(v = list(vars_i),
            mean_n = mean(n)) %>%
  mutate(c_vars = map(v, ~Reduce(intersect, .x))) %>%
  left_join(icd, by = c("p" = "Variable"))

text_ipsych <- res_lasso %>%
  inner_join(pars) %>%
  mutate(c = map(cv.lasso, ~ coef(.x, s=.x$lambda.1se)),
         na = map(c, ~.x@Dimnames[[1]]),
         v = map(c, ~as.vector(.x)),
         d = map2(na, v, ~tibble(n = as.factor(.x), v = .y))) %>%
  select(p, d, diagnosis) %>%
  unnest(d) %>%
  filter(!n %in% c("(Intercept)", "is_2012", paste0("PC", 1:20))) %>%
  group_by(p, n, diagnosis) %>%
  summarise(w_m = mean(v),
            w_sd = sd(v)) %>%
  filter(w_m != 0) %>%
  mutate(ind = map_chr(n, ~case_when(.x %in% c("sexM", "birth_year") ~ as.character(.x),
                                        !.x %in% c("sexM", "birth_year") ~ sumstats$`Reported trait`[sumstats$id == .x][1])),
         pheno = case_when(p == "skizo2015I" ~ "SCZ",
                           p == "bipol2015I" ~ "BD",
                           p == "affek2015I" ~ "AFF",
                           p == "autism2015I" ~ "ASD",
                           p == "adhd2015I" ~ "ADHD",
                           p == "anorek2016I" ~ "Ano",
                           str_detect(diagnosis, "substance") ~ "Substance abuse disorders",
                           diagnosis == "Fxxxx" ~ "Any psychiatric disorder",
                           p == "postpartum2015I" ~ "Post-partum disorder",
                           str_detect(diagnosis, "stress") ~ "Neurotic, stress-related and somatoform",
                           TRUE ~ diagnosis)) %>%
  ungroup() 

#saveRDS(text_ipsych, "results/weights_lasso_ipsych.rds")
