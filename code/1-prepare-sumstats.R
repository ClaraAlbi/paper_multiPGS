library(dplyr)
library(purrr)
library(threadr)
library(bigsnpr)
library(ggplot2)

# GWAS CATALOG

gwas_catalog <- read_tsv("files/gwas_catalog_v1.0.2-studies_r2020-09-09.tsv")
gwas_sumstats <- read_csv("files/list_gwas_summary_statistics_14_Sep_2020.csv")
gwas_ancestry <- read_tsv("files/gwas-catalog-v1.0.3-ancestries-r2020-09-09.tsv")
efo_map <- read_tsv("files/gwas_catalog_trait-mappings_r2020-09-09.tsv") %>%
  group_by(`EFO URI`, `EFO term`, `Parent term`) %>%
  summarise(subphens = n()) %>%
  mutate(Cat = case_when(str_detect(`Parent term`, "measurement") ~ "Measurement",
                         str_detect(`Parent term`, "disease|disorder|Cancer") ~ "Disease",
                         `Parent term` %in% c("Biological process", "Other trait") ~ "Other"))

sumstats <- gwas_sumstats %>%
  filter(str_detect(`Data access`, "API")) %>% # Formatted sumstats (theoretically)
  left_join(gwas_catalog, by = c("Study accession" = "STUDY ACCESSION")) %>% #More info
  left_join(gwas_ancestry, by = c("Study accession" = "STUDY ACCESSION", intersect(names(gwas_catalog), names(gwas_ancestry))[-3])) %>%
  filter(`BROAD ANCESTRAL CATEGORY` == "European") %>% # Containing european samples %>%
  mutate(`STAGE` = case_when(`NUMBER OF INDIVIDUALS` == 58655 ~ "replication",
                             TRUE ~ `STAGE`)) %>% # Fix error in reporting the study
  filter(`BROAD ANCESTRAL CATEGORY` == "European") %>% # Containing european samples
  pivot_wider(c(1:29, 32), names_from = `STAGE`, values_from = `NUMBER OF INDIVIDUALS`) %>%
  filter(!is.na(initial)) %>%
  filter(!str_detect(`INITIAL SAMPLE DESCRIPTION`, "founder")) %>% # Founder studies
  group_by(`Reported trait`) %>%
  slice(which.max(ymd(`Publication date`))) %>% # Latest catalog entry
  separate(`INITIAL SAMPLE SIZE`, into = c("Nca_i", "Nco_i"), sep = "cases", remove = FALSE) %>%
  separate(`REPLICATION SAMPLE SIZE`, into = c("Nca_r", "Nco_r"), sep = "cases", remove = FALSE) %>%
  mutate_at(vars(Nca_i, Nco_i, Nca_r, Nco_r), parse_number)  %>%
  mutate(Nca_i = case_when(is.na(Nco_i) ~ NA_real_,TRUE ~ Nca_i),
         Nca_r = case_when(is.na(Nco_r) ~ NA_real_,TRUE ~ Nca_r)) %>%
  left_join(efo_map, by = c("MAPPED_TRAIT_URI" = "EFO URI")) %>%
  arrange(Cat, `Parent term`, `EFO term`) %>%
  mutate(files_o = map(`FTP Path`, ~ unlist(list_files_ftp(.x), recursive = TRUE)),
         ss_f = map(files_o, ~str_detect(.x, "/$|harmonised|README|md5sum|xlsx|readme|Readme|ReadMe|Read-me|Read_Me|read.me|tbi|doc|md5|pdf|FieldDescription|original|rtf|Chinese|Female|female|FEMALE|Women|women|WOMEN|male|Male|MALE|men|Men|MEN|1000G", negate = TRUE)),
         total = map(ss_f, sum),
         ss = map2(files_o, ss_f, ~.x[.y]),
         files = map(`FTP Path`, ~ unlist(list_files_ftp(paste0(.x, "/harmonised")), recursive = TRUE)),
         formatted = map_chr(files, ~.x[1])) %>%
  filter(initial > 10000) %>%
  filter(total == 1 | !is.na(formatted)) %>%
  mutate(file = case_when(!is.na(formatted) ~ formatted,
                          TRUE ~ as.character(ss)),
         ext = map_chr(file, ~tail(str_split(.x, "\\.")[[1]], n = 1))) %>%
  filter(ext %in% c("csv", "gz", "meta", "tbl", "tsv", "txt", "zip")) %>%
  select(1:8, 36, 21, 22, 35, 24, 25, 37:40, file, ext) %>%
  mutate(source = "gwas_catalog", folder = NA) %>%
  rename(N = initial, Ncase = Nca_i, Ncontrol = Nco_i) %>%
  select(1:11, 19:22)


# GWAS ATLAS

gwas_atlas <- data.table::fread("files/gwasATLAS_v20191115.txt.gz") %>%
  filter(Population %in% c("EUR", "UKB2 (EUR)")) %>%
  #filter(!PMID %in% catalog_traits$`PubMed ID`) %>%
  filter(N > 1e4 & Nsnps > 250e3) %>%
  filter(Year > 2014) %>%
  group_by(Trait) %>%
  mutate(n = n()) %>%
  filter(n < 2) %>%
  ungroup() %>%
  mutate(f = map(File, ~str_split(.x, "/", simplify = TRUE)),
         f2 = map_chr(f, ~.x[length(.x)]),
         ext = map_chr(f2, ~tail(str_split(.x, "\\.")[[1]], n = 1))) %>%
  filter(ext %in% c("gz", "zip", "csv", "txt")) %>%
  mutate(`First Author` = "",
         `Journal` = "",
         id = paste0("GA", id),
         source = "GWASAtlas",
         folder = "gwas_atlas",
         `Publication date` = as.character(Year)) %>%
  rename(`PubMed ID` = PMID,
         `Study accession` = id,
         `Trait(s)` = Trait,
         `Reported trait` = uniqTrait,
         `EFO term` = Domain,
         `Parent term` = ChapterLevel,
         subphens = SubchapterLevel,
         file = File,
         Title = Consortium) %>%
  select(1:2,4,6, 10:11,13:15,33:38)


# Internal / NCRR TRAITS
ncrr <- read_csv2("files/ncrr_traits.csv") %>% rename(N = initial, Ncase = Nca_i, Ncontrol = Nco_i)

v1 <- ncrr %>% 
  filter(source %in% c("dropbox", "XFiles") | `First Author` == "Jones SE") %>%
  mutate(`Study accession` = case_when(!source %in% c("dropbox", "XFiles") ~ `Study accession`,
                                       TRUE ~ paste0("phen" ,row_number()))) %>%
  filter(source != "gwas_catalog") %>%
  select(1:11,19:22)

# Join all

v2 <- bind_rows(v1, catalog_traits, gwas_atlas) %>%
  mutate(ext = map_chr(file, ~tail(str_split(.x, "\\.")[[1]], n = 1))) %>%
  rename(id = `Study accession`) %>%
  select(1, id, 2, 4:15) 

write_tsv(v2, "data/ncrr_gwascatalog_traits.rds")
