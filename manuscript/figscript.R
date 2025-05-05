
# Setup ----

library(tidyverse)
library(ggh4x)
library(patchwork)

parametadata <- read_csv("metadata/parametadata.csv") %>%
  # filter(IS_type=="noIS") %>%
  filter(polarity=="pos") %>%
  filter(samp_type=="Smp") %>%
  mutate(startime=factor(startime, levels=c("Morn", "Eve"))) %>%
  mutate(depth=factor(depth, levels=c("Surf", "Deep"))) %>%
  mutate(amendment=factor(amendment, levels=c("Amm", "Nit", "Arg", "GMP"))) %>%
  arrange(amendment, depth, timepoint_num, startime) %>%
  mutate(timepoint=fct_inorder(timepoint))
peak_areas <- read_csv("targeted/all_peak_areas_w_isos_incl_13C_for_inorg.csv") %>%
  left_join(parametadata %>% distinct(filename, shortname))
all_stans <- read_csv("metadata/all_stans.csv")
samp_areas <- parametadata %>%
  filter(samp_type=="Smp") %>%
  distinct(shortname, depth) %>%
  left_join(peak_areas) %>%
  mutate(shortname=fct_inorder(shortname)) %>%
  filter(!compound_name%in%c("Melamine", "Uridine", "Asterina-330?", "Palythene/Usujirene?",
                             "Mycosporine-2-glycine?", "Porphyra-334?", "Itaconic acid?",
                             "Inosine", "Dexpanthenol", "Choline sulfate", "5-Hydroxyectoine")) %>%
  filter(!(depth=="Deep" & compound_name%in%c("Also Palythene/Usujirene?", "Shinorine?",
                                              "Mycosporine-glycine?", "Palythine?"))) %>%
  # filter(!str_detect(shortname, "Nit-Deep-Eve-T10")) %>%
  select(-depth)
labeled_fracs <- samp_areas %>%
  mutate(labeled_frac=1-area/sum(area), .by = c("shortname", "compound_name")) %>%
  select(compound_name, iso_name, shortname, labeled_frac) %>%
  filter(str_detect(iso_name, "13C0, 15N0")) %>%
  select(-iso_name)
n_labeled_fracs <- samp_areas %>%
  filter(str_detect(iso_name, "13C0")) %>%
  filter(!str_detect(iso_name, "15N0") | str_detect(iso_name, "13C0, 15N0")) %>%
  mutate(labeled_frac=1-area/sum(area), .by = c("shortname", "compound_name")) %>%
  select(compound_name, iso_name, shortname, labeled_frac) %>%
  filter(str_detect(iso_name, "13C0, 15N0")) %>%
  select(-iso_name)
log_coef_df <- labeled_fracs %>%
  left_join(parametadata %>% distinct(shortname, amendment, depth, timepoint_num)) %>%
  nest(data=c(shortname, timepoint_num, labeled_frac)) %>%
  mutate(log_fit=map(data, ~glm(.x$labeled_frac~.x$timepoint_num, family="quasibinomial")),
         coefs=map(log_fit, function(x)broom::tidy(x))) %>%
  unnest(coefs) %>%
  select(compound_name, amendment, depth, term, estimate) %>%
  pivot_wider(names_from = term, values_from = estimate) %>%
  rename(intercept=`(Intercept)`, slope=`.x$timepoint_num`) %>%
  mutate(half_max=-intercept/slope, scale=1/slope) %>%
  mutate(frac_at_T24 = 1/(1+exp(-(intercept+slope*24))))
splividis <- function(n, ...){
  vir <- viridis::viridis(n, ...)
  if(n%%2==0){
    vir[as.numeric(matrix(1:n, ncol=2, byrow = TRUE))]
  } else {
    index_vec <- suppressWarnings(as.numeric(matrix(1:n, ncol=2, byrow = TRUE)))
    vir[head(index_vec, -1)]
  }
}

