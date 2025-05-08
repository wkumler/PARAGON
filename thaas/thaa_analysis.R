

thaa_mzs <- read_csv("thaas/AminoAcid_masses_paragon1.csv") %>%
  set_names(c("compound_name", "cmpd_mz", "cmpd_rt")) %>%
  mutate(compound_name=str_replace(compound_name, "Phenylalaine", "Phenylalanine"))

ms_files <- list.files("thaas/mzMLs", pattern = "mzML", full.names = TRUE)
msdata <- grabMSdata(ms_files, grab_what=c("EIC", "metadata"), mz = thaa_mzs$cmpd_mz, ppm = 20)
msdata$BPC %>% qplotMS1data(color_col="filename")

metathaa <- data.frame(filename=basename(ms_files)) %>%
  mutate(samp_type=str_extract(filename, "Blk|Smp|Std|Poo")) %>%
  mutate(amendment=str_extract(filename, "Amm|Nit|Arg|GMP")) %>%
  mutate(timepoint=str_extract(filename, "T\\d+")) %>%
  mutate(timepoint=factor(timepoint, levels=c("T0", "T1", "T3", "T10", "T26", "T73"))) %>%
  mutate(depth=str_extract(filename, "Surf|Deep")) %>%
  mutate(depth=factor(depth, levels=c("Surf", "Deep"))) %>%
  mutate(tripl=str_extract(filename, "[A-F](?=\\.mzML)")) %>%
  mutate(startime=str_extract(filename, "Morn|Eve")) %>%
  mutate(startime=case_when(
    samp_type=="Smp" & tripl%in%c("D", "E", "F")~"Morn",
    TRUE~startime
  )) %>%
  mutate(startime=factor(startime, levels=c("Morn", "Eve"))) %>%
  mutate(tripl=case_when(
    tripl=="D"~"A",
    tripl=="E"~"B",
    tripl=="F"~"C",
    TRUE~tripl
  )) %>%
  left_join(msdata$metadata %>% distinct(filename, timestamp))
write_csv(metathaa, "thaas/metathaa.csv")



thaa_mzs %>%
  pwalk(function(compound_name, cmpd_mz, cmpd_rt){
    gp <- msdata$MS1[mz%between%pmppm(cmpd_mz, ppm = 10)] %>%
      filter(rt%between%c(cmpd_rt-0.5, cmpd_rt+0.5)) %>%
      qplotMS1data() +
      geom_vline(xintercept = cmpd_rt, color="red") +
      ggtitle(compound_name)
    print(gp)
  },.progress = TRUE)

compound_bounds <- tribble(
  ~compound_name, ~rtstart, ~rtend,
  "Histidine", 1.3, 1.45,
  "Arginine", 1.88, 2.02,
  "Serine", 2.12, 2.32,
  "Glycine", 2.38, 2.52,
  "Aspartic acid", 2.53, 2.68,
  "Glutamic acid", 2.78, 2.9,
  "Threonine", 2.95, 3.1,
  "Alanine", 3.3, 3.45,
  "Proline", 3.7, 3.8,
  "Lysine", 4.05, 4.18,
  "Tyrosine", 4.38, 4.48,
  "Methionine", 4.52, 4.62,
  "Valine", 4.58, 4.7,
  "Isoleucine", 5.18, 5.26,
  "Leucine", 5.26, 5.35,
  "Phenylalanine", 5.35, 5.5,
  "Citrulline, 2H4", 2.65, 2.75,
  "Ornithine, 2H6", 3.82, 3.92,
  "Taurine, 2H4", 2.13, 2.23
) %>%
  left_join(thaa_mzs %>% select(-cmpd_rt))

# Check for integration bounds quality using all data
compound_bounds %>%
  pwalk(function(compound_name, rtstart, rtend, cmpd_mz){
    gp <- msdata$MS1[mz%between%pmppm(cmpd_mz, ppm = 10)] %>%
      filter(rt%between%c(rtstart-0.2, rtend+0.2)) %>%
      qplotMS1data() +
      geom_vline(xintercept = c(rtstart, rtend), color=c("green", "red")) +
      ggtitle(compound_name) +
      scale_x_continuous(breaks = seq(0, 7, 0.1))
    print(gp)
  },.progress = TRUE)


stan_map <- all_stans %>% 
  distinct(compound_name, n_N, n_C) %>%
  mutate(compound_name=str_remove(compound_name, "^L-"))
iso_mzs <- compound_bounds %>%
  distinct(compound_name, cmpd_mz) %>%
  left_join(stan_map) %>%
  mutate(n_C=case_when(
    compound_name=="Phenylalanine"~9,
    compound_name=="Taurine, 2H4"~NA,
    TRUE~n_C
  )) %>%
  mutate(n_N=case_when(
    compound_name=="Phenylalanine"~1,
    compound_name=="Taurine, 2H4"~NA,
    TRUE~n_N
  )) %>%
  select(compound_name, cmpd_mz, n_C, n_N) %>%
  replace_na(list(n_N=0, n_C=0)) %>%
  pmap(function(...){
    row_data <- list(...)
    expand_grid(n_C=0:row_data$n_C, n_N=0:row_data$n_N) %>%
      mutate(compound_name=row_data$compound_name, cmpd_mz=row_data$cmpd_mz)
  }) %>%
  bind_rows() %>%
  mutate(iso_mz=cmpd_mz+1.003355*n_C+0.997035*n_N) %>%
  mutate(iso_name = paste0(compound_name, ", 13C", n_C, ", 15N", n_N))

all_iso_areas <- iso_mzs %>%
  left_join(compound_bounds) %>%
  select(iso_name, rtstart, rtend, iso_mz) %>%
  pmap(function(iso_name, rtstart, rtend, iso_mz){
    msdata$MS1[mz%between%pmppm(iso_mz, ppm = 10)] %>%
      filter(rt%between%c(rtstart, rtend)) %>%
      slice_max(int, by=c(filename, rt)) %>%
      summarise(area=trapz(rt, int), .by = filename) %>%
      mutate(iso_name=iso_name)
  },.progress = TRUE) %>%
  bind_rows()


# iso_mzs %>%
#   filter() %>%
#   left_join(compound_bounds)
# msdata$MS1[mz%between%pmppm(289.1229, ppm = 10) | mz%between%pmppm(286.1192, ppm = 10)] %>%
#   filter(str_detect(filename, "T0")) %>%
#   mutate(mz=round(mz)) %>%
#   filter(rt%between%c(3.5, 4)) %>%
#   qplotMS1data() +
#   facet_wrap(~mz, ncol=1)

all_iso_areas %>%
  filter(str_detect(iso_name, "Arginine")) %>%
  ggplot() +
  geom_col(aes(x=filename, y=area, fill=iso_name)) +
  facet_wrap(~amendment, ncol=1, scales="free_x")
all_iso_areas %>%
  filter(str_detect(iso_name, "Arginine")) %>%
  left_join(metathaa) %>%
  drop_na() %>%
  ggplot() +
  geom_col(aes(x=tripl, y=area, fill=iso_name), color="black") +
  facet_nested(startime+depth~amendment+timepoint)
all_iso_areas %>%
  filter(str_detect(iso_name, "Glutamic acid")) %>%
  left_join(metathaa) %>%
  drop_na() %>%
  ggplot() +
  geom_col(aes(x=tripl, y=area, fill=iso_name), color="black") +
  facet_nested(startime+depth~amendment+timepoint)

all_iso_areas %>%
  filter(str_detect(iso_name, "2H")) %>%
  left_join(metathaa) %>%
  arrange(timestamp) %>%
  mutate(filename=fct_inorder(filename)) %>%
  ggplot() +
  geom_col(aes(x=filename, y=area, fill=iso_name), color="black") +
  facet_wrap(~iso_name, scales="free_y", ncol=1)



all_iso_areas %>%
  filter(str_detect(iso_name, "13C0, 15N0")) %>%
  mutate(iso_name=str_remove(iso_name, ", 13C0, 15N0")) %>%
  left_join(metathaa) %>%
  filter(!(amendment=="Amm" & startime=="Morn" & depth=="Deep" & tripl=="B")) %>%
  drop_na() %>%
  ggplot() +
  geom_col(aes(x=tripl, y=area, fill=iso_name), color="black", position="fill") +
  facet_nested(startime+depth~amendment+timepoint)




# Quant
all_iso_areas %>%
  filter(iso_name%in%c("Citrulline, 2H4, 13C0, 15N0", "Ornithine, 2H6, 13C0, 15N0", "Taurine, 2H4, 13C0, 15N0")) %>%
  complete(filename, iso_name) %>%
  full_join(metathaa %>% distinct(filename, timestamp)) %>%
  arrange(timestamp) %>%
  mutate(filename=fct_inorder(filename)) %>%
  # ggplot() + geom_col(aes(x=filename, y=area)) + facet_wrap(~iso_name, ncol=1, scales="free_y")
  mutate(area=ifelse(is.na(area), (lag(area)+lead(area))/2, area)) %>%
  pivot_wider(names_from = iso_name, values_from=area) %>%
  # with(cor(`Citrulline, 2H4, 13C0, 15N0`, `Ornithine, 2H6, 13C0, 15N0`))
  with(cor(`Citrulline, 2H4, 13C0, 15N0`, `Taurine, 2H4, 13C0, 15N0`))

is_areas <- all_iso_areas %>%
  filter(iso_name%in%c("Citrulline, 2H4, 13C0, 15N0", 
                       # "Ornithine, 2H6, 13C0, 15N0", 
                       "Taurine, 2H4, 13C0, 15N0")) %>%
  mutate(norm_factor=area/mean(area), .by = iso_name) %>%
  select(-area) %>%
  pivot_wider(names_from = iso_name, values_from=norm_factor) %>%
  mutate(norm_factor=mean(c(`Citrulline, 2H4, 13C0, 15N0`, 
                            # `Ornithine, 2H6, 13C0, 15N0`,
                            `Taurine, 2H4, 13C0, 15N0`), na.rm=TRUE), 
         .by = filename) %>%
  select(filename, norm_factor) %>%
  arrange(desc(norm_factor))
all_iso_areas %>%
  left_join(is_areas) %>%
  mutate(norm_area=area/norm_factor, .by = iso_name) %>%
  filter(iso_name=="Arginine, 13C0, 15N0") %>%
  full_join(metathaa %>% distinct(filename, timestamp)) %>%
  arrange(timestamp) %>%
  mutate(filename=fct_inorder(filename)) %>%
  ggplot() +
  geom_col(aes(x=filename, y=norm_area))

stan_curve_data <- all_iso_areas %>%
  filter(str_detect(iso_name, "13C0, 15N0")) %>%
  left_join(metathaa %>% distinct(filename, samp_type)) %>%
  filter(samp_type=="Std") %>%
  mutate(conc=str_extract(filename, "(?<=Std_).*(?=uM)")) %>%
  mutate(conc=str_replace(conc, "-", ".")) %>%
  mutate(conc=as.numeric(conc)) %>%
  mutate(rep=str_extract(filename, "\\d(?=.mzML)")) %>%
  left_join(is_areas) %>%
  mutate(area=area/norm_factor)

stan_curve_data %>%
  ggplot(aes(x=conc, y=area)) +
  geom_point(aes(color=rep)) +
  geom_smooth(method=lm, formula=y~x) +
  geom_hline(yintercept = 0) +
  facet_wrap(~iso_name, scales="free_y", ncol=3)

cal_curves <- stan_curve_data %>%
  filter(!str_detect(iso_name, "2H")) %>%
  select(filename, iso_name, area, conc) %>%
  nest(data=c(filename, area, conc)) %>%
  mutate(lm_fit=map(data, ~broom::tidy(lm(.x$area~.x$conc)))) %>%
  unnest(lm_fit) %>%
  mutate(compound_name=str_remove(iso_name, ", .*")) %>%
  select(compound_name, term, estimate) %>%
  filter(term!="(Intercept)") %>%
  select(compound_name, slope=estimate)


quant_data <- all_iso_areas %>%
  left_join(is_areas) %>%
  mutate(norm_area=area/norm_factor, .by = iso_name) %>%
  select(filename, iso_name, norm_area) %>%
  mutate(compound_name=str_remove(iso_name, ", .*")) %>%
  left_join(cal_curves) %>%
  mutate(conc_um_in_vial=norm_area/slope) %>%
  mutate(env_conc_nm=conc_um_in_vial/100*1000) %>%
  filter(iso_name!="Proline, 13C2, 15N1") %>%
  filter(iso_name!="Proline, 13C3, 15N1") %>%
  filter(iso_name!="Proline, 13C4, 15N1") %>%
  filter(compound_name!="Lysine")
write_csv(quant_data, "thaas/quant_data.csv")



# Plots!
quant_data %>%
  filter(str_detect(iso_name, "13C0, 15N0")) %>%
  mutate(iso_name=str_remove(iso_name, ", 13C0, 15N0")) %>%
  left_join(metathaa) %>%
  filter(!(amendment=="Amm" & startime=="Morn" & depth=="Deep" & tripl=="B")) %>%
  drop_na() %>%
  ggplot() +
  geom_col(aes(x=tripl, y=env_conc_nm, fill=iso_name), color="black") +
  facet_nested(startime+depth~amendment+timepoint)
quant_data %>%
  filter(str_detect(iso_name, "Glutamic acid")) %>%
  left_join(metathaa) %>%
  drop_na() %>%
  ggplot() +
  geom_col(aes(x=tripl, y=env_conc_nm, fill=iso_name), 
           position="fill",
           color="black", width=1) +
  facet_nested(startime+depth~amendment+timepoint)
quant_data %>%
  filter(str_detect(iso_name, "Arginine")) %>%
  left_join(metathaa) %>%
  drop_na() %>%
  ggplot() +
  geom_col(aes(x=tripl, y=env_conc_nm, fill=iso_name), 
           position="fill",
           color="black", width=1) +
  facet_nested(startime+depth~amendment+timepoint)

pdf("thaas/all_cmpd_labels.pdf", width=9, height = 5)
for(cmpd_i in setdiff(unique(quant_data$compound_name), c("Citrulline", "Ornithine", "Taurine"))){
  print(cmpd_i)
  gp <- quant_data %>%
    filter(compound_name==cmpd_i) %>%
    left_join(metathaa) %>%
    drop_na() %>%
    ggplot() +
    geom_col(aes(x=tripl, y=env_conc_nm, fill=iso_name), 
             position="fill",
             color="black", width=1) +
    facet_nested(depth+startime~amendment+timepoint) +
    ggtitle(cmpd_i)
  print(gp)
}
dev.off()




quant_data %>%
  # filter(compound_name=="Histidine") %>%
  filter(compound_name=="Tyrosine") %>%
  left_join(metathaa) %>%
  drop_na() %>%
  ggplot(aes(x=tripl, y=env_conc_nm, fill=iso_name, label=iso_name)) +
  geom_col(position="fill", color="black", width=1) +
  facet_nested(depth+startime~amendment+timepoint)
plotly::ggplotly()
# Tyrosine, 13C4, 15N1 from GMP?
# Phenylalanine, 13C4, 15N1 from GMP?
# Histidine, 13C5, 15N3 from GMP?



quant_data %>%
  filter(compound_name%in%c("Phenylalanine", "Tyrosine")) %>%
  mutate(iso_type=str_extract(iso_name, "13C.*")) %>%
  left_join(metathaa) %>%
  drop_na() %>%
  filter(amendment=="GMP") %>%
  ggplot(aes(x=tripl, y=env_conc_nm, fill=iso_type, label=iso_name)) +
  geom_col(position="fill", color="black", width=1) +
  facet_nested(depth+startime~compound_name+timepoint)
# 4 carbons from PPP (E4P), 15N from highly labeled pool of glutamate?

iso_mzs %>% filter(iso_name=="Phenylalanine, 13C4, 15N1")
msdata$MS1[mz%between%pmppm(341.1453, ppm = 10) | mz%between%pmppm(336.1348, ppm = 10)] %>%
  filter(str_detect(filename, "GMP")) %>%
  mutate(mz=round(mz)) %>%
  filter(rt%between%c(5.2, 5.5)) %>%
  qplotMS1data() +
  facet_wrap(~mz, ncol=1, scales="free_y")
iso_mzs %>% filter(iso_name=="Tyrosine, 13C4, 15N1")
msdata$MS1[mz%between%pmppm(357.1402, ppm = 10) | mz%between%pmppm(352.1297, ppm = 10)] %>%
  filter(str_detect(filename, "GMP")) %>%
  mutate(mz=round(mz)) %>%
  filter(rt%between%c(4, 4.5)) %>%
  qplotMS1data() +
  facet_wrap(~mz, ncol=1, scales="free_y")




iso_mzs %>% filter(iso_name=="Histidine, 13C5, 15N3")
msdata$MS1[mz%between%pmppm(334.1332, ppm = 10) | mz%between%pmppm(326.1253, ppm = 10)] %>%
  filter(str_detect(filename, "GMP")) %>%
  mutate(mz=round(mz)) %>%
  filter(rt%between%c(1, 2)) %>%
  qplotMS1data() +
  facet_wrap(~mz, ncol=1, scales="free_y")
