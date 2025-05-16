
library(tidyverse)
# ms_files <- c(
#   "Z:/1_QEdata/2021/2021_Will_PARAGON/210813_Will_NIT-Expt_PARAGON_HILIC/positive",
#   "Z:/1_QEdata/2021/2021_Will_PARAGON/210916_Will_Amm-Expt_PARAGON_HILIC/positive",
#   "Z:/1_QEdata/2022/220715_Will_GMP-Expt_PARAGON_HILIC/positive",
#   "Z:/1_QEdata/2022/220729_Will_Arg-Expt_PARAGON_HILIC/positive"
# ) %>%
#   lapply(list.files, full.names=TRUE) %>%
#   unlist() %>%
#   str_subset("noIS") %>%
#   str_subset("-neg", negate = TRUE) %>%
#   # str_subset("Poo.*Full") %>%
#   str_subset("DDA", negate = TRUE)
# file.copy(ms_files, paste0("mzMLs/"))
ms_files <- list.files("mzMLs", full.names = TRUE)


# library(RaMS)
# msdata <- grabMSdata(ms_files %>% str_subset("Poo.*Full"))
# msdata$TIC %>%
#   mutate(amendment=str_extract(filename, "Amm|Arg|GMP|Nit")) %>%
#   qplotMS1data(color_col="filename") +
#   facet_wrap(~amendment, ncol=1)
# msdata$MS1[mz%between%pmppm(118.0862, 10)] %>%
#   slice_max(int, by=c(filename, rt)) %>%
#   mutate(amendment=str_extract(filename, "Amm|Arg|GMP|Nit")) %>%
#   qplotMS1data(color_col = "amendment")


library(xcms)
register(BPPARAM = SnowParam(progressbar = TRUE, workers = 10, tasks = length(ms_files)))
msnexp <- readMSData(ms_files, msLevel. = 1, mode = "onDisk")

cwp <- CentWaveParam(
  ppm=5,
  peakwidth = c(20, 80),
  prefilter=c(5, 1e4),
  snthresh=5,
  verboseColumns=TRUE,
  extendLengthMSW=TRUE,
  integrate=2
)
msnexp_withpeaks <- findChromPeaks(msnexp, cwp)
saveRDS(msnexp_withpeaks, "msnexp_withpeaks.rds")



msnexp_withpeaks <- readRDS("msnexp_withpeaks.rds")
register(BPPARAM = SerialParam(progressbar = TRUE))
pdp <- PeakDensityParam(
  sampleGroups = str_remove(ms_files, "-[A-C](?=-)|_\\d(?=\\.mzML)"), 
  bw=12,
  minFraction = 1,
  binSize = 0.001
)
msnexp_grouped <- groupChromPeaks(msnexp_withpeaks, pdp)

register(BPPARAM = SnowParam(progressbar = TRUE, workers = 5, tasks = length(ms_files)))
msnexp_filled <- fillChromPeaks(msnexp_grouped, FillChromPeaksParam())
saveRDS(msnexp_filled, "msnexp_filled.rds")




msnexp_filled <- readRDS("msnexp_filled.rds")
peak_df <- msnexp_filled %>%
  chromPeaks() %>%
  as.data.frame() %>%
  mutate(filename=basename(ms_files)[sample]) %>%
  select(mz, mzmin, mzmax, rt, rtmin, rtmax, into, filename) %>%
  rownames_to_column("peak_id")
feature_df <- featureDefinitions(msnexp_filled) %>%
  as.data.frame() %>%
  select(mzmed, rtmed, npeaks, peakidx) %>%
  rownames_to_column("feature")
oaoo_df <- feature_df %>%
  # filter(feature=="FT0001") %>%
  unnest(peakidx) %>%
  left_join(peak_df %>% mutate(peakidx=row_number()), by = join_by(peakidx)) %>%
  arrange(desc(into)) %>%
  slice(1, .by = c(feature, filename)) %>%
  arrange(feature, filename) %>%
  complete(feature, filename) %>%
  # filter(feature=="FT0003") %>%
  mutate(amendment=str_extract(filename, "Amm|GMP|Arg|Nit")) %>%
  arrange(amendment, feature, desc(into)) %>%
  fill(mzmed, rtmed, into)
oaoo_df %>%
  count(feature) %>%
  distinct(n)
saveRDS(oaoo_df, "one_and_only_one_peak_per_feature_file.rds")



# oaoo_df <- readRDS("one_and_only_one_peak_per_feature_file.rds")
stats_calc <- oaoo_df %>%
  mutate(timepoint=as.numeric(str_extract(filename, "(?<=T)\\d+"))) %>%
  mutate(timepoint=ifelse(timepoint==0, 0.3, timepoint)) %>%
  drop_na(timepoint) %>%
  select(feature, mzmed, rtmed, filename, into, amendment, timepoint) %>%
  nest(data=c(filename, into, timepoint)) %>%
  mutate(lm_test=map(data, ~broom::tidy(lm(log10(.x$into)~log10(.x$timepoint))), .progress = TRUE)) %>%
  unnest(lm_test) %>%
  filter(term=="log10(.x$timepoint)") %>%
  mutate(p_adj=p.adjust(p.value, method = "fdr"))
saveRDS(stats_calc, file = "lmtest_stats_calc.rds")



# msnexp_filled <- readRDS("msnexp_filled.rds")
# peak_df <- msnexp_filled %>%
#   chromPeaks() %>%
#   as.data.frame() %>%
#   mutate(filename=basename(ms_files)[sample]) %>%
#   select(mz, mzmin, mzmax, rt, rtmin, rtmax, into, filename) %>%
#   rownames_to_column("peak_id")
peak_quality_scores <- peak_df %>%
  split(.$filename) %>%
  # head(1) -> peak_df_i
  # head(5) %>%
  BiocParallel::bplapply(function(peak_df_i){
    library(tidyverse)
    library(RaMS)
    source("squallms_functions.R")
    msdata_i <- grabMSdata(paste0("mzMLs/", peak_df_i$filename[1]), verbosity = 0)
    peak_qscores_i <- peak_df_i %>%
      filter(mz%between%pmppm(118.0865+0.997035, 10)) %>%
      # slice(1) -> row_data
      # head(100) %>%
      pmap(function(...){
        row_data <- list(...)
        rt_bounded_eic <- msdata_i$MS1[rt%between%c(row_data$rtmin/60, row_data$rtmax/60)]
        unlab_eic <- rt_bounded_eic[mz%between%pmppm(row_data$mz-0.997035, 10)] %>%
          .[,mz:=NULL] %>%
          .[,filename:=NULL] %>%
          .[order(-int),.SD[1], by=rt] %>%
          .[order(rt)]
        peak_qscore <- qscoreCalculator(unlab_eic$rt, unlab_eic$int)
        lab_eic <- rt_bounded_eic[mz%between%pmppm(row_data$mz, 10)] %>%
          .[,mz:=NULL] %>%
          .[,filename:=NULL] %>%
          .[order(-int),.SD[1], by=rt] %>%
          .[order(rt)]
        merged_eic <- merge(unlab_eic, lab_eic, by="rt")
        if(nrow(merged_eic)>5){
          iso_cor <- cor(merged_eic$int.x, merged_eic$int.y)
          unlab_area <- trapz(unlab_eic$rt, unlab_eic$int)
          lab_area <- trapz(lab_eic$rt, lab_eic$int)
        } else {
          iso_cor <- NA
          unlab_area <- NA
          lab_area <- NA
        }
        list(peak_id=row_data$peak_id, beta_snr=peak_qscore[1], beta_cor=peak_qscore[2], 
             iso_cor=iso_cor, unlab_area=unlab_area, lab_area=lab_area)
      }) %>%
      bind_rows()
  }, BPPARAM = SnowParam(progressbar = TRUE, workers = 10, tasks = length(ms_files))) %>%
  bind_rows() %>%
  arrange(peak_id)
saveRDS(peak_quality_scores, "peak_quality_scores.rds")















feature_ids <- featureDefinitions(msnexp_filled) %>%
  as.data.frame() %>%
  rownames_to_column("feature") %>% 
  select(feature, mzmed, rtmed, peak_id=peakidx) %>% 
  unnest(peak_id)
peak_ids <- msnexp_filled %>%
  chromPeaks() %>%
  as.data.frame() %>%
  mutate(filename=basename(ms_files)[sample]) %>%
  rownames_to_column("peak_id") %>%
  mutate(peak_id=as.numeric(str_remove(peak_id, "^CP"))) %>%
  arrange(peak_id) %>%
  select(peak_id, mz, rt, into, filename)
cleaned_combined_peaks <- peak_quality_scores %>%
  drop_na() %>%
  mutate(peak_id=as.numeric(str_remove(peak_id, "^CP"))) %>%
  arrange(peak_id) %>%
  left_join(feature_ids) %>%
  filter(!is.na(feature)) %>%
  left_join(peak_ids) %>%
  select(feature, mzmed, rtmed, filename, peak_id, mz, rt, into, everything())




lmtest_output <- cleaned_combined_peaks %>%
  arrange(feature) %>%
  select(feature, filename, into) %>%
  mutate(samp_type=str_extract(filename, "Smp|Blk|Poo")) %>%
  filter(samp_type=="Smp") %>%
  mutate(timepoint=as.numeric(str_extract(filename, "(?<=T)\\d+"))) %>%
  mutate(timepoint=ifelse(timepoint==0, 0.3, timepoint)) %>%
  mutate(amendment=str_extract(filename, "Amm|GMP|Arg|Nit")) %>%
  mutate(depth=str_extract(filename, "Surf|Deep")) %>%
  nest(data=-feature) %>%
  filter(sapply(data, nrow)>2) %>%
  mutate(lm_test=map(data, ~broom::tidy(lm(log10(.x$into)~log10(.x$timepoint))), .progress = TRUE)) %>%
  unnest(lm_test) %>%
  filter(term=="log10(.x$timepoint)") %>%
  mutate(p_adj=p.adjust(p.value, method = "fdr"))


iso_df <- cleaned_combined_peaks %>%
  summarise(mzmed=unique(mzmed),
            rtmed=unique(rtmed),
            med_iso_cor=median(iso_cor), 
            med_snr=median(beta_snr),
            med_cor=median(beta_cor),
            .by=feature) %>%
  left_join(lmtest_output %>% select(feature, estimate, p_adj)) %>%
  mutate(peak_qual=ifelse(med_snr>-50*med_cor+50, "Good", "Bad")) %>%
  mutate(is_iso=ifelse(med_iso_cor>0.8, "Is iso", "Not iso")) %>%
  mutate(time_trend=ifelse(p_adj<0.001& estimate>0, "Increases", "No change"))
iso_df %>%
  filter(!is.na(time_trend)) %>%
  count(peak_qual, is_iso, time_trend)
best_labeled <- iso_df %>%
  filter(peak_qual=="Good" & is_iso=="Is iso" & time_trend=="Increases") %>%
  arrange(feature) %>%
  select(feature, mzmed, rtmed, med_iso_cor, med_snr, med_cor, p_adj)


# Annotate best ones quickly with standards list, to be finessed later manually
chosen_stans <- read_csv("https://github.com/IngallsLabUW/Ingalls_Standards/raw/refs/heads/master/Ingalls_Lab_Standards.csv") %>%
  filter(z>0) %>%
  filter(Column=="HILIC") %>%
  summarise(compound_name=paste(Compound_Name, collapse = "/"), .by = mz) %>%
  rowwise() %>%
  mutate(N15_mz_lower=pmppm(mz+0.997035, 10)[1]) %>%
  mutate(N15_mz_upper=pmppm(mz+0.997035, 10)[2]) %>%
  select(compound_name, N15_mz_lower, N15_mz_upper)
best_labeled %>%
  left_join(chosen_stans, join_by(between(x$mzmed, y$N15_mz_lower, y$N15_mz_upper))) %>%
  select(-N15_mz_lower, -N15_mz_upper) %>%
  mutate(compound_name=ifelse(duplicated(compound_name)&!is.na(compound_name), "Isocytosine", compound_name)) %>%
  mutate(compound_name=ifelse((mzmed-lag(mzmed))%between%c(0.996, 0.998), "iso_prev?", compound_name)) %>%
  mutate(rtmed=rtmed/60) %>%
  select(feature, mzmed, rtmed, compound_name) %>%
  write_csv("best_labeled_anno.csv")


# Repeat analysis across amendments to get additional numbers
lmtest_xamend <- cleaned_combined_peaks %>%
  arrange(feature) %>%
  select(feature, filename, into) %>%
  mutate(samp_type=str_extract(filename, "Smp|Blk|Poo")) %>%
  filter(samp_type=="Smp") %>%
  mutate(timepoint=as.numeric(str_extract(filename, "(?<=T)\\d+"))) %>%
  mutate(timepoint=ifelse(timepoint==0, 0.3, timepoint)) %>%
  mutate(amendment=str_extract(filename, "Amm|GMP|Arg|Nit")) %>%
  mutate(depth=str_extract(filename, "Surf|Deep")) %>%
  nest(data=-c(feature, amendment)) %>%
  filter(sapply(data, nrow)>2) %>%
  mutate(lm_test=map(data, ~broom::tidy(lm(log10(.x$into)~log10(.x$timepoint))), .progress = TRUE)) %>%
  unnest(lm_test) %>%
  filter(term=="log10(.x$timepoint)") %>%
  mutate(p_adj=p.adjust(p.value, method = "fdr"))
iso_df_xamend <- cleaned_combined_peaks %>%
  mutate(amendment=str_extract(filename, "Amm|GMP|Arg|Nit")) %>%
  summarise(mzmed=unique(mzmed),
            rtmed=unique(rtmed),
            med_iso_cor=median(iso_cor), 
            med_snr=median(beta_snr),
            med_cor=median(beta_cor),
            .by=c(feature, amendment)) %>%
  left_join(lmtest_xamend %>% select(feature, amendment, estimate, p_adj)) %>%
  mutate(peak_qual=ifelse(med_snr>-50*med_cor+50, "Good", "Bad")) %>%
  mutate(is_iso=ifelse(med_iso_cor>0.8, "Is iso", "Not iso")) %>%
  mutate(time_trend=ifelse(p_adj<0.001& estimate>0, "Increases", "No change"))

iso_df_xamend %>%
  bind_rows(iso_df %>% mutate(amendment="Overall")) %>%
  filter(!is.na(time_trend)) %>%
  count(peak_qual, is_iso, time_trend, amendment) %>%
  pivot_wider(names_from = amendment, values_from=n, values_fill = 0) %>%
  print()
iso_df_xamend %>%
  bind_rows(iso_df %>% mutate(amendment="Overall")) %>%
  filter(!is.na(time_trend)) %>%
  ggplot() +
  geom_bar(aes(x=amendment)) +
  facet_wrap(~peak_qual+is_iso+time_trend, scales="free_y", ncol=4)



labeled_fracs <- cleaned_combined_peaks %>%
  summarise(med_iso_cor=median(iso_cor), 
            med_snr=median(beta_snr),
            med_cor=median(beta_cor),
            .by=c(feature)) %>%
  filter(med_snr>-50*med_cor+50) %>%
  filter(med_iso_cor>0.8) %>%
  select(feature) %>%
  left_join(cleaned_combined_peaks %>% select(feature, filename, lab_area, unlab_area), multiple = "all") %>%
  mutate(labeled_frac=lab_area/(lab_area+unlab_area)) %>%
  select(-lab_area, -unlab_area)

cmpd_order <- labeled_fracs %>%
  mutate(samp_type=str_extract(filename, "Smp|Blk|Poo")) %>%
  filter(samp_type=="Smp") %>%
  mutate(timepoint=str_extract(filename, "T\\d+")) %>%
  mutate(amendment=str_extract(filename, "Amm|GMP|Arg|Nit")) %>%
  mutate(depth=str_extract(filename, "Surf|Deep")) %>%
  summarise(mean_label=mean(labeled_frac, na.rm=TRUE), .by=feature) %>%
  arrange(desc(mean_label)) %>%
  pull(feature)

labeled_fracs %>%
  mutate(samp_type=str_extract(filename, "Smp|Blk|Poo")) %>%
  filter(samp_type=="Smp") %>%
  mutate(timepoint=str_extract(filename, "T\\d+")) %>%
  mutate(amendment=str_extract(filename, "Amm|GMP|Arg|Nit")) %>%
  mutate(depth=str_extract(filename, "Surf|Deep")) %>%
  mutate(startime=str_extract(filename, "Morn|Eve")) %>%
  mutate(tripl=str_extract(filename, "(?<=-)[A-C](?=-)")) %>%
  mutate(feature=factor(feature, levels=cmpd_order)) %>%
  ggplot() +
  geom_raster(aes(x=timepoint, y=tripl, fill=labeled_frac)) +
  scale_fill_gradientn(colours = viridis::viridis(10, alpha=seq(0, 1, length.out=10), direction = -1),
                       limits = c(0, 1)) +
  facet_wrap(feature~amendment+depth+startime)
