library(tidyverse)
library(arrow)
library(RaMS)

IS_type <- "noIS"

all_stans <- read_csv("metadata/all_stans.csv") %>%
  mutate(polarity=factor(polarity, levels=c("pos", "neg"))) %>%
  arrange(polarity) %>%
  filter(compound_type!="Internal Standard")

for(polarity in c("pos", "neg")){
  for(amendment in c("Amm", "Nit", "GMP", "Arg")){
    parametadata <- read_csv("metadata/parametadata.csv") %>% 
      filter(polarity==!!polarity) %>%
      filter(amendment==!!amendment) %>%
      filter(IS_type==!!IS_type)
    output_folder <- paste(amendment, polarity, IS_type, "output", sep = "_") %>%
      paste0("untargeted/", ., "/")
    best_feat_data <- read_csv(paste0(output_folder, "peak_data.csv")) %>%
      left_join(parametadata) %>%
      filter(samp_type=="Smp") %>%
      group_by(feature) %>%
      summarize(mzmed=unique(feat_mzmed), rtmed=unique(feat_rtmed)/60, n=n(),
                med_cor=median(beta_cor, na.rm=TRUE), 
                med_snr=median(beta_snr, na.rm=TRUE)) %>%
      filter(med_cor>0.75) %>%
      select(-med_cor, -med_snr)
    
    pdf(paste0(output_folder, "anno_eics.pdf"), width = 8, height = 4)
    all_stans %>%
      filter(polarity==!!polarity) %>%
      mutate(mz=ifelse(compound_name=="Carnitine", mz[which(compound_name=="Threonine betaine?")], mz)) %>%
      group_by(mz) %>%
      summarise(cmpd_lab=paste0(compound_name, " (", mix, " ,", rt, ")", collapse = ", ")) %>%
      pwalk(function(...){
        row_data <- list(...)
        mass_feats <- best_feat_data %>% filter(mzmed%between%pmppm(row_data$mz, 5))
        chosen_eic <- open_dataset("tmzMLs/pqds") %>%
          filter(polarity==!!polarity) %>%
          filter(amendment==!!amendment) %>%
          filter(samp_type=="Poo" | samp_type=="Std") %>%
          filter(mz%between%pmppm(row_data$mz, 5)) %>%
          filter(rt%between%c(2, 20)) %>%
          dplyr::collect()
        gp_plot <- chosen_eic %>%
          left_join(parametadata, by = join_by(filename, samp_type, amendment, polarity)) %>%
          mutate(samp_type=str_extract(filename, "Mix\\d|Deep|Surf|Poo|Blk")) %>%
          filter(!str_detect(filename, "Mat")) %>%
          mutate(plot_facet=ifelse(samp_type%in%c("Mix1", "Mix2"), "Mixes", "Samples")) %>%
          ggplot() +
          geom_line(aes(x=rt, y=int, group=filename, color=samp_type)) +
          facet_wrap(~plot_facet, scales = "free_y", ncol=1) +
          geom_vline(aes(xintercept=rtmed), color="black", data = mass_feats) +
          geom_text(aes(x=rtmed, y=Inf, label=feature), 
                    angle=90, hjust=1, vjust=1, data = mass_feats) +
          ggtitle(row_data$cmpd_lab) +
          scale_x_continuous(breaks = 2:20) +
          theme_bw() +
          theme(legend.position = c(0.9, 0.8), legend.title = element_blank(),
                axis.title = element_blank())
        print(gp_plot)
      }, .progress = TRUE)
    dev.off()
  }
}
