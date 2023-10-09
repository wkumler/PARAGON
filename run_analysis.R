
library(tidyverse)

# Generate metadata ----
if(!file.exists("metadata/parametadata.csv")){
  knitr::knit(input = "metadata/parametadata_control.Rmd", 
              output = "metadata/parametadata_control.R", 
              tangle = TRUE)
}


# Generate XCMS outputs ----
if(!file.exists("untargeted/Arg_neg_noIS_output/best_feat_data.csv")){
  for(amendment in c("Amm", "Nit", "GMP", "Arg")){
    for(polarity in c("pos", "neg")){
      source("untargeted/paruntargeted_control.R")
    }
  }
}


# Generate parquet files for fast EIC extraction ----
if(!dir.exists("tmzMLs/pqds")){
  library(data.table)
  library(RaMS)
  library(arrow)
  parametadata <- fread("metadata/parametadata.csv")[
    IS_type=="noIS",c("filepath", "filename", "samp_type", "amendment", "polarity")
  ]
  parametadata[, .N, by=c("samp_type", "amendment", "polarity")]
  
  for(samp_type_i in unique(parametadata$samp_type)){
    for(amendment_i in unique(parametadata$amendment)){
      for(polarity_i in unique(parametadata$polarity)){
        ms_files <- parametadata[samp_type==samp_type_i][amendment==amendment_i][polarity==polarity_i]$filepath
        # ms_files <- head(ms_files, 3)
        msdata <- grabMSdata(ms_files, grab_what = "MS1")
        parq_dir <- paste0("tmzMLs/pqds",
                           "/samp_type=", samp_type_i, 
                           "/amendment=", amendment_i, 
                           "/polarity=", polarity_i)
        if(!dir.exists(parq_dir))dir.create(parq_dir, recursive = TRUE)
        parq_name <- paste0(parq_dir, "/part-0.parquet")
        print(parq_name)
        write_parquet(msdata$MS1, parq_name)
      }
    }
  }
}


# Generate annotation chromatogram pdfs ----
if(!file.exists("untargeted/Arg_neg_noIS_output/anno_eics.pdf")){
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
}

