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

# library(tidyverse)
# open_dataset("tmzMLs/pqds") %>%
#   filter(amendment=="Amm") %>%
#   filter(polarity=="neg") %>%
#   filter(samp_type=="Smp") %>%
#   filter(mz%between%pmppm(150.041585, 5)) %>%
#   filter(rt%between%c(7.5, 9)) %>%
#   collect() %>%
#   left_join(parametadata) %>%
#   ggplot() +
#   geom_line(aes(x=rt, y=int, group=filename, color=samp_type))
