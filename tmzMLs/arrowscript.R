library(data.table)
library(RaMS)
library(arrow)

parametadata <- fread("metadata/parametadata.csv")[
  IS_type=="noIS"|is.na(IS_type),c("filepath", "filename", "samp_type", "amendment", "polarity")
]
parametadata[, .N, by=c("samp_type", "amendment", "polarity")]

for(samp_type_i in unique(parametadata$samp_type)){
  for(amendment_i in unique(parametadata$amendment)){
    for(polarity_i in unique(parametadata$polarity)){
      ms_files <- parametadata[samp_type==samp_type_i][amendment==amendment_i][polarity==polarity_i]$filepath
      # ms_files <- head(ms_files, 3)
      if(length(ms_files)==0)next
      msdata <- grabMSdata(ms_files, grab_what = "MS1")
      parq_dir <- paste0("tmzMLs/pqds",
                         "/samp_type=", samp_type_i, 
                         "/amendment=", amendment_i, 
                         "/polarity=", polarity_i)
      if(!dir.exists(parq_dir))dir.create(parq_dir, recursive = TRUE)
      parq_name <- paste0(parq_dir, "/part-0.parquet")
      print(parq_name)
      write_parquet(msdata$MS1[order(mz)], parq_name)
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




# library(DBI)
# ms_files <- parametadata[polarity=="pos"]$filepath
# ms_file_list <- split(ms_files, 1:10)
# msduck <- dbConnect(duckdb::duckdb(), "tmzMLs/duckdb")
# dbWriteTable(msduck, "MS1_pos", grabMSdata(ms_file_list[[1]], grab_what = "MS1")$MS1[order(mz)])
# pbapply::pblapply(ms_file_list[-1], function(files_i){
#   ms1_data <- grabMSdata(files_i, grab_what = "MS1", verbosity = 0)$MS1[order(mz)]
#   dbAppendTable(msduck, name = "MS1_pos", value = ms1_data)
# })
# dbDisconnect(msduck)

# library(DBI)
# msduck <- dbConnect(duckdb::duckdb(), "tmzMLs/duckdb")
# lower_bound <- 118.0864
# upper_bound <- 118.0866
# # lower_bound <- 90.0552
# # upper_bound <- 90.0558
# dbq <- paste("SELECT * FROM MS1_pos WHERE mz BETWEEN", lower_bound, "AND", upper_bound, 
#              "AND rt BETWEEN 10 AND 12")
# qplotMS1data(dbGetQuery(msduck, dbq))
