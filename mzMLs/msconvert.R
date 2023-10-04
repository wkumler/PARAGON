
folder_paths <- c(
  "Z:/1_QEdata/2021/2021_Will_PARAGON/210813_Will_NIT-Expt_PARAGON_HILIC",
  "Z:/1_QEdata/2021/2021_Will_PARAGON/210916_Will_Amm-Expt_PARAGON_HILIC",
  "Z:/1_QEdata/2022/220715_Will_GMP-Expt_PARAGON_HILIC",
  "Z:/1_QEdata/2022/220729_Will_Arg-Expt_PARAGON_HILIC"
)

sapply(folder_paths, function(path_i){
  sapply(c("positive", "negative"), function(polarity){
    message(paste("Converting", basename(path_i), polarity))
    short_pol <- gsub("[ia]tive", "", polarity)
    match_pattern <- paste0("wIS.raw|noIS-", short_pol, ".*.raw|_Std_.*raw|wIS-(Full|Half)")
    files_to_convert <- list.files(path_i, pattern = match_pattern, full.names = TRUE)
    mscommand <- paste(
      "msconvert",
      paste(files_to_convert, collapse = " "),
      "--mzML",
      paste0('-o ', short_pol),
      '--filter "peakPicking true 1-"',
      paste('--filter "polarity', polarity, '"')
    )
    system(mscommand)
    
    msmsmatch_pattern <- paste0("_DDA", short_pol, ".*.raw")
    msms_to_convert <- list.files(path_i, pattern = msmsmatch_pattern, full.names = TRUE)
    msmscommand <- paste(
      "msconvert",
      paste(msms_to_convert, collapse = " "),
      "--mzML",
      paste0('-o ', short_pol, "/MSMS"),
      '--filter "peakPicking true 1-"',
      paste('--filter "polarity', polarity, '"')
    )
    system(msmscommand)
  })
})
