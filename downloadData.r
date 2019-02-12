#' 
#' Downloading all data sets rquired for running the pipeline



url ="http://amp.pharm.mssm.edu/static/hdfs/harmonizome/data/lincscmapchemical/gene_attribute_edges.txt."

setwd("/home/memon/projects/msdrp/")


if(!http_error(url) == TRUE){
  download.file(url,destfile = './data/test.txt.gz')
} else {
  print("The url is outdated, please update!")
}

downloadData <- function(url, folder) {
  filename <- basename(url)
  base <- tools::file_path_sans_ext(filename)
  ext <- tools::file_ext(filename)
  
  file_exists <- grepl(base, list.files(folder), fixed = TRUE)
  
  if (any(file_exists))
  {
    filename <- paste0(base, " (", sum(file_exists), ")", ".", ext)
  }
  
  download.file(url, file.path(folder, filename), mode = "wb", method = "libcurl")
}

download_without_overwrite(
  url = "https://raw.githubusercontent.com/nutterb/redcapAPI/master/README.md",
  folder = "[path_to_folder]")