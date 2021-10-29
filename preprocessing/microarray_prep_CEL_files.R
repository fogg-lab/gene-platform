library(tidyverse)

# Custom package
library(rutils)

# Helper functions
untar_gse <- function(tarfile) {
    extract_dest <- paste0(dirs$data_dir, "/raw_GEO_data/untarred/", str_extract(tarfile, "GSE[0-9]+"))
    untar(tarfile, exdir = extract_dest)
}

# Main script
dirs <- rutils::get_dev_directories("../dev_paths.txt")
tarfiles <- list.files(paste0(dirs$data_dir, "/raw_GEO_data/tarred"), pattern = "*.tar", full.names = TRUE)

length(tarfiles)

plyr::laply(tarfiles, untar_gse)

