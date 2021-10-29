library(tidyverse)
if (!("devtools" %in% (installed.packages() %>% as_tibble() %>% pull(Package)))) {
    install.packages("devtools")
}

library(devtools)

tryCatch(
    expr = {
        devtools::uninstall()
    },
    finally = {
        devtools::document()
        devtools::load_all()
        devtools::install()
    }
)
