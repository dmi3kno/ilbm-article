library(prereviewr)
library(rmarkdown)
library(tidyverse)

uuid <- prereviewr::get_preprints(search="the tenets") %>% pull(uuid)

reviews <- prereviewr::get_preprint_reviews(uuid, full=TRUE)
review1 <- reviews$drafts[[1]][[1]][["contents"]]
write_lines(review1, "PREreviews/reviewer1.html")
rmarkdown::pandoc_convert(here::here("PREreviews", "reviewer1.html"), to="markdown_strict",
                          output = here::here("PREreviews", "reviewer1.Rmd"))
