cur_dir <- getwd()
setwd("/home/chris/starr_lab/ocpmi/notebook")

rmarkdown::render("index.Rmd")
rmarkdown::render("notes.Rmd")
rmarkdown::render("results.Rmd")
rmarkdown::render("methods.Rmd")
rmarkdown::render("references.Rmd")
rmarkdown::render_site("index.Rmd")

fl <- Sys.glob("*.html")

ifelse(file.exists(fl), file.remove(fl))

setwd(cur_dir)
