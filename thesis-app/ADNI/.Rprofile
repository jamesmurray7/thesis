cat("\f")

library(gmvjoint)
options(
	prompt = glue::glue_col("{green ADNI> }"),
	continue = "  "
)

save.dir <- ifelse(grepl("macOS", utils::sessionInfo()$running),
                   '../macresults.nosync',
                   '/data/c0061461/THESIS/ADNI')
save.dir.file.path <- function(x, LHS = save.dir) file.path(LHS, x)


.startmessage <- function(){
  cat("\f")
  message("ADNI - Thesis Application")
  cat("Current existing directories:\n")
  fs::dir_tree(save.dir)
  cat("\n")
  assign("out", get(load(save.dir.file.path("ADNI.RData"))))
  adni <<- out
}

source('../../thesis-sims/theme_csda.R')

.startmessage()
surv.data <- adni[!duplicated(adni[,'id']), ]
