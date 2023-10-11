cat("\f")
library(gmvjoint)
data(PBC, package = "gmvjoint")
PBC$id <- as.numeric(as.character(PBC$id))
options(
	prompt = glue::glue_col("{green PBC> }"),
	continue = "  ",
	stringsAsFactors = FALSE
)

save.dir <- ifelse(grepl("macOS", utils::sessionInfo()$running),
                   '../macresults.nosync',
                   '/data/c0061461/THESIS/PBC')
save.dir.file.path <- function(x, LHS = save.dir) file.path(LHS, x)

# Create some extra datasets
PBCredx <- na.omit(PBC[,c("id", "survtime", "status", "drug", "sex",
                          "age", "histologic", "time", "hepatomegaly", "spiders",
                          "alkaline", "SGOT", "platelets", "prothrombin", 
                          "serBilir", "albumin")])
PBCredx$serBilir <- log(PBCredx$serBilir)
PBCredx$SGOT <- log(PBCredx$SGOT)

BL <- PBCredx[!duplicated(PBCredx$id), c("id", "histologic")]
BL$histologic <- factor(BL$histologic, 1:4)

PBCredx <- merge(dplyr::select(PBCredx, -histologic), 
                 BL, "id")
surv.data <- PBCredx[!duplicated(PBCredx[,'id']), ]
rm(BL)

.startmessage <- function(){
  cat("\f")
  message("PBC - Thesis Application")
  cat("Current existing directories:\n")
  fs::dir_tree(save.dir, recurse = FALSE)
  cat("\nData preview:\n")
  print(head(PBCredx))
  cat("\n(Exploratory plots previously created in ../../Explanatory).\n")
}

source('../../thesis-sims/theme_csda.R')

.startmessage()

