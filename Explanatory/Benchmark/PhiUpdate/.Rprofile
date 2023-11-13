library(gmvjoint)
library(microbenchmark)
utils::data(PBC, package = "gmvjoint")
PBC$id <- as.numeric(as.character(PBC$id))
options(
	prompt = glue::glue_col("{green PhiUpdate> }"),
	continue = "  ",
	stringsAsFactors = FALSE,
	rstudio.help.showDataPreview = FALSE
)

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
PBCredx$histologic2 <- ifelse(PBCredx$histologic %in% c("3", "4"), 1, 0)

unpack.microbenchmark <- function(x){
  stopifnot(inherits(x, "microbenchmark"))
  data.frame(what = x$expr, time = x$time/1e9)
}

cat("\f")

.startmessage <- function(){
  cat("\f")
  message("Benchmarks - Phi Update")
  print(head(PBCredx))
}

source('../../../thesis-sims/theme_csda.R')

.startmessage()
