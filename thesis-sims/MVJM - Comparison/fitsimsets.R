rm(list=ls())
source(".Rprofile")
load.dir <- save.dir.file.path("data")
save.dir <- save.dir.file.path("fits")
N <- 100
files <- dir(load.dir)[3] # change if needed
for(f in files){
  assign("sim.sets", get(load(file.path(load.dir, f))))
  K <- as.numeric(gsub("^K", "", gsub("\\.RData", "", f)))
  out.file.gmv <- save.dir.file.path(paste0("K", K,"_gmv.RData"), save.dir)
  out.file.jML <- save.dir.file.path(paste0("K", K,"_jML.RData"), save.dir)
  cli::cli_alert_info("Starting fits on {f}")
  cli::cli_inform("Starting gmvjoint fits, this will be saved in {out.file.gmv}...")
  fits <- fit.all(sim.sets, method = "gmvjoint", K = K)
  save(fits, file = out.file.gmv)
  cli::cli_alert_success("Done\n")
  rm(fits)
  cli::cli_inform("Starting joineRML fits, this will be saved in {out.file.jML}...")
  fits <- fit.all(sim.sets, method = "joineRML", K = K)
  save(fits, file = out.file.jML)
  cli::cli_alert_success("Done\n")
}
gc()