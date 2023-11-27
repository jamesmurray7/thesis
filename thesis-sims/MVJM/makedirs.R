# #######################################
# Initialise all directories in save.dir#
# #######################################

make.fits.and.data.subdirs <- function(dirname, replaceAll = FALSE){
  current <- save.dir.file.path(dirname) # "/K, /n, ..., /vechD"
  cat(current, "\n")
  for(f in c("data", "fits")){
    f.to.create <- file.path(save.dir, dirname, f)
    if(dir.exists(f.to.create) && !replaceAll){
      cli::cli_alert_warning("The directory {f.to.create} already exists")
      cat("Type 'y' to replace the current directory with an empty one.")
      xx <- tolower(readline())
      if(xx == 'y'){
        fs::dir_delete(f.to.create)
        fs::dir_create(f.to.create)
        cli::cli_alert_success("{f.to.create}")
      }
    }else{
      fs::dir_create(f.to.create)
      cli::cli_alert_success("{f.to.create}")
    }
  }
}

# file.show("info.txt")

to.create <- c("K", "n", "r", "omega", "beta", "vechD", "C")

for(i in seq_along(to.create)){
  current <- to.create[i]
  cli::cli_alert_info("Creating directory {current}")
  if(dir.exists(save.dir.file.path(current))){
    cli::cli_alert_warning("The directory {save.dir.file.path(current)} already exists")
    cat("Type 'y' to replace the current directory with an empty one.")
    xx <- tolower(readline())
    if(xx == 'y'){
      fs::dir_delete(save.dir.file.path(current))
      fs::dir_create(save.dir.file.path(current))
    }else{
      cli::cli_alert_success("Skipped creation for {save.dir.file.path(current)}")
    }
  }else{
    fs::dir_create(save.dir.file.path(current))
  }
  cli::cli_alert("Creating 'fits' and 'data' directory...")
  make.fits.and.data.subdirs(current)
  cli::cli_rule()
}


