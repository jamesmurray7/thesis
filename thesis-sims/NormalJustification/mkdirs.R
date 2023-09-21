existing.dirs <- gsub("\\/\\/","/",
                      dir(save.dir, pattern = "Justification$", full.names = T))
dirs.to.create <- sapply(.valid.families, family.dir.name, USE.NAMES = F)
dirs.to.create <- save.dir.file.path(dirs.to.create)
dirs.to.create <- setdiff(dirs.to.create, existing.dirs)

if(length(dirs.to.create) == 0L){
  cat(cli::rule(col = "mediumseagreen"))
  rr <- paste(rep("-", cli::console_width()/2.5),collapse="",sep="")
  message(rr, "> All directories created already!")
  cat(cli::rule(col = "mediumseagreen"))
}else{
  for(f in dirs.to.create){
    invisible(fs::dir_create(f))
    cat(sprintf("\n%s created\n", f))
  }
}
