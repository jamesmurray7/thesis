# #######################################
# Initialise all directories in save.dir#
# #######################################

does.this.exist <- function(path){
  as.logical(fs::dir_exists(save.dir.file.path(path)))
}

names.of.dirs <- c("data", "fits")

sapply(names.of.dirs, function(f){
  if(does.this.exist(f)){
    warning(sQuote(save.dir.file.path(f)), " alredy exists!")
  }else{
    fs::dir_create(save.dir.file.path(f))
  }
})
