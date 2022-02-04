#' ###########################################################
#' A short, inconsequential program to quickly print a 
#' .png which just says --- Chapter 1 ---- and so on
#' in order to break up Figures folder a little on Overleaf
#' ###########################################################

printchapter <- function(i){
  png(filename = paste0('chapter', i, '.png'), width = 200, height = 100, units = 'mm', res = 1e3)
  plot(1,1,cex = 0,xaxt='n',yaxt='n',xlab='',ylab='',bty='n')
  text(1,1,paste0('--- Chapter ', i, ' ---'), cex=5)
  dev.off()
}

for(i in 1:5) printchapter(i)
