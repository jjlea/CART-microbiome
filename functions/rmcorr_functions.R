my_rmcorr <- function(a,b,s){
  tmp <- data.frame(a=as.numeric(a), b=as.numeric(b), stringsAsFactors = F)
  tmp$s <- s
  re <- rmcorr(s, a, b, tmp)
  res_list <- list(r=re$r,p=re$p)
}
