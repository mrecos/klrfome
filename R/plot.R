#' K_corrplot
#'
#' @param K - Matrix
#' @param train_data - data.frame
#' @param clusters - integer
#'
#' @return - a plot
#' @export
#'
K_corrplot <- function(K,train_data,clusters=4){
  # should be a method instead of a function probably
  colnames(K) <- names(train_data)
  rownames(K) <- names(train_data)
  col3 <- colorRampPalette(c("red", "white", "blue"))
  corrplot::corrplot(K,tl.cex = 0.5, tl.col = "black",
                     order="hclust", col=col3(100), cl.lim=c(0,1),
                     addrect = clusters)
}
