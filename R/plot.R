#' K_corrplot
#' 
#' `K_corrplot()` is a wrapper around `corrplor::corrplot()` that returns a corrplot for visualizing the similarity matrix of site and background bags.
#' 
#' This function is a wrapper of `corrplot::corrplot()` with defaults and hierarchical clustering order. The inputs are the similarity kernel matirx `K`, the `train_data` used to creat `K`, and the number of clusters to display. The `train_data` is only used to procure the labels for site and background bags.
#'
#' @param K [K] - Similarity Kernel Matrix
#' @param train_data [list] - training data
#' @param clusters [scalar] - Number of clusters
#'
#' @return - a correlation matrix object
#' @export
#' @importFrom corrplot corrplot
#'
#' @examples
#' \dontrun{
#'##### Logistic Mean Embedding KRR Model
#' #### Build Kernel Matrix
#' K <- build_K(train_data, sigma = sigma, dist_metric = dist_metric, progress = FALSE)
#' #### Train KLR model
#' train_log_pred <- KLR(K, train_presence, lambda, 100, 0.001, verbose = 2)
#' ### Plot K Matrix
#' K_corrplot(K, train_data, clusters = 4)
#'}
#'
K_corrplot <- function(K,train_data,clusters=4){
  # should be a method instead of a function probably
  if (!requireNamespace("corrplot", quietly = TRUE)) {
    stop("corrplot needed for this function to work. Please install it.",
         call. = FALSE)
  }
  colnames(K) <- names(train_data)
  rownames(K) <- names(train_data)
  col3 <- colorRampPalette(c("red", "white", "blue"))
  corrplot::corrplot(K,tl.cex = 0.5, tl.col = "black",
                     order="hclust", col=col3(100), cl.lim=c(0,1),
                     addrect = clusters)
}
