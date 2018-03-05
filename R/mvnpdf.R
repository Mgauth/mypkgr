# mvnpdf1 <- function(x,mean,varcovM,Log=FALSE){
#   # x : matrice p*n
#   n <- dim(varcovM)[2]
#   p <- dim(varcovM)[1]
#   y <- NULL
#   for (j in 1:n){
#     res <- ((1/(2*pi)^p)*(det(varcovM))^0.5)*exp(-0.5*t(x[,j]-mean)%*%solve(varcovM)%*%(x[,j]-mean))
#     y <- c(y,res)
#     }
#   l <- list()
#
#   if (Log==TRUE){
#     y <- log(y)
#   }
#   l[[1]] <- x
#   l[[2]] <- y
#   return(l)
# }
#



######

#' Densite d une loi normale multivariee
#'
#' @param x matrice taille p*n (n observations et p dimensions)
#' @param mean vecteur de moyennes
#' @param varcovM une matrice de variance-covariance
#' @param Log booleen pour retourner le log
#'
#' @return renvoie une liste contenant la matrice x ainsi qu un vecteur
#' des images des points de x par la fonction de densite de la variable
#' aleatoire de loi normale multivariee consideree.
#' @export
#'
#' @examples
#' mvnpdf(x=matrix(1.96), log=FALSE)
#' dnorm(1.96)
#'
mvnpdf <- function(x, mean =  rep(0, nrow(x)),
                   varcovM = diag(nrow(x)), Log = TRUE) {
  n <- ncol(x)
  p <- nrow(x)
  x0 <- x - mean
  Rinv <- solve(varcovM)
  LogDetvarcovM <- log(det(varcovM))

  y <- NULL
  for (j in 1:n) {
    yj <- - p/2 * log(2*pi) - 0.5 * LogDetvarcovM -
      0.5 * t(x0[, j]) %*% Rinv %*% x0[, j]
    y <- c(y, yj)
  }

  if (!Log) {
    y <- exp(y)
  }

  return(list(x=x, y=y))
}
