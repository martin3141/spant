# ============================================ #
# Functions for solving least squares problems #
#                                              #
# Author: Yong Wang (yongwang@auckland.ac.nz)  #
#         Department of Statistics             #
#         University of Auckland               #
#         New Zealand                          #
# ============================================ #

#---------------------------------- #
# Nonnegative least squares (NNLS): #
#                                   #
#        Minimize     ||ax - b||    #
#        subject to   x >= 0        # 
# --------------------------------- #

nnls = function(a, b) {
  if(!is.vector(b)) b = drop(b)
  if(!is.matrix(a)) stop("a not matrix")
  m = as.integer(dim(a)[1])
  n = as.integer(dim(a)[2])
  if(length(b) != m) stop("length(b) != ncol(a)")
  storage.mode(a) = "double"
  storage.mode(b) = "double"
  x = double(n)                       # only for output
  rnorm = double(1)                   # only for output
  w = x                               # n-vector of working space
  zz = b                              # m-vector of working space
  index = integer(n)                  # n-vector index, only for output
  mode = integer(1)                   # success-failure flag; = 1, success
  .Fortran("nnls",r=a,m,m,n,b=b,x=x,rnorm=rnorm,w,zz,index=index,
           mode=mode,PACKAGE="spant")[c("x","r","b","index","rnorm","mode")]
}

pnnls = function(a, b, k=0, sum=NULL) {
  if(!is.vector(b)) b = drop(b)
  if(!is.matrix(a)) stop("a not matrix")
  m = as.integer(dim(a)[1])
  n = as.integer(dim(a)[2])
  if(!is.null(sum)) {
    if(sum <= 0) stop("Argument 'sum' must be positive or NULL")
    if(k<n) a[,(k+1):n] = a[,(k+1):n] * sum - b
    else stop("k == ncol(a) (null simplex)")
    a = rbind(a, c(rep(0,k), rep(1, n-k)))
    b = c(rep(0, m), 1)
    m = as.integer(m+1)
  }
  if(length(b) != m) stop("length(b) != ncol(a)")
  storage.mode(a) = "double"
  storage.mode(b) = "double"
  x = double(n)                       # only for output
  rnorm = double(1)                   # only for output
  w = x                               # n-vector of working space
  zz = b                              # m-vector of working space
  index = integer(n)                  # n-vector index, only for output
  mode = integer(1)                   # success-failure flag; = 1, success
  k = as.integer(k)
  r = .Fortran("pnnls",r=a,m,m,n,b=b,x=x,rnorm=rnorm,w,zz,index=index,
      mode=mode,k=k,PACKAGE="spant")
  r$r = r$r[1:min(m,n),]
  if(!is.null(sum)) {
    t = sum(r$x[(r$k+1):n])
    r$x = r$x / t
    r$x[(r$k+1):n] = r$x[(r$k+1):n] * sum
    r$rnorm = sqrt( pmax((r$rnorm/t)^2 - (1 - 1/t)^2, 0) )
  }
  r[c("x","r","b","index","rnorm","mode","k")]
}