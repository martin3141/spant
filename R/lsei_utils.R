## Utility functions

## Index-finding

indx = function(x, v) {
  storage.mode(x) = storage.mode(v) = "double"
  m = length(x)
  n = length(v)
  ind = integer(m)
  .C("indx", x, m, v, n, ind=ind, PACKAGE="spant")["ind"]$ind
}

## # Much slower than indx().
## 
## indx2 = function(x, v) {
##   ox = order(x)
##   vx = c(v,x)
##   o = order(vx)
##   indo = c(rep(1, length(v)), rep(0, length(x)))[o]
##   indx = double(length(x))
##   indx[ox] = cumsum(indo)[indo == 0] 
##   indx
## }

matMaxs = function(x, dim=1) {
  if(length(x) == 0) return(NULL)
  x.mode = storage.mode(x)
  n = nrow(x)
  m = ncol(x)
  v = if(dim == 1) double(n) else double(m)
  storage.mode(dim) = storage.mode(n) = storage.mode(m) = "integer"
  x[x == Inf] = 1e308
  x[x == -Inf] = -1e308
  storage.mode(x) = "double"
  v = .C("matMaxs", x, n, m, v=v, dim, PACKAGE="spant")["v"]$v
  v[v > 9e307] = Inf
  v[v < -9e307] = -Inf
  storage.mode(v) = x.mode
  v
}
