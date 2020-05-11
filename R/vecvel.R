#============================================================
# Function vecvel() -- Microsaccade Toolbox 0.9
# (R-language Version)
# Authors: Ralf Engbert, Petra Sinn, Konstantin Mergenthaler, 
# and Hans Trukenbrod
# Date: February 20th, 2014
#============================================================
vecvel <- function(x,SAMPLING=500,TYPE=2) {
  d <- dim(x)
  N <- d[1]
  v <- matrix(rep(0,2*N),ncol=2)
  
  if ( TYPE==2 ) {
    v[3:(N-2),] <- SAMPLING/6*(x[5:N,] + x[4:(N-1),] - x[2:(N-3),] - x[1:(N-4),])
    v[2,] = SAMPLING/2*(x[3,] - x[1,])
    v[(N-1),] = SAMPLING/2*(x[N,] - x[(N-2),])   
  }  else  {
    v[2:(N-1),] <- SAMPLING/2*(x[3:N,] - x[1:(N-2),])
  }
  return(v)
}