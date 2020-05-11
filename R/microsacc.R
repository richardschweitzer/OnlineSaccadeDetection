#============================================================
# Function microsacc() -- Microsaccade Toolbox 0.9
# (R-language Version)
# Authors: Ralf Engbert, Petra Sinn, Konstantin Mergenthaler, 
# and Hans Trukenbrod
# Date: February 20th, 2014
# MODIFIED BY RICHARD SCHWEITZER to make the function deal with NaNs
#============================================================
#  INPUT:
#  x[,1:2]		position vector
#  VFAC			  relative velocity threshold
#  MINDUR		  minimal saccade duration
#
#  OUTPUT:
#  $table[,1:7]		[(1) onset, (2) end, (3) peak velocity, 
#                 (4) horizontal component, (5) vertical component,
#			            (6) horizontal amplitude, (6) vertical amplitude ] 	
#  $radius		    parameters of elliptic threshold
#---------------------------------------------------------------------
microsacc <- function(x,VFAC=5,MINDUR=3,SAMPLING=500) {
  # Compute velocity
  v <- vecvel(x,SAMPLING=SAMPLING)
  
  # Compute threshold
  medx <- median(v[,1], na.rm = TRUE)
  msdx <- sqrt( median((v[,1]-medx)^2, na.rm = TRUE) )
  medy <- median(v[,2], na.rm = TRUE) 
  msdy <- sqrt( median((v[,2]-medy)^2, na.rm = TRUE) )
  if ( msdx<1e-10 ) {
    msdx <- sqrt( mean(v[,1]^2, na.rm = TRUE) - (mean(v[,1], na.rm = TRUE))^2 )
    if ( msdx<1e-10 ) {
      stop("msdx<realmin in microsacc")
    }
  }  
  if ( msdy<1e-10 ) {
    msdy <- sqrt( mean(v[,2]^2, na.rm = TRUE) - (mean(v[,2], na.rm = TRUE))^2 )
    if ( msdy<1e-10 ) {
      stop("msdy<realmin in microsacc")
    }
  }  
  radiusx <- VFAC*msdx
  radiusy <- VFAC*msdy
  radius <- c(radiusx,radiusy)
  
  # Apply test criterion: elliptic treshold
  test <- (v[,1]/radiusx)^2 + (v[,2]/radiusy)^2
  indx <- which(!is.na(test) & test>1)
  
  # Determine saccades
  N <- length(indx) 
  nsac <- 0
  sac <- NULL
  dur <- 1
  a <- 1
  k <- 1
  
  # Loop over saccade candidates
  while ( k<N ) {
    if ( indx[k+1]-indx[k]==1 ) {
      dur <- dur + 1
    } else {
      # Minimum duration criterion (exception: last saccade)
      if ( dur>=MINDUR ) {
        nsac <- nsac + 1
        b <- k
        sac <- rbind(sac,c(indx[a],indx[b],rep(0,5)))
      }
      a <- k+1
      dur <- 1
    }
    k <- k + 1
  }
  
  # Check minimum duration for last microsaccade
  if  ( dur>=MINDUR ) {
    nsac <- nsac + 1
    b <- k
    sac <- rbind(sac,c(indx[a],indx[b],rep(0,5)))
  }
  
  if ( nsac>0 ) {
    # Compute peak velocity, horiztonal and vertical components
    for ( s in 1:nsac ) {
      # Onset and offset for saccades
      a <- sac[s,1] 
      b <- sac[s,2]
      idx <- a:b
      # Saccade peak velocity (vpeak)
      vpeak <- max( sqrt( v[idx,1]^2 + v[idx,2]^2 ) )
      sac[s,3] <- vpeak
      # Saccade vector (dx,dy)
      dx <- x[b,1]-x[a,1] 
      dy <- x[b,2]-x[a,2] 
      sac[s,4:5] <- c(dx,dy)
      # Saccade amplitude (dX,dY)
      minx <- min(x[idx,1])
      maxx <- max(x[idx,1])
      miny <- min(x[idx,2])
      maxy <- max(x[idx,2])
      ix1 <- which.min(x[idx,1])
      ix2 <- which.max(x[idx,1])
      iy1 <- which.min(x[idx,2])
      iy2 <- which.max(x[idx,2])
      dX <- sign(ix2-ix1)*(maxx-minx)
      dY = sign(iy2-iy1)*(maxy-miny)
      sac[s,6:7] = c(dX,dY)
    }
    sac <- list(table=sac,radius=radius)
  }  else  {
    sac <- NULL
  }
  return(sac)
}
