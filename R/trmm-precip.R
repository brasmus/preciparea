## TRMM data
## Time series in area with precipitation
## Time series with total precipitation
## Time series with mean precipitation for precipitation area
##
## A trend in the area with precipitation means that the mean intensity changes
## A_e E = A_p P
## Compare annual cycle in precipitation area with intensity and E.

require(esd)

## Indices for integrated evaporation from ERAINT reanalysis
## used to prepare data for data(Es)
evaporation <- function(fname='~/Downloads/ERAINT_t2m.nc',R.s=18.02) {
  require(esd)
  print('50S-50N')
  t2m <- retrieve(fname,lat=c(-50,50))
  es <- C.C.eq(mask(t2m,land = TRUE))
  ## p = rho R_s T -> rho = p / (R_s * T)
  rho <- es/(R.s * (t2m + 273.15))
  es.50S50N <- aggregate.area(es,FUN='mean')
  rho.50S50N <- aggregate.area(rho)
  rm('t2m','es','rho'); gc(reset=TRUE)
  
  print('90S-50s')
  t2m <- retrieve(fname,lat=c(-90,-50))
  es <- C.C.eq(mask(t2m,land = TRUE))
  ## p = rho R_s T -> rho = p / (R_s * T)
  rho <- es/(R.s * (t2m + 273.15))
  es.90S50S <- aggregate.area(es,FUN='mean')
  rho.90S50S <- aggregate.area(rho)
  rm('t2m','es','rho'); gc(reset=TRUE)
  
  print('50N-90N')
  t2m <- retrieve(fname,lat=c(50,90))
  es <- C.C.eq(mask(t2m,land = TRUE))
  ## p = rho R_s T -> rho = p / (R_s * T)
  rho <- es/(R.s * (t2m + 273.15))
  es.50N90N <- aggregate.area(es,FUN='mean')
  rho.50N90N <- aggregate.area(rho)
  rm('t2m','es','rho'); gc(reset=TRUE)
  
  print('combine')
  Es <- combine.stations(es.50S50N,es.90S50S,es.50N90N)
  Rho <- combine.stations(rho.50S50N,rho.90S50S,rho.50N90N)
  #plot(Es,map.show=FALSE)
  save(Es,file='Es.rda')
  save(Rho,file='Rho.rda')
  invisible(Es)
}


## Function to estimate the precipitation area for monthly mean precipitation
# TRMM: 180W-180E; MERRA-2: 180W-180E; ERAINT: 0-360E

monthP2area <- function(fname='~/TRMM-monthly/trmm-precip-mon.nc',param='precip',x0=1,dx=0.01,lons=c(-180,180)) {
  for (j in seq(lons[1],lons[2]-10,by=10)) {
    if (j < 0) dj <- -dx else dj <- dx
    print(c(j,j+10-dj))
    y <- retrieve(fname,param=param,lon=c(j,j+10-dj),lat=c(-50,50))
    if (j==lons[1]) {
      Pa <- zoo(aggregate.area(y,FUN='area',x0=x0)) 
      A <- coredata(aggregate.area(y,FUN='area'))[1]
    } else {
      Pa <- Pa + zoo(aggregate.area(y,FUN='area',x0=x0))
      A <- A + coredata(aggregate.area(y,FUN='area'))[1]
    }
    #plot(Pa)
  }
  if (!is.null(attr(y,'source'))) attr(Pa,'source') <- attr(y,'source') 
  if (!is.null(attr(y,'model_id'))) attr(Pa,'model_id') <- attr(y,'model_id') 
  attr(Pa,'total area') <- A
  attr(Pa,'threshold') <- x0
  attr(Pa,'model_res') <- paste(diff(lon(y))[1],diff(lat(y))[1],sep='x')
  save(Pa,file='month2Parea.rda')
  rm('y'); gc(reset=TRUE)
  invisible(Pa)
}



## A function that reads a netCDF in latitude strips and accumulates the amount for computers with small memory
## Returns the area sum for three latitude bands: 90S-50S, 50S-50N, and 50N-90N
# TRMM: 180W-180E; MERRA-2: 180W-180E; ERAINT: 0-360E
globalsum <- function(fname,param=NULL,FUN='sum',dlon=30,dx=0.01,lons=c(-180,180)) {
  for (y in seq(lons[1],lons[2]-dlon,by=dlon)) {
    if (y + dlon < 0) dj <- -dx else dj <- dx
    if (is.null(param)) x <- retrieve(fname,lon=c(y,y+dlon-dj)) else {
                        x <- retrieve(fname,param=param[1],lon=c(y,y+dlon-dj))
      if (length(param)> 1) 
        for (i in 2:length(param)) 
          coredata(x) <- coredata(x) + coredata(retrieve(fname,param[i],lon=c(y,y+dlon-dj)))
    } 
    
    if (y == lons[1]) {
      X.90s50s <- aggregate.area(subset(x,is=list(lon=lons,lat=c(-90,-50))),FUN=FUN)
      X.50s50n <- aggregate.area(subset(x,is=list(lon=lons,lat=c(-50,50))),FUN=FUN)
      X.50n90n <- aggregate.area(subset(x,is=list(lon=lons,lat=c(50,90))),FUN=FUN)
    } else {
      X.90s50s <- X.90s50s + aggregate.area(subset(x,is=list(lon=lons,lat=c(-90,-50))),FUN=FUN)
      X.50s50n <- X.50s50n + aggregate.area(subset(x,is=list(lon=lons,lat=c(-50,50))),FUN=FUN)
      X.50n90n <- X.50n90n + aggregate.area(subset(x,is=list(lon=lons,lat=c(50,90))),FUN=FUN)
    }
    print(c(y,y+dlon-dj, round(c(mean(X.50s50n), mean(X.90s50s),mean(X.50n90n))/1.0e12)))
  }
  X <- merge(X.50s50n,X.90s50s,X.50n90n)
  invisible(X)
}

## Prepare/extract information from ERAINT
moistflux <- function(fname='~/Downloads/ERAINT_Qs.nc') {
  for (it in seq(1979,2016,by=10)) {
    print(it)
    x <- retrieve('~/Downloads/ERAINT_Qs.nc',it=c(it,it+9))
    qs <- zoo(aggregate.area(x,FUN='sum')/1.0e6)
    if (it==1979) Qs <- qs else Qs <- c(Qs,qs)
  }
  save(Qs,file='Qs.rda')
  attr(Qs,'unit') <- 'Gton/(day*m^2)'
  invisible(Qs)
}

## Prepare/extract information from ERAINT
columnH2O <- function(fname='~/Downloads/ERAINT_Q-columntotal.nc') {
  #require(esd)
  #h2o <- retrieve(fname,lat=c(-50,50))
  #h2o.50S50N <- aggregate.area(h2o)
  #rm('h2o'); gc(reset=TRUE)
  #
  #h2o <- retrieve(fname,lat=c(-90,-50))
  #h2o.90S50S <- aggregate.area(h2o)
  #rm('h2o'); gc(reset=TRUE)
  #
  #h2o <- retrieve(fname,lat=c(50,90))
  #h2o.50N90N <- aggregate.area(h2o)
  #rm('h2o'); gc(reset=TRUE)
  
  #H2O <- combine.stations(h2o.50S50N,h2o.90S50S,h2o.50N90N)
  #plot(Es,map.show=FALSE)
  # original units = "kg m**-2"
  H2O <- globalsum(fname=fname)*1.0e6 ## global sum returns kg/km^2
  attr(H2O,'unit') <- rep('kg',3)
  save(H2O,file='H2O.rda')
  invisible(H2O)
}



## function that reads the MODIS data and estimate indices for precip area, total precip,...
## prepare for data(Parea)
retrievemodis <- function(path='~/TRMM',
                          pattern='3B42_Daily',
                          param='precipitation',
                          nx=1440,ny=400,a = 6378, x0 = 1) {
  require(ncdf4)  
  fnames <- list.files(path=path,pattern=pattern,full.names=TRUE)
  print(fnames)
  nt <- length(fnames)
  Ap <- rep(NA,nt); Af <- Ap; Am <- Ap; Pt <- Ap; mu <- Ap; Pm <- Ap; nv <- Ap; Al <- Ap; Pl <- Ap
  t <- rep(NA,nt)
  
  if (file.exists('retrievemodis.tmp.rda')) load('retrievemodis.tmp.rda') else i1 <- 1
  
  for (i in i1:nt) {
    fname <- fnames[i]
    print(fname)
    ncid <- nc_open(fname)
    x <- ncvar_get(ncid,varid=param)
    x <- t(x)
    dim(x) <- c(1,nx*ny)
    tim <- ncatt_get(ncid,varid = 0,attname='EndDate')
    t[i] <- tim$value
    lon <- ncvar_get(ncid,'lon')
    lat <- ncvar_get(ncid,'lat')
    nc_close(ncid)
    attr(x,'longitude') <- lon
    attr(x,'latitude') <- lat
    x <- as.field(zoo(x,order.by=t[i]),lon=lon,lat=lat,param='precip',unit='mm/day',src='TRMM')
    
    if (i==1) {
      ## Create a mask for land areas.
      test <- subset(retrieve('air.mon.mean.nc'),it=1)
      test <- regrid(test,is=list(lon=lon,lat=lat))
      A.tot <- as.numeric(aggregate.area(test,FUN='area'))[1]
      print(paste('A.tot fraction of planetary area: ',round(100*A.tot/(4*pi*a^2)),'%'))
      ## Sanity check
      if ( (A.tot < 2*pi*a^2) | (A.tot > 4*pi*a^2)) {
        print('Total estimated area between 50S-50N is less than 30S-30N or larger than planetary area')
        browser()
      }
      ocean <- mask(test,land=TRUE)
      land <- mask(test,land=FALSE)
      ocean <- coredata(ocean)
      ocean[is.finite(ocean)] <- 1
      coredata(test) <- ocean
      A.ocean <- as.numeric(aggregate.area(test,FUN='area'))[1]
      land <- coredata(land)
      land[is.finite(land)] <- 1
      coredata(test) <- land
      A.land <- as.numeric(aggregate.area(test,FUN='area'))[1]
    }
    
    xl <- x; coredata(xl) <- coredata(xl)*land
    xo <- x; coredata(xo) <- coredata(xo)*ocean
    
    ## Precipitation area: land + ocean
    Ap[i] <- aggregate.area(x,FUN='area',x0=x0)
    ## Total precipitation
    Pt[i] <- aggregate.area(x,FUN='sum')
    nv[i] <- sum(is.finite(x))
    
    ## Over oceans
    coredata(x) <- xo
    Pm[i] <- aggregate.area(xo,FUN='sum')
    Am[i] <- aggregate.area(xo,FUN='area',x0=x0)
    ## Over land
    Pl[i] <- aggregate.area(xl,FUN='sum')
    Al[i] <- aggregate.area(xl,FUN='area',x0=x0)
    
    ## Wet mean: discard all grid boxes with no/negligible rainfall
    coredata(x)[x < x0] <- NA
    mu[i] <- aggregate.area(x,FUN='mean')
    if ( (mu[i] < 1) | (mu[i]> 100) ) browser() ## Check
    print(c(100*Ap[i]/A.tot,mu[i]))
    
    ## Save temporarily the results every 350 iteration
    if (i %% 350==0) {
      i1 <- i+1
      save(i1,Ap, Am, Al, mu, Pt, Pm, Pl, nv,t,A.tot,A.land,A.ocean,
           x0,land,ocean,file='retrievemodis.tmp.rda')
    }
  }
  
  file.remove('retrievemodis.tmp.rda')
  
  X <- cbind(Ap, Am, Al, mu, Pt, Pm, Pl, nv)
  colnames(X) <- c('A.precip','A.precip.marine','A.precip.land',
                   'mu','P.tot','P.marine','P.land','nv')
  Y <- zoo(X,order.by=as.Date(t))
  attr(Y,'longitudes') <-range(lon)
  attr(Y,'latitudes') <-range(lat)
  attr(Y,'unit') <- c('km^2','km^2','km^2','mm','kilo-ton','kilo-ton','kilo-tonton','count')
  attr(Y,'longname') <-c('Precipitation area',
                         'Precipitation area over ocean','Precipitation area over land',
                         'wet-day mean preciptation',
                         'total precipitation','total precipitation over oceans',
                         'total precipitation over land','number of valid data')
  attr(Y,'src') <- 'TRMM'
  attr(Y,'precip.threshold') <- x0
  attr(Y,'total area') <- A.tot
  attr(Y,'ocean area') <- A.ocean
  attr(Y,'land area') <- A.land
  attr(Y,'url')='http://disc2.gesdisc.eosdis.nasa.gov/data/TRMM_L3/TRMM_3B42_Daily.7/'
  attr(Y,'history') <- 'trmm-precip.R'
  attr(Y,'author') <-'rasmus.benestad@met.no'
  Parea <- Y
  fname <- 'Parea.rda'
  save(Parea,file=fname)
  return(Y)
}

# wget --load-cookies ~/.urs_cookies --save-cookies ~/.urs_cookies --keep-session-cookies -r -c -nH -nd -np -A nc4 "https://disc2.gesdisc.eosdis.nasa.gov/data/TRMM_L3/TRMM_3B42_Daily.7/1998/01"


## function that reads the MODIS data and estimate indices for precip area, total precip,...
## prepare for data(Parea.merra)
retrieveMERRA <- function(path='~/MERRA',
                          pattern='MERRA2_',reg50s50n=TRUE,
                          param=c('PRECCUCORR','PRECLSCORR','PRECSNOCORR'),
                          a = 6378, x0 =1) {
  require(esd)
  require(ncdf4)  
  fnames <- list.files(path=path,pattern=pattern,full.names=TRUE)
  fnames <- fnames[nchar(fnames)==56]
  print(fnames)
  nt <- length(fnames)
  
  for (i in 1:nt) {
    fname <- fnames[i]
    print(fname)
    if (nchar(fname)==56) {
      ncid <- nc_open(fname)
      torg <- ncatt_get(ncid,'time',attname='begin_date')
      t1 <- ncatt_get(ncid,0,attname='RangeBeginningDate')
      t2 <- ncatt_get(ncid,0,attname='RangeEndingDate')
      #tim <- ncvar_get(ncid,'time')  ## unit = "minutes since 1998-06-01 00:30:00" -> days - not right!
      origin <- as.Date(paste(trunc(torg$value/10000),
                            trunc(100*trunc(torg$value/100)-10000*trunc(torg$value/10000))/100,
                            torg$value-100*trunc(torg$value/100),sep='-'))
      tim <- seq(as.Date(t1$value),as.Date(t2$value),length=length(tim))
      for (j in 1:length(param)) {
        x <- ncvar_get(ncid,param)*86400
        if (j==1) X <- x else X <- X + x
      }
      lon <- ncvar_get(ncid,'lon')
      lat <- ncvar_get(ncid,'lat')
      d <- dim(X)
      nc_close(ncid)
    }
    
    ## Need to reorder the matrix to create a field object 
    ## From aggregate
    xy <- rep(lon,d[2]) + 1000*sort(rep(lat,d[1]))
    dim(X) <- c(d[1]*d[2],d[3])
    str(X); str(xy); str(lon); str(lat)
    x <- X[order(xy),]
    y <- as.field(zoo(t(as.matrix(x)),order.by=as.Date(tim,origin='1998-06-01')),
                  param='precip',unit='mm/day',
                  lon=lon,lat=lat,src='MERRA-2')
    class(y) <- c('field','day','zoo')
    ## Same region as TRMM
    Y <- y
    if (reg50s50n) y <- subset(y,is=list(lon=c(-180,360),lat=c(-50,50)))
    
    z <- y; zc <- y
    
    if (i==1) {
      Tim <- tim
      ## Create a mask for land areas.
      test <- subset(retrieve('air.mon.mean.nc'),it=1)
      test <- regrid(test,is=list(lon=lon,lat=lat))
      # map(all,type='fill')  # sanity check!
      A.tot <- as.numeric(aggregate.area(y,FUN='area'))[1]
      print(paste('A.tot: ',A.tot,'<=',as.numeric(aggregate.area(Y,FUN='area'))[1]))
      if (A.tot > as.numeric(aggregate.area(Y,FUN='area'))[1]) browser()
      A.ocean <- as.numeric(aggregate.area(mask(y,land=TRUE),FUN='area')[1])
      A.land <- as.numeric(aggregate.area(mask(y,land=FALSE),FUN='area')[1])
      print(paste('A.tot fraction of planetary area: ',round(100*A.tot/(4*pi*a^2)),'%'))
      
      ## Total precipitation
      Pt <- zoo(aggregate.area(y,FUN='sum'))
      print(paste('Pt: ',round(mean(Pt,na.rm=TRUE)),'<=',
                  round(mean(zoo(aggregate.area(Y,FUN='sum')),na.rm=TRUE)),' (global)'))
      if (mean(Pt,na.rm=TRUE) > mean(aggregate.area(Y,FUN='sum'),na.rm=TRUE)) browser()
      ## Precipitation sum over oceans
      Pm <- zoo(aggregate.area(mask(y,land=TRUE),FUN='sum'))
      ## Precipitation sum over land
      Pl <- zoo(aggregate.area(mask(y,land=FALSE),FUN='sum'))
      ## Precipitation area
      Ap <- zoo(aggregate.area(y,FUN='area',x0=x0))
      print(paste('Ap: ',round(mean(Ap,na.rm=TRUE)),'<=',
                  round(mean(aggregate.area(Y,FUN='area',x0=x0),na.rm=TRUE)),' (global)'))
      if (mean(Ap,na.rm=TRUE) > mean(aggregate.area(Y,FUN='area',x0=x0),na.rm=TRUE)) browser()
      ## Precipitation area over oceans
      Am <- zoo(aggregate.area(mask(y,land=TRUE),FUN='area',x0=x0))
      ## Precipitation area over land
      Al <- zoo(aggregate.area(mask(y,land=FALSE),FUN='area',x0=x0))
     
      ## Precipitation intensity
      yp <- y; ypc <- coredata(yp); ypc[ypc < x0] <- NA; ypc -> coredata(yp)
      mu <- zoo(aggregate.area(yp,FUN='mean'))
    } else {
      Tim <- c(Tim,tim)
      ## Total precipitation
      Pt <- c(Pt,zoo(aggregate.area(y,FUN='sum')))
      print(paste('Pt: ',round(mean(aggregate.area(y,FUN='sum'),na.rm=TRUE)),'<=',
                         round(mean(aggregate.area(Y,FUN='sum'),na.rm=TRUE)),' (global)'))
      ## Precipitation sum over oceans
      Pm <- c(Pm,zoo(aggregate.area(mask(y,land=TRUE),FUN='sum')))
      ## Precipitation sum over land
      Pl <- c(Pl,zoo(aggregate.area(mask(y,land=FALSE),FUN='sum')))
      ## Precipitation area
      Ap <- c(Ap,zoo(aggregate.area(y,FUN='area',x0=x0)))
      print(paste('Ap: ',round(mean(aggregate.area(y,FUN='area',x0=x0),na.rm=TRUE)),'<=',
                         round(mean(aggregate.area(Y,FUN='area',x0=x0),na.rm=TRUE)),' (global)'))
      ## Precipitation area over oceans
      Am <- c(Am,zoo(aggregate.area(mask(y,land=TRUE),FUN='area',x0=x0)))
      ## Precipitation area over land
      Al <- c(Al,zoo(aggregate.area(mask(y,land=FALSE),FUN='area',x0=x0)))
      ## Precipitation intensity
      yp <- y; ypc <- coredata(yp); ypc[ypc < x0] <- NA; ypc -> coredata(yp)
      mu <- c(mu,zoo(aggregate.area(yp,FUN='mean')))
    }
    plot(merge(Ap,Pt))
  }
  
  Y <- merge(Ap, Am, Al, mu, Pt, Pm, Pl)
  colnames(Y) <- c('A.precip','A.precip.marine','A.precip.land',
                   'mu','P.tot','P.marine','P.land')
  attr(Y,'longitudes') <-range(lon(y))
  attr(Y,'latitudes') <-range(lat(y))
  attr(Y,'unit') <- c('km^2','km^2','km^2','mm','kilo-ton','kilo-ton','kilo-ton')
  attr(Y,'longname') <-c('Precipitation area',
                         'Precipitation area over ocean','Precipitation area over land',
                         'wet-day mean preciptation',
                         'total precipitation','total precipitation over oceans',
                         'total precipitation over land','number of valid data')
  attr(Y,'src') <- 'MERRA'
  attr(Y,'precip.threshold') <- x0
  attr(Y,'total area') <- A.tot
  attr(Y,'ocean area') <- A.ocean
  attr(Y,'land area') <- A.land
  attr(Y,'url')='https://disc.gsfc.nasa.gov/uui/datasets/M2TUNXLFO_V5.12.4/summary?keywords=%22MERRA-2%22'
  attr(Y,'history') <- 'trmm-precip.R'
  attr(Y,'author') <-'rasmus.benestad@met.no'
  Parea.merra <- Y
  fname <- 'Parea.merra.rda'
  save(Parea.merra,file=fname)
  return(Y)
}


## Function that aggregates daily TRMM data to monthly mean values
## The original data was not stored in the CF convention and the monthly mean data
## were not found for ready download in the netCDF format.
TRMM2monthly <- function(path='~/TRMM',
                         pattern='3B42_Daily',
                         param='precipitation',yr1=1998,yr2=2016,
                         nx=1440,ny=400) {
  require(ncdf4)  
  fnames <- list.files(path=path,pattern=pattern,full.names=TRUE)
  #print(fnames)
  nt <- length(fnames)
  nyrs <- 1 + yr2 - yr1
  nym <- 12*nyrs
  for (i in 1:nym) {
    yyyy <- yr1 + (i-1)%/% 12
    mm <- i %% 12; if (mm==0) mm <- 12
    if (mm < 10) mm <- paste(0,mm,sep='')
    ymfnames <- fnames[grep(paste('3B42_Daily.',yyyy,mm,sep=''),fnames)]
    #print(ymfnames)
    X <- NULL
    nf <- 0
    
    for (fname in ymfnames) {
      ncid <- nc_open(fname)
      print(fname)
      x <- ncvar_get(ncid,param)
      if (is.null(X)) { 
        X <- t(x) 
        tim <- ncatt_get(ncid,varid = 0,attname='EndDate')
        time <- as.Date(tim$value)
        lon <- ncvar_get(ncid,'lon')
        lat <- ncvar_get(ncid,'lat')
        nc_close(ncid)
      } else X <- X + t(x)
      nf <- nf + 1
    }
    rm('ncid','tim'); gc(reset=TRUE)
    X <- X/nf
    print(c(yyyy,mm,nf))
    dim(X) <- c(nx,ny,1)
    #image(lon,lat,X[,,1]); lines(geoborders)
    X[!is.finite(X)] <- -99
    
    dimLon <- ncdim_def( "longitude", "degrees_east", lon )
    dimLat <- ncdim_def( "latitude", "degrees_north", lat )
    dimtim <- ncdim_def( "time", "days since 1970-01-01", as.numeric(time), unlim=TRUE )
    varPr <- ncvar_def("precip", "mm/day", list(dimLon,dimLat,dimtim), -99, 
                       longname="monthly_mean_precipitation", prec="float")
    ## Save the monthly mean data in a new directory for later concatination by NCO:
    if (!file.exists("~/TRMM-monthly/")) dir.create("~/TRMM-monthly/")
    ncnew <- nc_create( paste("~/TRMM-monthly/trmm-precip",yyyy,'-',mm,".nc",sep=""), varPr )
    
    # Write some values to this variable on disk.
    ncvar_put( ncnew, varPr, X )
    nc_close(ncnew)
    rm('x','X','ncnew','varPr','lon','lat', 'dimLon','dimLat','dimtim'); gc(reset=TRUE)
    
  }
  
  #attr(Z,'longitude') <- lon
  #attr(Z,'latitude') <- lat
  #Y <- as.field(Z,lon=lon,lat=lat,param='precip',unit='mm/day',src='TRMM')
  #save(Y,file=paste('TRMM-3B42_monthly',yr1,'-',yr2,'.nc',sep=''))
  #return(Y)
}


PareaCMIP <- function(path='~/CMIP5.monthly/rcp45/',pattern='pr_Amon_ens_',param='pr',x0=4,dx=0.01,lons=c(-180,180)) {
  fnames <- list.files(path=path,pattern=pattern,full.names=TRUE)
  A.tot <- rep(NA,length(fnames))
  model.id <- rep('NA',length(fnames)); model.res <- model.id
  for (i in 1:length(fnames)) {
    print(fnames[i])
    y <- monthP2area(fname=fnames[i],param=param,x0=x0,dx=dx,lons=lons)
    if (i==1) Y <- y else Y <- merge(Y,y)
    A.tot[i] <- attr(y,'total area')
    model.id[i] <- attr(y,'model_id')
    model.res[i] <- attr(y,'model_res')
    plot(100*Y/A.tot[1:i],ylab='%',plot.type='single',col=rgb(0,0,0,0.2)); grid()
  }
  A.tot -> attr(Y,'total area')
  x0 -> attr(Y,'threshold')
  model.id -> attr(Y,'model_id')
  model.res -> attr(Y,'model_res')
  Parea.cmip5 <- Y
  save(Parea.cmip5,file='Parea.cmip5.rda')
  invisible(Y)
}


readGDCN <- function(name,type="PRCP") {
  # Reads the ASCII files with daily GDCN data:
  cnames <- c(paste("day.",1:31,sep=""),paste("flags.day.",1:31,sep=""))
  reshufle <- rep(0,62); ii <- 0
  for (i in seq(1,62,by=2)) {
    ii <- ii + 1
    reshufle[i:(i+1)] <- seq(ii,62,by=31)
  }
  cnames <- cnames[reshufle]
  x <- read.fwf(name,widths=c(3,8,4,2,4,rep(c(5,2),31)),
                col.names=c("country","stnr","year","month","type",cnames))
  ipick <- is.element(x$type,type)
  if (sum(ipick)>0) {
    dat <- as.matrix(x[ipick,seq(6,67,by=2)])*0.1
    dat[dat < -99] <- NA
    attr(dat,"Data source") <- "GDCN" 
    attr(dat,"year") <- x$year[ipick]
    attr(dat,"month") <- x$month[ipick]
    attr(dat,"URL") <- "http://www.ncdc.noaa.gov/oa/climate/research/gdcn/gdcn.html"
    attr(dat,"Station_number") <- x$stnr[1]
    attr(dat,"Country_code") <- x$country[1]
    attr(dat,"original data name") <- name
    attr(dat,"flags") <- x[ipick,seq(6,66,by=2)]
    attr(dat,"history") <- "read with readGDCN - R-script."
    attr(dat,"Observations") <- rownames(table(x$type))
  } else dat <- NA
  invisible(dat)
}


## Analysis of the number of instantaneous precipitation recorded by raingauge data
## The GDCN data is stored in a very clumsy way, which makes the reading cumbersome...
raingauges <- function(x0=1,tim=seq(as.Date('1960-01-01'),as.Date('2016-01-01'),by=1),
                       max.stations=NULL,
                       N.min=40*360,path="~/GDCN/") {
  # Main routine reading in, analysing, plotting, and storing the
  # quantiles from the GDCN data set.
  
  stationlist <- list.files(path,pattern=".dly",full.names = TRUE)
  load('preciparea/data/gdcn.inv.rda')
  #data(gdcn.inv,envir = environment())
  
  nt <- length(tim); 
  if (!is.null(max.stations)) N <- max.stations else N <- length(stationlist)
  totP <- rep(0,nt); fw <- totP; Ng <- totP
  lat <- rep(NA,N); lon <- lat; alt <- lat; stnr <- lat; cntr <- lat; nval <- lat 
  start <- lat; end <- lat; loc <- rep('',N)
  data("geoborders")
  plot(geoborders$x,geoborders$y,type="l")
  grid()
  
  # Check if there are unfinished results
  if (file.exists("raingaugedata.temp.rda")) {
    print("read intermediate results from previous run")
    load("raingaugedata.temp.rda")
  } else {
    i1 <- 1
    i <- 1
  }
  
  # Loop through all station files
  
  for (ii in i1:N) {
    X <- readGDCN(stationlist[ii])
    imatch <- is.element(gdcn.inv$stnr,attr(X,"Station_number"))
    doit <- TRUE 
    ## If specified, only include stations with minimum valid data
    if ( (!is.null(N.min)) & (length(dim(X))==2) ) 
      if (sum(is.finite(X[attr(X,'year') %in% year(tim),])) < N.min) doit <- FALSE
    if ( (sum(imatch)>0) & doit ) {
      lat[i] <- gdcn.inv$lat[imatch]
      lon[i] <- gdcn.inv$lon[imatch]
      alt[i] <- gdcn.inv$alt[imatch]
      loc[i] <- gdcn.inv$location[imatch]
      stnr[i] <- attr(X,"Station_number")
      cntr[i] <- attr(X,"Country_code")
      nval[i] <- sum(is.finite(X))
      start[i] <- min(attr(X,'year'))
      end[i] <- max(attr(X,'year'))
      
      ## TEST:
      ## Use the temperature data with a clear annual cycle to mak sure that the
      ## right matrix transpose is used...
      ## X <- readGDCN('GDCN/42500012386.dly',type='TMAX')
      
      ## Match the date of X with the time frame
      #startdates <- paste(attr(X,'year'),attr(X,'month'),'01',sep='-')
      #next.mon <- seq(as.Date(startdates[length(startdates)]), length=2, by='1 month')[2] 
      #last.date <- seq(next.mon, length=2, by='-1 day')[2] 
      #period <- seq(as.Date(startdates[1]),last.date,by=1)
      #print(range(period))
      
      ## The data is stored as if all months have 31 days... Needs to take this into account
      mysrt <- order(rep(attr(X,'year'),31))
      cdates <- paste(rep(attr(X,'year'),31)[mysrt],rep(attr(X,'month'),31)[mysrt],1:31,sep='-')
      valdate <- rep(FALSE,length(cdates))
      for (itd in 1:length(cdates)) {
        testdate <- try(as.Date(cdates[itd]))
        if (inherits(testdate,"try-error")) valdate[itd] <- FALSE else valdate[itd] <- TRUE
      }
      period <- as.Date(cdates[valdate])
      
      ## Remove the 'fake dates' (day 31 in months with 28, 29 or 30 days only)
      print(sum(valdate))
      x <- c(t(X)); x <- x[valdate]
      if (length(x) != length(period)) print('Warning - something strange happened!')
      x[x >=999] <- NA
      ok <- is.finite(x)
      if (sum(ok)==0) print('Strange results') # browser()
      x <- x[ok]; period <- period[ok]
      #plot(zoo(x,order.by=tim))
      
      iP <- is.element(tim,period)
      ix <- is.element(period,tim)
      if ( (sum(iP)==0) | (sum(ix)==0) | sum(iP) != sum(ix) ) print('Strange results') # browser()
      
      ## Sum the precipitation
      print(rowSums(X))
      totP[iP] <- totP[iP] + x[ix] 
      if (sum(!is.finite(totP))>0) print('Strange results') # browser()
      ## Detect rainy days
      x[x < x0] <- 0; x[x >= x0] <- 1
      fw[iP] <- fw[iP] + x[ix]
      if (sum(!is.finite(fw))>0) print('Strange results') # browser()
      ## Keep track of number of data points
      x[is.finite(x)] <- 1
      Ng[iP] <- Ng[iP] + x[ix]
      if (sum(!is.finite(Ng))>0) print('Strange results') # browser()
      
      ## Plot the points to show progress.
      points(lon[i],lat[i],pch=19,cex=0.6,col=rgb(0.5,0,0,0.2))
      print(paste(i,length(list),stationlist[i],cntr[i],stnr[i],lat[i],
                  lon[i],alt[i],sum(is.finite(X))))
      i <- i + 1
    } else if (length(dim(X))==2) 
      print(c(ii,i,N,range(attr(X,'year')),sum(is.finite(X[attr(X,'year') %in% year(tim),])),N.min)) else
      print(c(ii,i,N,range(attr(X,'year')),N.min))
    
    if (i%%300==0) {
      # Intermediate save:
      i1 <- ii+1
      save(file="readGDCN.temp.rda",totP,fw,Ng,i1,i,lon,lat,alt,stnr,cntr,nval)
    }
  } # End-of-loop: station files
  if (file.exists("readGDCN.temp.rda")) file.remove("readGDCN.temp.rda")
  Y <- cbind(round(100*fw/Ng,2),Ng,totP,fw)
  print(dim(Y))
  colnames(Y) <- c('fraction','number.total','total.Precip','number.Precip')
  raingaugedata <- zoo(Y,order.by=tim)
  attr(raingaugedata,"longitude") <- lon
  attr(raingaugedata,"latitude") <- lat
  attr(raingaugedata,"altitude") <- alt
  attr(raingaugedata,"station_id") <- stnr
  attr(raingaugedata,"country") <- cntr
  attr(raingaugedata,"n.valid") <- nval
  attr(raingaugedata,"start") <- start
  attr(raingaugedata,"end") <- end
  attr(raingaugedata,"description") <- "raingauges (GCDN)"
  attr(raingaugedata,"history") <- match.call()
  attr(raingaugedata,"x0") <- x0
  attr(raingaugedata,"source") <- "GDCN"
  #Save the main results:
  save(file="raingaugedata.rda",raingaugedata)
  
  invisible(raingaugedata)
}


## Function to plot time series and trend
timeseries <- function(x,is=1,denominator='total area',ylab=NULL,...) {
  par(bty='n')
  if (is.character(denominator)) {xd <- attr(x,denominator)/100; if (is.null(ylab)) ylab=paste(denominator,'(%)')} else
                           {xd <- denominator; if (is.null(ylab)) ylab=unit(x)[is]}
  y <- x[,is]/xd
  plot(y,main=attr(x,'longname')[is],plot.type='single',col=c('red','blue','grey'),
       sub=paste(attr(x,'src'),paste(abs(range(lat(x))),collapse='S-'),'N: ',min(year(y)),
                 ' - ',max(year(y)),sep=''),
       ylab=ylab,xlab="",...) 
  for (i in 1:length(is)) {
    lines(trend(y[,i]))
    t1 <- trend(y[,i])[1]; t2 <- trend(y[,i])[length(index(y))]
    lines(c(as.Date('1990-01-01'),index(y)[1]),rep(t1,2),col=rgb(0.5,0.5,0.5,0.5),lwd=2)
    lines(c(as.Date('1990-01-01'),index(y)[length(index(y))]),rep(t2,2),col=rgb(0.5,0.5,0.5,0.5),lwd=2)
    text(as.Date('1997-08-01'),t1,round(t1),cex=0.7)
    text(as.Date('1997-08-01'),t2,round(t2),cex=0.7)
  }
  grid()
}

## Function to plot annual cycle
AC <- function(x,plot=TRUE,ylab=NULL,...) {
  y <- as.monthly(as.station(x),FUN='mean')
  if (!plot) return(y)
  if (is.null(ylab)) ylab <- attr(x,'unit')
  z <- aggregate(y,month,FUN='mean')
  plot(z,col=c('red','blue','grey'),lwd=3,errorbar=FALSE,
       ylab=ylab,xlab='month',new=FALSE,...)
  grid()
  invisible(z)
}

## Function to plot histogram/pdf
PDF <- function(x,is=1,plot=TRUE,...) {
  z <- coredata(trend(x[,is],result = "residual"))
  h <- hist(z,col='grey',freq=FALSE,main='',ylab='',
       xlab=paste(attr(x,'longname')[is],' (',attr(x,'unit')[is],')',sep=''),plot=plot,...)
  x <- seq(-max(abs(z)),max(abs(z)),length=100)
  if (plot) {
    lines(x,dnorm(x,mean=mean(z),sd=sd(z)),lwd=3,col=rgb(1,0,0,0.4))
    grid()
    par(new=TRUE,fig=c(0.75,0.98,0.75,0.98),mar=rep(0,4),cex.axis=0.7)
    qqnorm(z,main=''); qqline(z,col='red')
  }
  invisible(list(h,x))
}



## Illustration - show map
Fig1 <- function(path='~/TRMM',
                 pattern='3B42_Daily',
                 param='precipitation',
                 nx=1440,ny=400,a = 6378, x0 =1) {
  require(ncdf4)  
  fname <- list.files(path=path,pattern=pattern,full.names=TRUE)[1]
  print(fname)
 
  ncid <- nc_open(fname)
  x <- ncvar_get(ncid,param)
  x <- t(x)
  x[x < x0] <- NA; x[x >= x0] <- 1
  dim(x) <- c(1,nx*ny)
  tim <- ncatt_get(ncid,varid = 0,attname='EndDate')
  time <- tim$value
  lon <- ncvar_get(ncid,'lon')
  lat <- ncvar_get(ncid,'lat')
  nc_close(ncid)
  attr(x,'longitude') <- lon
  attr(x,'latitude') <- lat
  x <- as.field(zoo(x,order.by=time),lon=lon,lat=lat,param='precip',unit='NA',src='TRMM')
  map(x,type='fill',colbar=FALSE)
  invisible(x)
}


Fig2 <- function() {
  data("Parea")
  timeseries(Parea)
}

## Show the time evolution of the precipitation area
Fig3 <- function(col=c('black',rgb(0.2,0.2,0.2,0.4),rgb(0.2,0.2,0.6,0.4),rgb(0.1,0.1,0.5,0.5))) {
  data("Parea")
  data("Parea.month")
  data("Parea.eraint.month")
  data("tpa.eraint")
  par(bty='n')
  x <- merge(aggregate(100*Parea[,1]/attr(Parea,'total area'),year,FUN='mean'),
             aggregate(100*Parea.month/attr(Parea.month,'total area'),year,FUN='mean'),
             aggregate(100*Parea.eraint.month/attr(Parea.eraint.month,'total area'),year,FUN='mean'),
             aggregate(50*(tpa.eraint[,2]+tpa.eraint[,3])/tpa.eraint[,5],year,FUN='mean'),all=TRUE)
  plot(x,xlab='',ylab='%',main='Rainfall area',
       plot.type='single',col=col,lwd=c(3,2,2,3))
  n <- dim(x)[2]
  stats <- rep('',n)
  for (i in 1:n) {
    lines(trend(x[,i]),lty=2,col=col[i])
    stats[i] <- paste(round(trend.coef(x[,i]),1),'%/decade (p-val=',
                      round(trend.pval(coredata(x[,i])),3),')',sep='') 
  }
  lines(x[,1],lwd=3)
  legend(1980,22.6,paste(c('TRMM day','TRMM month','ERAINT month','ERAINT day'),stats),
         col=col,lwd=c(3,2,2,3),bty='n',cex=0.8)
  grid()
}

## Similar as Fig2 but for CMIP5 simulations
Fig4 <- function() {
  require(esd)
  data("Parea.cmip5")
  coredata(Parea.cmip5) <- 100*t(t(coredata(Parea.cmip5)/attr(Parea.cmip5,'total area')))
  ## Use the model resolution to set the colour scheme for each model
  data("IPCC.AR5.Table.9.A.1")
  rbind(IPCC.AR5.Table.9.A.1$Model.Name,as.character(IPCC.AR5.Table.9.A.1$Horizontal.Grid))
  nm1 <- tolower(gsub('.','',gsub('-','',IPCC.AR5.Table.9.A.1$Model.Name),fixed=TRUE))
  nm2 <- tolower(gsub('.','',gsub('-','',attr(Parea.cmip5,'model_id')),fixed=TRUE))
  x <- cmipgcmresolution()[charmatch(nm2,nm1)]
  ## Colour for the non-specified GCMs.
  col <- rep(rgb(0,0,0,0.2),length(x))
  dx <- as.numeric(rownames(table(x))); n <- length(dx)
  pal <- rgb(1-(1:n)/n,0.1,(1:n)/n,0.3)
  for ( i in 1:n ) {col[is.element(x,dx[i])] <- pal[i]}
  
  par(bty='n')
  plot(aggregate(Parea.cmip5,year,FUN='mean'),main='Simulated monthly mean precipitation area',ylab='%',plot.type='single',
       sub='RCP4.5, CMIP5',col=col,xlab='')
  grid()
  legend(1870,36,as.character(round(dx,1)),col=pal,horiz = TRUE,pch=19,bty='n')
  
  change <- apply(aggregate(Parea.cmip5,year,FUN='mean'),2,trend.coef)
  print('Trend estimates for Fig4 in % over 2000-2100')
  print(summary(as.numeric(10*change))); print(length(change))
}


