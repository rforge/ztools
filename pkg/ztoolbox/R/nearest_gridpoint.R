#-------------------------------------------------------------------------
# finds nearest point to lon/lat on grid, which can be of several types 
#-------------------------------------------------------------------------
nearest_gridpoint <- function(lon, lat, grid, measure = c("euclidian","haversine")){
     
     stopifnot( is.data.frame(grid) | is.matrix(grid) & colnames(grid) %in% c("lon","lat") )
     cn <- colnames(grid)
     c.i.lon <- which(cn == "lon")
     c.i.lat <- which(cn == "lat")
     
     measure <- match.arg(measure)
     
     
     switch(measure,
          euclidian = {
               euclid.dist <- function(x1,x2,y1,y2){
                    sqrt((x2-x1)^2+(y2-y1)^2)
               }
               index <- which.min(euclid.dist(lon,grid[,c.i.lon],lat,grid[,c.i.lat]))      
          },
          haversine = {
               # calculates geodesic distance between two points specified by radian latitude/longitude using the
               # Haversine formula (hf)
               haversine <- function(x1,x2,y1,y2) {
                    R <- 6371 # mean earth radius [km]
                    dx <- x2 - x1
                    dy <- y2 - y1
                    a <- sin(dy/2)^2 + cos(y1) * cos(y2) * sin(dx/2)^2
                    #c <- 2 * asin(min(1,sqrt(a)))
                    c <- 2 * asin(apply(cbind(1,sqrt(a)),1,min))
                    d = R * c
                    return(d) # Distance in km
               }
               index <- which.min(haversine(lon,grid[,c.i.lon],lat,grid[,c.i.lat])) 
          }
     )
     
     if( length(index) > 1 )
          warning("there is more than one nearest gridpoint")
     
     res <- list(idx = index, lon = grid[index,c.i.lon], lat = grid[index,c.i.lat])
     return(res)
    
}





