#-------------------------------------------------------------------------
# convenient selection of colors based on HCL color space
#-------------------------------------------------------------------------
zcolors <- function(n, parameter = c("std","rain","snow","temperature","terrain"), type=c("qualitative","sequential","diverge","warm","cold"), plot = FALSE, ...){
    
     stopifnot( is.numeric(n) )
     parameter <- match.arg(parameter)
     type <- match.arg(type)
     
     pal <- define_colors(n, parameter, type, ...)
     
     if(plot)
          image(matrix(1:n),col=pal,main=paste(parameter,type),
                axes=FALSE,xlab="",ylab="",srt=45)
     
     return(pal)
}



show_palettes <- function(n = 10, parameter = c("std","rain","snow","temperature","terrain"), type = c("qualitative","sequential","diverge","warm","cold")){

     palettes <- list()
     for(p in parameter){
          for(t in type){
               pal <- define_colors(n, p, t)
               palettes[[paste0(p," ",t)]] <- pal
          }
     }
     
     rows <- length(names(palettes)) 
     
     # bottom, left, top, right
     par(mfrow=c(rows,1), mar = c(0.1,0.4,2,0.5), oma=c(0,0,2,0))
     for(i in 1:rows){
          image(x=matrix(1:n),y=0.1,col=palettes[[i]],
                main=names(palettes)[i],
                axes=FALSE,xlab="",ylab="",srt=45)
     }
     mtext("ZAMG color palettes", line=0, side=3, outer=TRUE)
}


define_colors <- function(n, p, t, ...){
     pal <- NULL
     switch(p,
            std = switch(t,
                         qualitative = { pal <- qualitative_hcl(n,"Dark3", ...) } ,
                         sequential = { pal <- sequential_hcl(n, "Purple-Yellow", p1 = 1.3, c2 = 20, rev = TRUE, ...) } ,
                         diverge = { pal <- diverge_hcl(n,"Cyan-Magenta", ...) } 
            ),
            rain = switch(t,
                          sequential = { pal <- sequential_hcl(10, palette = "Purple-Yellow", rev = TRUE, c1 = 70, cmax = 100, l2 = 80, h2 = 500) } ,
                          diverge = { pal <- diverge_hcl(n, h = c(260, 0), c = 100, l = c(30, 90), power = 0.7, ...) } 
            ),
            snow = switch(t,
                          sequential = { pal <- sequential_hcl(n, "Blues") } ,
                          diverge = { pal <- diverge_hcl(n, h = c(260, 0), l = c(30, 90), power = 0.7, ...) } 
            ),
            temperature = switch(t,
                                 sequential = { pal <- sequential_hcl(n, 'Purple-Yellow', p1 = 1.3, c2 = 20, ...) } ,
                                 diverge = { pal <- diverge_hcl(n, h = c(255, 12), c = c(50, 80), l = c(20, 97), power = c(1, 1.3), ...) } ,
                                 warm = { pal <- sequential_hcl(n, h = 10, c = c(65, 100, NA), l = c(20, 97), power = 1.1, rev = TRUE, ...) } ,
                                 cold = { pal <- sequential_hcl(n, h = 245, c = c(50, 75, NA), l = c(20, 98), power = 0.8, rev = TRUE, ...) } 
            ),
            terrain = switch(t,
                             sequential = { pal <- sequential_hcl(n, "Terrain", ...) } 
                             
            )
     )
     return(pal)

}


