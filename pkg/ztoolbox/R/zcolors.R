#-------------------------------------------------------------------------
# convenient selection of colors based on HCL color space
#-------------------------------------------------------------------------
zcolors <- function(n, parameter = c("std","rain","snow","temperature","terrain"), type=c("qualitative","sequential","diverge"), plot = FALSE, ...){
    
     stopifnot( is.numeric(n) )
     parameter <- match.arg(parameter)
     type <- match.arg(type)
     
     colors <- define_colors(n, ...)
     pal <- colors[[parameter]][[type]]
     
     if(plot)
          image(matrix(1:n),col=pal,main=paste(parameter,type),
                axes=FALSE,xlab="",ylab="",srt=45)
     
     return(pal)
}



show_palettes <- function(n = 10){
     colors <- define_colors(n)
     i <- 1
     palettes <- list()
     for(p in names(colors)){
          for(t in names(colors[[p]])){
               pal <- colors[[p]][[t]]
               palettes[[paste0(p," ",t)]] <- pal
          }
     }
     
     rows <- length(names(palettes)) 
     
     # bottom, left, top, right
     #opar <- par()      # copy current settings
     par(mfrow=c(rows,1), mar = c(0.1,0.4,2,0.5), oma=c(0,0,2,0))
     for(i in 1:rows){
          image(x=matrix(1:n),y=0.1,col=palettes[[i]],
                main=names(palettes)[i],
                axes=FALSE,xlab="",ylab="",srt=45)
     }
     mtext("ZAMG color palettes", line=0, side=3, outer=TRUE)
     #par(opar)  
}


define_colors <- function(n, ...){
     
     colors <- list(
          std = list(
               qualitative = qualitative_hcl(n,"Dark3", ...), 
               sequential = sequential_hcl(n, "Lajolla", ...),
               diverge = diverge_hcl(n,"Cyan-Magenta", ...) 
          ),
          rain = list(
               qualitative = NULL,
               sequential = sequential_hcl(n, "Purple-Yellow", ...),
               diverge = diverge_hcl(n, h = c(260, 0), c = 100, l = c(30, 90), power = 0.7)
          ),
          snow = list(
               qualitative = NULL,
               sequential = sequential_hcl(n, "Blues", ...),
               diverge = diverge_hcl(n, h = c(260, 0), c = 100, l = c(30, 90), power = 0.7)
         ),
          temperature = list(
               qualitative = NULL,
               sequential = sequential_hcl(n, "Purple-Yellow", p1 = 1.3, c2 = 20),
               diverge = diverge_hcl(n, h = c(260, 0), c = 100, l = c(30, 90), power = 0.7)
          ),
          terrain = list(
               qualitative = NULL,
               sequential = sequential_hcl(n, "Terrain", ...),
               diverge = NULL
          )
     )
     return(colors)
}


