library(grid)
library(reshape)
library(RColorBrewer)

#=====================================================================
# EXAMPLE:
# grid.newpage()
# a = cbind(rnorm(50,3), rnorm(50,6), rnorm(50,-3), rnorm(50,2))
# plot_all(N = 4, diag_nr = c(100, 200, 300, 400), diag_names = c("a", "b", "c","d"),
#          n_matr = matrix(1:16, ncol=4), fc.matr = a, p.scale = c(-12,12), draw.dots = 'no',
#          pos = data.frame(a = c(0.1,0.5,0.9),
#                           b = c("a", "b", "c"), stringsAsFactors = FALSE))


# INPUT:

# N = ncol = nrow
# diag = names and values to be plotted diagonal
# n_matr = matrix ov Nr of overlapping genes (order must be same as diag)
# fc.matr = matrix of logFC
# p.scale = lims of the plotting area




#=====================================================================
# Functions

# plot function for the diagonal
grid_diag <- function(N, diag_names, diag_nr){
    for(i in 1:N){
        pushViewport( viewport(layout.pos.col= i+1, layout.pos.row= i+1))
        grid.roundrect(gp=gpar(fill="grey", alpha=0.6))
        grid.text(diag_names[i] , y = 0.65,  gp=gpar(cex=1.5))
        grid.text(diag_nr[i] , y = 0.3,  gp=gpar(cex=2, fontface='italic'))
        upViewport()
    }
}
# plot text in the lower triangles
grid_low <- function( lo.lay, n.matr, col.matr = NULL){
    for( i in 1:dim(lo.lay)[1]){
        row = lo.lay[i,1]
        col = lo.lay[i,2]
        pushViewport( viewport(layout.pos.col= col+1, layout.pos.row= row+1))
        if( is.null(col.matr) ){
            grid.roundrect()
        }else{
            grid.rect(gp=gpar(fill=col.matr[col,row], alpha=0.6))
        }
        grid.text(n.matr[col,row],  gp=gpar(cex=2.5) )
        upViewport()
    }
}
# function to bin the data and plot..

plot_bins <- function(x, y, p.scale, bins){
    colPlate <- colorRampPalette(c('lightsteelblue1 ','navy'))
                                        # cut the vectors
    x = cut(x, seq(p.scale[1], p.scale[2], length.out = bins))
    y = cut(y, seq(p.scale[1], p.scale[2], length.out = bins))
                                        # 2D matrix
    mat = as.data.frame.matrix(table(x, y))
                                        # order according to 2D coordinate system
    mat = mat[rev(1:bins),]
    dims = dim(mat)
                                        # logaritmic cutting points and color
                                        # breaks = c(0, exp(seq(0,7, length.out = bins-1)))
    breaks = seq(0, 24, length.out = bins)^2
    thecol = c("white", colPlate(bins-1))
                                        # transform to color  
    mat.col = thecol[cut(as.matrix(mat), breaks = breaks)]
                                        # retransform to matrix
    mat.col <- matrix(mat.col, nrow=dims[1], ncol=dims[2])
                                        # plot
    grid.raster(mat.col, interpolate=TRUE, width=1, height=1)  
}

# plot the points in the upper triangle
grid_upper <- function( up.lay, fc.matr, p.scale, draw.dots ){
    for( i in 1:dim(up.lay)[1]){
        row = up.lay[i,1]
        col = up.lay[i,2]
                                        # calculate corelation
        cor.r = round(cor(x=fc.matr[,col], y=fc.matr[,row]), digits =2)
        pushViewport( dataViewport( xscale=p.scale, yscale=p.scale,
                                   layout.pos.col= col+1, layout.pos.row= row+1))
                                        #grid.roundrect()
        if( draw.dots == 'yes' ){
            grid.points(x=fc.matr[,col], y=fc.matr[,row], pch=16,
                        size= unit(0.3, "char"), gp=gpar(col = "blue", alpha=0.4) )
        }
        if( draw.dots == 'no' ){
            plot_bins(x=fc.matr[,row], y=fc.matr[,col], p.scale, bins =50)
        }
        grid.roundrect()
        grid.lines(x = unit(c(0, 1), "npc"),  y = unit(c(0, 1), "npc"), gp=gpar(col ="red", lty=3))
        grid.lines(x = unit(c(0, 1), "npc"),  y = unit(c(0.5, 0.5), "npc"),
                   gp=gpar(col ="black", lty=2))
        grid.lines(x = unit(c(0.5, 0.5), "npc"),  y = unit(c(0, 1), "npc"),
                   gp=gpar(col ="black", lty=2))
        grid.text(paste("r = ",cor.r), x=0.75, y=0.1, gp=gpar(cex=1))
        upViewport()
    }
}

# plot the axis for the dot plots
grid_ax <- function(N, pos){
    for(i in 1:(N-1) ){
                                        # x-axis
        pushViewport( viewport(layout.pos.col= i+2, layout.pos.row= 2))
        grid.xaxis(at = pos[,1] , label = pos[,2],edits = gEdit(gPath="labels", vjust = 1.2, hjust=0.5),
                   main =FALSE, gp = gpar(cex = 0.75))
        upViewport()
                                        # y-axis
        pushViewport( viewport(layout.pos.col= N+1, layout.pos.row= i+1))
        grid.yaxis(at = pos[,1] , label = pos[,2], edits = gEdit(gPath="labels", rot=90, vjust = 0.5, hjust=0.5),
                   main = FALSE, gp = gpar(cex = 0.75))
        upViewport()
    }
}
# alltogether
plot_all <- function(N, diag_names, diag_nr, n_matr, fc.matr, p.scale, draw.dots = 'no',
                     pos, col.matr = NULL  ){
    lay.mat <- matrix(1, nrow=N, ncol=N)
    # upper triangel
    up.lay <- melt(upper.tri(lay.mat))
    up.lay <- up.lay[up.lay[,3] ==TRUE,]
    # lower triangle
    lo.lay <- melt(lower.tri(lay.mat))
    lo.lay <- lo.lay[lo.lay[,3] ==TRUE,]
    # page layout
    my.layout <- grid.layout( N+2, N+2, width=c(0.01,rep(1/N,N),0.05),
                             height=c(0.05 ,rep(1/N,N), 0.01))
    pushViewport(viewport(layout = my.layout))
    #grid.show.layout(my.layout)
    # plot diagonal
    grid_diag( N, diag_names, diag_nr )
    # plot lower triangle
    if( is.null(col.matr ) ){
        grid_low( lo.lay, n_matr)
    }else{
        grid_low( lo.lay, n_matr, col.matr )
    }
    
    if(! is.null(fc.matr)){
        grid_upper( up.lay, fc.matr, p.scale, draw.dots )
        grid_ax( N, pos )
    }    
}




