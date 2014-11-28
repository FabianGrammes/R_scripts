 library(grid)

# function for normalization and scaling of the counts
ns.counts <- function(mat, col =c('dodgerblue2','white','red') ){   
    my.col <- colorRampPalette(col)
    # scale + cluster rows
    scaled <- t(scale(t(mat)))
    dis <- dist(scaled)
    dis <- as.dendrogram(hclust(dis))
    or <- order.dendrogram(dis)
    # reaorder the color map
    scaled <- scaled[or,]
    # color
    col.X <- my.col(100)[cut(scaled, 100)]
    col.X <- matrix(col.X, nrow=nrow(scaled),
                    ncol=ncol(scaled),
                    dimnames = list( rownames(scaled), colnames(scaled)))
    colnames(col.X) <- colnames(mat)
    return(list(col=col.X, ord=or))
}

# split the heatmap
s.split <- function(X, group){  
    XL <- list()
    for(i in seq(levels(group))){ 
        XL[[i]] <- X[,group == levels(group)[i]]
        names(XL)[i] <-  levels(group)[i]
    }
    return(XL)
}

# Grid.plotting functions
plot.heat <- function(x, r, c, rect=TRUE){
    pushViewport( viewport(layout.pos.col= c, layout.pos.row= r))
    grid.raster(x, interpolate=FALSE, width=1, height=1)
    if( rect ){
        grid.rect(gp=gpar(col="black", fill="transparent"))
    }
    upViewport()  
}
# plot the header
plot.text <- function(text, r, c, ...){
    pushViewport( viewport(layout.pos.col= c, layout.pos.row= r))
    grid.text(text, ...)
    upViewport()   
}
# plot x axis
plot.xax <- function(x, r, c, ...){
    lab <- colnames(x)
    ats <- seq(lab)/length(lab)-0.5/length(lab)
    pushViewport( viewport(layout.pos.col= c, layout.pos.row= r))
    grid.xaxis(at= ats, label= lab, edits = gEdit(gPath="labels", rot=90, hjust=0.6), ... )
    upViewport() 
}
# plot y axis
plot.yax <- function(r.names,r,c, ...){
    pushViewport( viewport(layout.pos.col= c, layout.pos.row= r))
    grid.yaxis(at = seq((1/length(r.names))/2,1-(1/length(r.names)/2),
                   length.out=length(r.names)),
               label = r.names, main= FALSE, gp = gpar(cex = 0.7))
    upViewport()
}
# Final heatmap function (all previous functions joined together)
grid.heat <- function(XL, r.names=NULL, xax.gpar ){
    nWidth <- unlist(lapply(XL, ncol))
    nWidth <- nWidth/sum(nWidth)
    width <- c()
    for(i in seq(nWidth)){
        w <- c(0.05, nWidth[i])
        width <- c(width, w)
    }
    if(is.null(r.names)){
        width <- c(width, 0.05)
    }else{
        width <- c(width, 0.35)
    }
    height <- c(0.075, 0.8, 0.1)
    h.layout <- grid.layout(nrow = length(height),
                            ncol = length(width),
                            widths = width,
                            heights = height)
    pushViewport(viewport(layout = h.layout))
    for(i in seq(XL)){
        plot.heat(XL[[i]], r = 2, c = i*2)
        plot.text(names(XL)[i], r = 1, c = i*2)
        plot.xax(XL[[i]], r =2, c = i*2, gp = xax.gpar)      
    }
    if(!is.null(r.names)){
       plot.yax(rev(r.names), c = length(width)-1, r=2) 
    }   
}


# Heatmap function for plotting 2 heatmaps over each other
grid.heat2 <- function(XL1, XL2, r.names1=NULL, r.names2=NULL, h.names=NULL, xax.gpar ){
    nWidth <- unlist(lapply(XL1, ncol))
    nWidth <- nWidth/sum(nWidth)
    width <- c()
    for(i in seq(nWidth)){
        w <- c(0.02, nWidth[i])
        width <- c(width, w)
    }
    if(is.null(r.names1)){
        width <- c(width, 0.05)
    }else{
        width <- c(width, 0.2)
    }
    height <- c(0.075, 0.2, 0.005, 0.075, 0.65, 0.005)
    h.layout <- grid.layout(nrow = length(height),
                            ncol = length(width),
                            widths = width,
                            heights = height)
    pushViewport(viewport(layout = h.layout))
    # plot XL1
    for(i in seq(XL1)){
        plot.heat(XL1[[i]], r = 2, c = i*2, rect = FALSE)
        plot.text("a", r = 1, c = 1, y = 0.8, x = 1, gp = gpar(cex = 1.9))
        plot.text(names(XL1)[i], r = 1, c = i*2, y = 0.2, x = 0.5)
        plot.text(h.names[1], r = 1, c = 2:6, y = 0.6, x = 0.5,
                  gp = gpar(fontace = "bold", cex = 1.7) )
        #plot.xax(XL[[i]], r =2, c = i*2, gp = xax.gpar)      
    }
    if(!is.null(r.names1)){
       plot.yax(rev(r.names1), c = length(width)-1, r=2) 
    }
    # plot XL2
    for(i in seq(XL2)){
        plot.heat(XL2[[i]], r = 5, c = i*2, rect =FALSE)
        plot.text("b", r = 4, c = 1, y = 0.8, x = 1, gp = gpar(cex = 1.9))
        plot.text(h.names[2], r = 4, c = 2:6, y = 0.6, x = 0.5,
                  gp = gpar(fontace = "bold", cex = 1.7) )
        plot.text(names(XL2)[i], r = 4, c = i*2, y = 0.2, x = 0.5)
        #plot.xax(XL[[i]], r =2, c = i*2, gp = xax.gpar)      
    }
    if(!is.null(r.names2)){
       plot.yax(rev(r.names2), c = length(width)-1, r=5) 
    }
}
