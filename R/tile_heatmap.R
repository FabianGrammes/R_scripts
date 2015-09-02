#=====================================================================
# Script to make a tile heatmap = diffreent tiles by sample and gene
# group .

# use like this:
# col.matrix <- ns.counts(dat)
# tile.heat(col.matrix, col.list, row.list)

# row.list: list of gene ids with gene sybmbols as names
#=====================================================================

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
    col.X <- matrix(col.X, nrow=nrow(scaled), ncol=ncol(scaled))
    colnames(col.X) <- colnames(mat)
    row.names(col.X) <- row.names(mat)
    return(list(col=col.X, ord=or))
}
# heatmap plotting
plot.heat <- function(x, r, c){
    pushViewport( viewport(layout.pos.col= c, layout.pos.row= r))
    grid.raster(x, interpolate=FALSE, width=1, height=1)
    grid.rect(gp=gpar(col="black", fill="transparent"))
    upViewport()  
}
# plot xaxis
plot.xax <- function(lab, r, c, ...){
    ats <- seq(lab)/length(lab)-0.5/length(lab)
    pushViewport( viewport(layout.pos.col= c, layout.pos.row= r))
    grid.xaxis(at= ats, label= lab, edits = gEdit(gPath="labels", rot=90) )
    upViewport() 
}
# plot yaxis 1st level (gene symbols)
plot.yax1 <- function(r.names,r,c, ...){
    pushViewport( viewport(layout.pos.col= c, layout.pos.row= r))
    grid.yaxis(at = seq((1/length(r.names))/2,1-(1/length(r.names)/2),
                   length.out=length(r.names)),
               label = r.names, main= FALSE, gp = gpar(cex = 0.7))
    upViewport()
}
grid.brace <- function(){
  grid.lines(x = c(0.15,0.15),
             y = c(0,1))
  grid.lines(x = c(0.05,0.15),
             y = c(1,1))
  grid.lines(x = c(0.05,0.15),
             y = c(0,0))
             }
mosaic.length <- function(Nlist, space){
    nLength <- unlist(lapply(Nlist, length))
    nLength <- nLength/sum(nLength)
    lenge <- c()
    for(i in seq(nLength)){
        l <- c(space, nLength[i])
        lenge <- c(lenge, l)
    }
    return(lenge)
}
tile.heat <- function(mat, col.list, row.list, space = c(0.1,0.1,0.1,0.1,0.01), cex2=1.2){
    # space = c(head, bottom, left1, left2, mosaic.space)
    # mosaic width
    width <- mosaic.length(col.list, space = space[5])
    # mosaic height
    height <- mosaic.length(row.list, space = space[5])
    # Add space for col/row.names
    width <- c(width,space[3:4])
    height <- c(height,space[2])
    height[1] <- space[1]
    # Layout
    h.layout <- grid.layout(nrow = length(height),
                            ncol = length(width),
                            widths = width,
                            heights = height)
    pushViewport(viewport(layout = h.layout))
    # plot Mosaic heatmaps
    for(c in seq(col.list)){
        for(r in seq(row.list)){
            ri = which(row.names(mat) %in% row.list[[r]])
            ci = which(colnames(mat) %in% col.list[[c]])
            plot.heat(mat[ri ,ci], r = r*2, c=c*2)
        }
    }
    # plot colnames
    for(c in seq(col.list)){
        pushViewport( viewport(layout.pos.col= c*2, layout.pos.row= 1))
        grid.text(names(col.list)[c], x=0.5, y=0.5,  gp=gpar(cex=1.2))
        upViewport()
    }
    # plot xaxis 
    for(c in seq(col.list)){
        plot.xax(col.list[[c]], r=length(row.list)*2, c = c*2)
    }
    # plot row.names level 1
    for(r in seq(row.list)){
        plot.yax1(names(row.list[[r]]), r= r*2, c = length(col.list)*2)
    }
    # plot row.names level 2
    for(r2 in seq(row.list)){
        pushViewport( viewport(layout.pos.col= length(col.list)*2+2, layout.pos.row= r2*2))
        grid.brace()
        grid.text(names(row.list)[r2], x=0.6, y=0.5,  gp=gpar(cex=cex2))
        upViewport()
    }
}
list.genes <- function(x){
    df <- get.genes(x)
    out <- df$gene_oldid
    names(out) <- df$gene_symbol
    return(out)
}

