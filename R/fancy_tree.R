

# Tip lables do not match in order
# -- function to list the alignment blocks --
find.algn.interval <- function(x){
    out <- list()
    out$range <- c(1, length(x))
    xx <- ! x %in% 24:25
    out$ival <- list()
    oi <- FALSE
    idx <- 0
    for( i in 1:length(xx) ){
        if( oi == FALSE & xx[i] == TRUE ){
            idx <- idx+1
            st <- i
        }
        if( oi == TRUE & xx[i] == FALSE ){
            end <- i-1
            out$ival[[idx]] <- c(st,end)
        }
        if( xx[i] == TRUE & i == length(x) ){
            end <- i
            out$ival[[idx]] <- c(st,end)
        }   
        oi <- xx[i]
    }
    return(out)
}

# SCREWS UP THE FUNCTION:
# find.algn.intervalV <- Vectorize(find.algn.interval, SIMPLIFY = FALSE)

library(grid)
library(gridBase)

# -- function to plot alignments, only 1 --
plot.algn <- function(x, fill, col = NA, lwd = 1){
    pushViewport( viewport( width = 0.95, height = 0.75, xscale = x$range) )
    # plot alignemnt blocks
    x <- x$ival
    for(i in 1:length(x)){
        xx <- x[[i]]
        xx <- rep(xx, each = 2)
        grid.polygon(x = unit(xx, units = 'native'),
                     y = c(0,1,1,0),
                     gp = gpar(fill = fill, col = fill, lwd = 0.5))
    }
    # bottom line
    grid.lines(x = unit(c(0,1), units = 'npc'),
               y = unit(c(0,0), units = 'npc'),
               gp = gpar(col = fill, lwd = lwd))
    # top line
    grid.lines(x = unit(c(0,1), units = 'npc'),
               y = unit(c(1,1), units = 'npc'),
               gp = gpar(col = fill, lwd = lwd))
    # pop
    popViewport()                                
}

plot.algn.wrapper <- function( algn, col, lwd ){
    nrow <- length(algn)
    lay <- grid.layout(nrow = length(algn),
                       ncol =1,
                       heights = rep(1/nrow, nrow))
    ralgn <- rev(algn)
    pushViewport( viewport( layout = lay))
    for( i in 1:length(algn)){
        pushViewport( viewport(layout.pos.row = i, layout.pos.col = 1)  )
        aa <- find.algn.interval(ralgn[[i]])
        plot.algn(aa, fill = col, lwd)
        #grid.text( label = names(ralgn)[i] )
        upViewport()
    }
    popViewport()
}


# - wrapping the functions -

plot.tb <- function( tree, algn, block.col, block.lwd, nl = FALSE, ... ){
    # --- 1st check
    if( ! all(tree$tip.label %in% names(algn)) )
        stop( 'ERR1: TREE and Alignemnt do not match')
    # --- reaorder alignment according to tree ---
    ord <- match(tree$tip.label, names(algn))
    algn <- algn[ ord ]
    # --- 2nd check
    if( ! all(tree$tip.label == names(algn)) )
        stop( 'ERR2: Names in TREE and Alignemnt do not match')
    # -- plot dimensions -- scaling factor for the tree y-axis
    yfact <- 1/(0.93 + 0.93/length(algn))
    # page layout
    glay <- grid.layout(nrow = 3,
                        ncol = 4,
                        widths = c(0.05, 0.6, 0.2, 0.05),
                        heights = c(0.05, 0.9, 0.05))
    pushViewport(viewport(layout = glay))
    # --- plot Dendrogram ---
    pushViewport( viewport(layout.pos.col = 2 , layout.pos.row = 2) )
    pushViewport( viewport(height = yfact, width =1))
    par( plt = gridPLT());par(new=TRUE)
    plot(tree, cex = 0.5, ...)
    upViewport()
    upViewport()
    # --- plot alignment bars ---
    pushViewport( viewport(layout.pos.col = 3, layout.pos.row = 2 ) )
    plot.algn.wrapper(algn, col = block.col, lwd = block.lwd)
    upViewport()
}

