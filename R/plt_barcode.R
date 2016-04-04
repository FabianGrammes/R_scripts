library(limma)
library(grid)
library(RColorBrewer)
library(Ssa.RefSeq.db)


## Function:
## mat = matrix that should be converted to ranks (t-statistics og lfc)
## gene.df = data frame with a column gene_id that contains the relevant genes
## importnat link column : set the column name that links mat and gene.df meanning the link.col in gene.df should contain row.names(mat)..

my.barcode <- function(mat, gene.df, name.col = 'gene', link.col = 'gene_id',
                       col = NULL, palette = 'Set1', parse.label = FALSE,
                       main = ""){

    ## Internal FUNs

    ##----- fixing color: repeats color paletts
    .generate.col <- function(n, palette, col.n = 8  ){
        if(n <= 9){
            col = brewer.pal( n = n,  palette)
        }else{
            col = list()
            for( i in 1:5 ){
                col[[i]] <- brewer.pal( n = col.n,  palette)
            }
            col <- unlist(col)
        }
        return(col)   
    }

    ## ------ plotting the bars via grid
    .grid.bar <- function(x, col = 'black', yscale){
        if(length(col) == 1)
            col <- rep(col, length(x))
        for(i in 1:length(x)){
            xi <- x[[i]]$data
            N <- length(xi)
            x.bar <- as.vector( matrix( nrow = 2, c((1:N-0.25),(1:N+0.25)), byrow= T) )
            y.bar <- rep(xi, each = 2)
            ## plot the bars
            grid.polyline(x = unit(x.bar, units = 'native'),
                          y = unit(y.bar, units = 'native'),
                          id = rep(1:N, each = 2),
                          gp=gpar(col = col[i], lwd=2.5))
            ## plot the dotted connection lines
            if(N > 1){
                grid.polyline(x = unit(x.bar[2:(length(x.bar)-1)], units = 'native'),
                              y = unit(y.bar[2:(length(y.bar)-1)], units = 'native'),
                              id = rep(1:(N-1), each = 2),
                              gp=gpar(col = col[i], lwd=1.5, lty= 3))
        }
        }
    }
    
    ## rank the matrix
    rmat <- apply(mat, 2, rank)
    yscale <- c(0, ceiling(max(rmat)))
    xscale <- c(0, ncol(rmat))
    N  <- c(1:ncol(rmat))
    
    ## Check
    #if(! all(gene.df$gene_id %in% row.names(rmat))){
    #    warning('Removed some gene_ids')
    #}
    ## subset to gene set and make a list
    gene.df <- gene.df[ gene.df[, link.col] %in% row.names(rmat) ,]   
    sub <- rmat[ row.names(rmat) %in% gene.df[,link.col], ]
    ## order by last column in the rank matrix, makes plotting connections way easier
    if(N ==1){
        sub <- as.matrix( sub[order(sub)] )
    }else{
        sub <- sub[ order(sub[,length(N)]), ]
    }
    ## enforce same order
    gene.df <- gene.df[ match( row.names(sub),gene.df[,link.col]),]

    nr <- nrow(sub)
    
    if( ! all( row.names(sub) == gene.df[,link.col] ))
        stop('ERROR: gene order does not match')
    sub.l <- list()
    for (i in row.names(sub)){
        sub.l[[i]] <- list()
        sub.l[[i]][['data']] <- sub[ row.names(sub) == i,] 
        sub.l[[i]][['name']] <- gene.df[,name.col][ gene.df[,link.col] == i]
    }

    if(is.null(col))
        col <- .generate.col( nr, palette )
    
    ## PLOTTING
    grid.newpage()
    lay <- grid.layout(nrow = 3, ncol = 4,
                       width = c(0.08, 0.55, 0.45, 0.01),
                       height = c(0.1,0.8,0.1))
    pushViewport( viewport( layout = lay) )
    ##-------------------------------------------------------------------------
    ## -- PLT: Bars
    pushViewport( viewport(xscale = xscale+0.5,
                           yscale = yscale ,
                           layout.pos.row = 2, layout.pos.col = 2))
    grid.rect(gp = gpar( col = 'grey80' ))
    grid.grill(v = unit(N, units = 'native'),
               h = NA,
               gp = gpar(lty = 3, col = 'grey'))
    .grid.bar(x = sub.l, col = col, yscale)
    ## ---- Y-axis
    grid.yaxis(at = c(yscale),
               label = c('Down','Up'),
               edits = gEdit(gPath="labels", hjust=0.9),
               gp = gpar(cex = 0.8))
    ## ----- X-axis
    grid.xaxis(at = c(1:ncol(rmat)),
               label = colnames(mat),
               main = TRUE,
               edits = gEdit(gPath="labels", rot = 90,vjust=0.5, hjust = 0.7),
               gp = gpar(cex = 0.6))
    upViewport()
    
    ## -- PLT: Gene Lables
    if( parse.label ){
        gene.df[, name.col ] <- sub( " transcript variant.*","", gene.df[, name.col ])
        gene.df[, name.col ] <- sub( "%2C","", gene.df[, name.col ])
    }
    label.pos <- seq(yscale[1] ,yscale[2], length.out = nrow(sub))
    pushViewport( viewport(xscale = c(0,1),
                           yscale = yscale ,
                           layout.pos.row = 2, layout.pos.col = 3))
    grid.text(label = gene.df[, name.col ], x = 0.2,
              y = unit(label.pos, units = 'native'),
              just = 'left', gp = gpar(cex = 0.7))

    ## ---- Connect lables to bars
    label.x <- rep(c(0, 0.18), length(sub.l))
    label.y <- as.vector( matrix( nrow = 2,
                                 c(unlist(lapply(sub.l, function(x){
                                     x[['data']][length(x[['data']])]}
                                                 )), label.pos),
                                 byrow = T))
    grid.polyline(x = label.x,
                  y = unit(label.y, units = 'native'),
                  id = rep(1:length(sub.l), each =2),
                  gp = gpar( col = col, lty = 2))     
    upViewport()
    ## -- PLT: title
    pushViewport( viewport( layout.pos.row = 1, layout.pos.col = 2))
    grid.text( label = main, gp = gpar( cex =2))
    upViewport()
}



##------------------------------------------------------------------------------



# my.barcode(hk.Tf$t,  gene.df = tt[[2]], name.col = 'product', palette = 'Accent',
#           parse.label = T, main = 'Th1')






