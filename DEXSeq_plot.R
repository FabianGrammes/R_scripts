
#===============================================================================
# Fabian Grammes 2015
#===============================================================================

# - plot.DEX function with extended functionality compared to DEXSeq:::plotDEXSeq
# - The main plotting the same, though a bit more fancy :)
# - Main difference is that it's possible to include inflrmation from a external
#   .gff3 file in the plotting of the transcripts. IT'S IMPORTANT to match the
#   information of the .gff3 file here and the DEXSeq generated .gff file.

# INPUT:
# --object: A DEXSeqResults object
# --geneID
# --fitExpToVar: The column with the group ID (same as in DEXSeq:::plotDEXSeq)
# --expression: TRUE/FALSE havn't done enything with it yet
# --FDR: Cutoff for the FDR (default = 0.01)
# --gff: .gff3 file as data.table; just use fread to read in gff, nothing more.

#===============================================================================


plot.DEX <- function(object, geneID, fitExpToVar, gff, expression = TRUE, FDR = 0.01, geneNAME = NULL ){

    # LOAD required packaged
    require(reshape)
    require(DEXSeq)
    require(grid)
    require(data.table)
    
    if( is(object, "DEXSeqResults")){
        rt<-which(object$groupID==geneID)
        sampleData <- as.data.frame( object@sampleData )
        genomicData <- object$genomicData
        each <- object$padj[rt]
        if(any(is.na(each))){
            each[is.na(each)] <- 1
        }
    }
    # ERASE CHUNK ?
    #if(sum(count) == 0){
    #    warning("No read counts falling in this gene, there is nothing to plot.")
    #    return()
    #}

    # - ARRANGE DATA -

    # -- ARRANGE THE DATA FOR THE EXON MODEL TRACK --
    if( length( start(unlist(genomicData))) > 0 ){
        model.track <- data.frame(start=start(genomicData[rt]),
                                  end=end(genomicData[rt]),
                                  chr=as.character( seqnames( genomicData[rt] ) ),
                                  strand=as.character( strand( genomicData[rt] ) ) )
        model.scale <- c(min(as.vector(model.track[,c('start','end')])),
                         max(as.vector(model.track[,c('start','end')])))
    }

    model.track$ticks <- apply(model.track[,1:2], 1, mean)

    # -- ARRANGE DATA FOR CDS TRACK --

    # ---- FUNCTION to import and clean the GFF data (data.table) -----
    # ===> CAN BE CUSTOMIZED TO DIFFERENT FORMATS <====
    clean.gff <- function( gff, id.var = 'Parent=', geneID ){
        sub <- gff[ grepl(geneID, gff$V9), ]
        if( nrow(sub) == 0 )
            stop('Filtering went completly wrong!!')
        sub[, gsub('=','',id.var) := sapply(sub$V9, ff, id.var = id.var, USE.NAMES = FALSE ) ]
        # remove non-relevant rows
        sub[,c('V2','V6','V8','V9'):=NULL]
        sub <- as.data.frame(sub)
        colnames(sub) <- c('chr', 'type', 'start', 'end', 'strand','parent')
        sub <- sub[!sub$type %in% c('mRNA','gene'),]
        sub <- sub[order(sub$parent, sub$start),]
        sub <- split( sub, sub$parent )
        sub <- sub[ order( as.numeric(gsub('\\S+\\.t', '', names(sub))) ) ]
        return(sub)
        # ---- Function for internal string filtering ----
        ff <- function(x, id.var){
            sp <- unlist(strsplit(x, split =';'))
            sp <- gsub(paste(id.var), '', sp[grepl('Parent', sp)])
            if(identical(sp, character(0)))
                sp <- NA
            return(sp)
        }
    }
    # --- FILTER THE GFF DATA ---
    tgff <- clean.gff( gff = gff, geneID = geneID)

    # -- ARANGE DATA FOR THE EXPRESSION TRACK --
    # ----- COPIED partly FROM DEXSeq: Part 3 = Calc normalized counts
    count <- t( t(object$countData[rt,])/sampleData$sizeFactor )
    count <- DEXSeq:::vst( count, object )
    count <- melt(count)
    count$exon <- as.factor(gsub("\\S+:", "", count$X1))
    count <- merge(count, sampleData[,c('sample.name',fitExpToVar)],
                   by.x = 'X2', by.y = 'sample.name', all.x = TRUE,
                   sort = FALSE)
    dim(count)
    
    numexons <- length(each)
    # ----- COPIED FROM DEXSeq: Part 1 = arrange the data for calculation of coefficients
    if( expression ){
        stopifnot(is( object, "DEXSeqResults"))
        mf <- object@modelFrameBM
        mf <- mf[as.vector( sapply( split( seq_len(nrow(mf)), mf$sample ), "[", seq_len( numexons ) ) ),]
        featuresInGene <- object$featureID[rt]
        mf$exon <- factor( rep( featuresInGene, nrow(sampleData) ) )
        counts <- object$countData[rt,]
        rownames(counts) <- gsub("\\S+:", "", rownames(counts))
        dispersions <- object$dispersion[rt]
        dispersions[is.na( dispersions )] <- 1e-8
        names(dispersions) <- object$featureID[rt]    
        for( i in seq_len(nrow(mf))){
            mf[i,"dispersion"] <- dispersions[as.character(mf[i,"exon"])]
            mf[i,"count"] <- counts[as.character(mf[i,"exon"]), as.character(mf[i,"sample"])]
        }
        mf <- droplevels( mf )  
    }
    # ----- COPIED FROM DEXSeq: Part 2 = Calc coefficients
    if(expression){  
        es <- DEXSeq:::fitAndArrangeCoefs(frm=as.formula(paste("count ~", fitExpToVar, "* exon")),
                                          balanceExons=TRUE,
                                          mf=mf)
    if(is.null(es)){
        warning(sprintf("glm fit failed for gene %s", geneID))
        return()
    }
        coeff <- as.matrix( t(
            DEXSeq:::getEffectsForPlotting(es, averageOutExpression=FALSE, groupingVar=fitExpToVar) )[featuresInGene,] )
        coeff <- exp(coeff)
        ylimn <- c(0, max(coeff, na.rm=TRUE))
        coeff <- DEXSeq:::vst( coeff, object )
        coeff <- melt(coeff)
        coeff <- split(coeff, coeff$group)
    }

    # - SETUP PLOT FUNCTIONS -

    # --- function to plot single exons ---
    plot.exons <- function(x, fill, col = NA){
        if(length(col) == 1 )
            col <- rep(col, nrow(x))
        if(length(fill) == 1)
            fill <- rep(fill, nrow(x))
        for(i in 1:nrow(x)){
            xx <- unlist(x[i,c('start','end'), drop = TRUE])
            xx <- rep(xx, each = 2)
            if(i > 1 ){
                gxl <- c(x[i-1,'end'], x[i,'start'], x[i,'start'])
                gxl[2] <- mean(gxl[1:2])
                grid.polyline(x = unit(gxl, units = 'native'),
                              y = c(0.5,0.7,0.5),
                              gp = gpar(col = 'black', lwd = 1.5))
            }
        }
        for(i in 1:nrow(x)){
            xx <- unlist(x[i,c('start','end'), drop = TRUE])
            xx <- rep(xx, each = 2)
            grid.polygon(x = unit(xx, units = 'native'),
                         y = c(0.1,0.9,0.9,0.1),
                         gp = gpar(fill = fill[i], col = col[i], lwd = 0.5))
        }
    }

    # --- function to loop the plot.exon for all transcripts ---
    loop.exons <- function(x, xscale){
        exon.lay <- grid.layout(nrow = length(x), ncol = 1,)
        pushViewport( viewport( layout = exon.lay))
        for(i in 1:length(x)){
            pushViewport( viewport(layout.pos.row = i, layout.pos.col = 1,
                                   xscale = xscale ) )
            fill <-  c('grey30', 'grey80', 'grey80')[as.factor(x[[i]]$type)]
            plot.exons(x[[i]], fill = fill, col = fill)
            upViewport()
        }
        popViewport()
    }

    # --- function to make a colorRampPalette with alpha
    addAlpha <- function(colors, alpha) {
        paste(colors, sprintf("%x", ceiling(255*alpha)), sep="")
    }

    # --- function to plot the GLM means as bars
    grid.bar <- function(x, col = 'black'){
        if(length(col) == 1)
            col <- rep(col, length(x))
        for(i in 1:length(x)){
            xi <- x[[i]]
            N <- nrow(xi)
            x.bar <- as.vector( matrix( nrow = 2, c((1:N-0.25),(1:N+0.25)), byrow= T) )
            y.bar <- rep(xi$value, each = 2)
            grid.polyline(x = unit(x.bar, units = 'native'),
                          y = unit(y.bar, units = 'native'),
                          id = rep(1:N, each = 2),
                          gp=gpar(col = col[i], lwd=2))
            # plot the dotted connection lines
            grid.polyline(x = unit(x.bar[2:(length(x.bar)-1)], units = 'native'),
                          y = unit(y.bar[2:(length(y.bar)-1)], units = 'native'),
                          id = rep(1:(N-1), each = 2),
                          gp=gpar(col = col[i], lwd=1.5, lty= 3))
        }
    }

    # - LAYOUT -
    # eventuall use a function to determine the height of the GFF track
    calc.size <- function( x ){
        if( x > 0.5 )
            x <-  0.5
        return(x)       
    }
    cds.size <- calc.size( length(tgff)/ 25)
    
    grid.newpage()  # ==> ERASE LATER <==
    lay <- grid.layout(nrow = 6, ncol = 3,
                       width = c(0.07, 0.8, 0.2),
                       height = c(0.1,0.4,0.1,0.08,cds.size, 0.05))
    pushViewport( viewport( layout = lay))

    # - SET PLOTTING COLORS -
    # -- DE color --
    model.fill <- c('white','firebrick1')[as.numeric(each <= FDR)+1]
    col.grid <- c('grey70','firebrick1')[as.numeric(each <= FDR)+1]
    model.bars <- c('#00FF7F', colorRampPalette(c('dodgerblue', 'deeppink1'))(3) )
    model.points <- addAlpha(model.bars, alpha = 0.5)
    model.points <- model.points[as.factor(count$group)]

    # set grid line type
    lty = 3
    
    # - PLOT THE DIFFERENT TRACKS -
    # -- PLOT EXPRESSION TRACK --
    # how many exons
    nlev <- length(levels(coeff[[1]]$exon))
    # get scales for the plot
    exp.scale.x <- c(0.5, length(levels(count$exon))+0.5)
    exp.scale.y <- c(0, ceiling(max(count$value)))
    exp.track.ticks <- c(1:nlev)
    pushViewport( viewport(xscale = exp.scale.x,
                           yscale = exp.scale.y,
                           layout.pos.row = 2, layout.pos.col = 2 ) )
    # GRILL
    grid.grill(v = unit(exp.track.ticks, units = 'native'),
               h = NA,
               gp = gpar(lty = lty, col = col.grid))
    # Y-axis
    grid.yaxis(at = seq(0,exp.scale.y[2], length.out = 5),
               label = seq(0,exp.scale.y[2], length.out = 5),
               edits = gEdit(gPath="labels", hjust=0.6),
               gp = gpar(cex = 0.8))
    # X-axis on top
    grid.xaxis(at = c(1:nlev),
               label = levels(coeff[[1]]$exon),
               main = FALSE,
               edits = gEdit(gPath="labels", rot = 90,vjust=0.5),
               gp = gpar(cex = 0.8))
              
    grid.points(x = jitter(as.numeric(count$exon)), y = count$value,
                size = unit(0.1, units = 'native'),
                pch = 19,
                gp = gpar( col = model.points))

    grid.bar(coeff, col = model.bars)
    # get relative x-coordinates 
    exp.x.npc <- convertX( unit(exp.track.ticks, 'native'), 'npc' )
    upViewport()


    # -- PLOT THE EXON MODEL TRACK --
    pushViewport( viewport(xscale = model.scale,
                           layout.pos.row = 4, layout.pos.col = 2 ) )
    # get relative x-coordinates 
    model.x.npc <- convertX( unit(model.track$ticks, 'native'), 'npc' )
    grid.grill(v = unit(model.x.npc, units = 'npc'),
               h = NA,
               gp = gpar(lty = lty, col = col.grid))
    plot.exons(model.track, fill = model.fill, col = 'black')
    upViewport()

    # --- TRACK LABEL ---
    pushViewport( viewport(layout.pos.row = 4, layout.pos.col = 3 ) )
    grid.text(label = 'Exon model', x = 0.05, just = 'left')
    upViewport()

    # -- PLOT CONNECTION LINES --
    pushViewport( viewport(layout.pos.row = 3, layout.pos.col = 2 ) )
    axis.conn <- as.vector( matrix( c(exp.x.npc, model.x.npc), nrow = 2,
                                   byrow = TRUE))
    grid.polyline(x = unit(axis.conn, units = 'npc'),
                  y = rep(c(1,0), length(axis.conn)/2),
                  id = rep(c(1:nlev), each =2),
                  gp = gpar(lty = lty, col = col.grid))
    upViewport()

    # -- PLOT GRILL LINES --
    pushViewport( viewport(layout.pos.row = 5, layout.pos.col = 2 ) )
    grid.grill(v = unit(model.x.npc, units = 'npc'),
               h = NA,
               gp = gpar(lty = lty, col = col.grid))
    upViewport()

   # -- PLOT THE CDS TRACK --
   pushViewport( viewport(layout.pos.row = 5, layout.pos.col = 2 ) )
   loop.exons(tgff, xscale = model.scale)
   tiks <- round(seq(model.scale[1], model.scale[2], length.out = 5))
   upViewport()

   # --- CDS AXIS ---
   pushViewport( viewport(layout.pos.row = 5, layout.pos.col = 2, xscale = model.scale ) )
   grid.xaxis(at = tiks, label = tiks,
              edits = gEdit(gPath="labels",vjust= -0.5),
              gp = gpar(cex = 0.8))
   upViewport()
              
   # --- TRACK LABEL ---
   pushViewport( viewport(layout.pos.row = 5, layout.pos.col = 3 ) )
   grid.text(label = names(tgff),
             x = 0.05,
             y = seq(1-1/length(tgff)/2, (1/length(tgff)/2), length.out = length(tgff)),
             just = 'left',
             gp = gpar(cex = 0.8))
   upViewport()

   # -- PLOT LEGEND
   pushViewport( viewport(layout.pos.row = 1:2, layout.pos.col = 3 ) )
   grid.text(label = paste(geneNAME, geneID, sep ='\n'),
             x = 0.5, y = 0.85,
             gp = gpar(cex = 1.2, fontface = 'bold'))
   grid.polyline(x = rep(c(0.15,0.2),length(coeff)),
                 y = rep(seq(0.6, 0.2, length.out = length(coeff)), each = 2),
                 id = rep(1:length(coeff), each =2),
                 gp = gpar(lwd =5, col = model.bars))
   grid.text(label = names(coeff),
             x = 0.25,
             y = seq(0.6, 0.2, length.out = length(coeff)),
             just = 'left',
             gp = gpar(cex = 1))
   upViewport()
   # - END -
}
