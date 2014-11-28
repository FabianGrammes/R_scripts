# Functions for GO enrichment using the GOstats package


# GOoverrepresentation function
my.go.over <- function(gsc, genes, universe, pval, onto, cutoff = 5, conditional = TRUE){
      Gparams =  GSEAGOHyperGParams(name="GO_BP",
          geneSetCollection = gsc,
          geneIds = genes,
          ontology = onto,
          universeGeneIds = universe,
          pvalueCutoff = pval,
          conditional= conditional,
          testDirection = "over")
      Gparams = hyperGTest(Gparams)
      out = summary(Gparams)
      out <- out[out$Count >= cutoff,]
      return(out)
  }

list.genes <- function(x){
    df <- get.genes(x)
    return(df)
}

# Function to list the significant genes by GO term as data frame
# express Data should be a data.frame with all data that should be
# merged with the GOs.
# merged by express.data$gene_id

# note th
get.go.genes <- function( gsc, GOs, significant, express.data=NULL ){
    require( plyr )
    go.list <- geneIds( gsc[ GOs ] )
    row.list <- lapply( go.list, function(x) x[ x %in% significant ] )
    gene.list <- lapply( row.list, list.genes ) 
    gene.list <- ldply( gene.list, data.frame )
    if( !is.null(express.data) ){
         out <- merge( gene.list, express.data, by = 'gene_id', all.x = TRUE)
         out <- out[ order( out$.id, out$gene_name ), ]
         return( out )
    }else{
        return( gene.list )
    }
}

p.barplot <- function( x, label, axe.label, col ){
    require(grid)
    require(gridBase)
    oldpar <- par()
    par(mar=c(2,0.5,0,1))
    plt <- barplot(log2(x), col=col, names.arg="", horiz = T, axes =F)
    axis(1,at=log2(rev(axe.label)),labels=format( rev( axe.label ), scientific=TRUE))
    ## Use grid to add the labels    
    vps <- baseViewports()
    pushViewport(vps$inner, vps$figure, vps$plot)
    #grid.rect( gp = gpar( col = "blue"))
    grid.text(label, x = unit(log2(x), "native"),
              y = unit(plt, "native"),
              just="right", hjust = -0.01)
    popViewport()
    par <- oldpar
}
