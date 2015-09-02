# GO GOstats fisher test
go.fun <- function(gsc, selected, universe, file = "none",
                   pval = 0.05, set.min = 5, set.max = 1000,
                   conditional = TRUE){
    require( GOstats )
      Gparams =  GSEAGOHyperGParams(name="GO_BP",
          geneSetCollection = gsc,
          geneIds = selected,
          ontology = "BP",
          universeGeneIds = universe,
          pvalueCutoff = pval,
          conditional= conditional,
          testDirection = "over")
      Gparams = hyperGTest(Gparams)
      out = summary(Gparams)
      out <- out[out$Count >= set.min & out$Size <= set.max ,]
      if( file == "none"){
          return( out )
      }else{
          write.csv( out , file = file, row.names=FALSE,
                    quote=FALSE)
          return( out )
      }
  }

# KEGG GOstats fisher test
kegg.fun <- function(gsc, selected, universe, file = "none",
                   pval = 0.05, set.min = 3, set.max = 1000){
    require( GOstats )
    require( CigSsa.db )
    Kparams =  GSEAKEGGHyperGParams(name="KEGG",
          geneSetCollection = gsc,
          geneIds = selected,
          universeGeneIds = universe,
          pvalueCutoff = pval,
          testDirection = "over")
    Kparams = hyperGTest(Kparams)
    out = summary(Kparams)
    out$Term <- NULL
    term <- get.KEGG( out$KEGGID, mode = "path" )
    out <- merge( out, term, by.x = 'KEGGID', by.y='k_path', all.x=TRUE )
    out <- out[ order(out$Pvalue ), ]
    out <- out[out$Count >= set.min & out$Size <= set.max ,]
    if( file == "none"){
        return( out )
    }else{
        write.csv( out , file = file, row.names=FALSE,
                    quote=FALSE)
        return( out )
    }
  }

list.genes <- function(x){
    df <- get.genes(x)
    return(df)
}

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

