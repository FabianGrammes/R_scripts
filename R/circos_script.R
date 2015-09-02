#--------------------------------------------------------------------
# 1) Ideogram data

# Columns: chr - ID LABEL START END COLOR
# ID = Identifier to link the data to the other files
# LABEL = Appears in the plot

#--------------------------------------------------------------------
# 2) Linkage data

# Columns: ID.1 START.1 END.1 ID.2 START.2 END.2 COLOR

#====================================================================

# converts R colors to circos frendly colors
circos.color <- function(x, alpha = NULL){
    if( is.null(alpha) ){
        col <- col2rgb( x )
    }else{
        col <- col2rgb( x, alpha=TRUE )
        col[4,] <- alpha
    }
    col <- apply( col, 2, paste, collapse ="," )
    col <- paste("color=(", col, ")", sep ="")
    return( col )
}

# This is somethong I can do later
#x <- scan("circos_source/circos.conf", what="", sep="\n", comment.char = "#")

# handlink karyogram data
get.karyogram <- function( file ){
    kary <-  read.table(file = file,    # update path
                        header = FALSE, stringsAsFactors = FALSE, sep = "\t")
    colnames(kary) <- c('chr','blank', 'label','id','start','end','color')
    return( kary )
}
update.karyogram <- function( kary, path = NULL ){
    if( is.null(path) ){
        path <- file.path( getwd(), "circos" )
    }
    if( !file.exists(path) ){
        stop( "Run setup circos first!")
    }
    write.table( file = file.path(path, "cigene3p6_chrom.karyotype"), sep="\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
}


setup.circos <- function( path = NULL, links , configuration = "std"){
    if(is.null(links)){
        print("Links data has to be given")
    }
    if( is.null(path) ){
        path <- file.path( getwd(), "circos" )
    }
    dir.create( path )
    if( configuration == "std" ){
        file.copy( "Circos_config.txt", path ) # fix file location
    }else{
        file.copy( configuration, path )
    }
    
    write.table(links, file = file.path(path, "Circos_links"), sep ="\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Do this optimization later
run.circos <- function( path = NULL ){
    if( is.null(path) ){
        path <- file.path( getwd(), "circos" )
    }
    if( !file.exists(path) ){
        stop("Can't find circos directory")
    }
    system("module load circos") # load the module
    system("circos -conf circos.conf 2>&1 > report.txt")
}


