
# FUN to converts R colors to circos frendly colors
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

# Function to setup the circos plot data:
# 1) Path to the folder with the circos data
# 2) Linkage data:
#    Columns: ID.1 START.1 END.1 ID.2 START.2 END.2 COLOR

setup.circos <- function( path = NULL, links ){
    if(is.null(links)){
        print("Links data has to be given")
    }
    if( is.null(path) ){
        print("Need path of the default data")
    }  
    dir.create( "circos_plot" )
    plot.path <- file.path( getwd(), "circos_plot")
    # copy karyogram + configuration to new location
    file.copy( file.path(path, "circos.conf") , plot.path ) # fix file location
    file.copy( file.path(path, "Circos_karyogram.txt"), plot.path )
    # write links data    
    write.table(links, file = file.path(path, "Circos_links.txt"), sep ="\t",
                row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Wrapper to run circos:
run.circos <- function( path = NULL ){
    if( is.null(path) ){
        path <- file.path( getwd(), "circos_plot" )
    }
    if( !file.exists(path) ){
        stop("Can't find circos_plot directory")
    }
    #system("module load circos") # load the module
    system("circos -conf circos.conf 2>&1 > report.txt")
}


