

pfam.con <- function(){
                                        # get the data subdirectory 
    db.wd <- system.file("extdata", "PFAM.sqlite", package = "PFAM.db") 
                                        # set up the connection
    drv <- dbDriver("SQLite")
    con.pf <- dbConnect(drv,dbname = db.wd)
    return(con.pf)
}


get.PFAM.DE <- function(x){
    x = paste(x,sep=', ', collapse="', '")
    x = paste("('", x, "')", sep="")
    dbGetQuery(con.pf, paste("SELECT *
FROM de
WHERE ac IN",x, sep =" "))
}
