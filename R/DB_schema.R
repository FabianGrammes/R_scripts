# Functions to dra DB diagrams in R

# The code for the functions and the whole idea of drawing a DB schema is kind of based on:
# Drawing Diagrams with R (Paul Murell, R.Journal Vol. 1/1, May 2009 
#---------------------------------------------------------------------
## Example use:
# boxG <- boxGrob(c("index","id", "probe_id", "acc"),
#                 x=0.1, y=0.5 )
# boxGS <- boxGrob(c("gene_info", "probe_id", "gene_name","gene_symbol","gene_E_Val"),
#                x=0.4, y=0.65 )

# grid.draw(boxG)
# grid.draw(boxGS)

# connect(boxG, boxGS, 2,1)
#---------------------------------------------------------------------

tableBox<-function(labels, x=0.5, y=0.5, clrs=c("gray55","white","gray90")){
        nlabel<-length(labels)
        tablevp<-viewport(x=x,y=y, width=max(stringWidth(labels))+unit(4,"mm"),height=unit(nlabel,"lines"),
                          just = c("center", "center"))
        pushViewport(tablevp)
        grid.roundrect()
        parity<-"odd"
        if ((nlabel %% 2)==0){parity<-"even"}
        if (nlabel>1){
                for (i in 1:nlabel-1){
                        if (i==(nlabel-1)){
                                fill<-clrs[1]
                        }else{
                                if ((nlabel %% 2)==0){
                                fill<-clrs[2:3][i %% 2 + 1]
                                } else {
                                fill<-clrs[3:2][i %% 2 + 1]
                                }
                        }
                        grid.clip(y=unit(i, "lines"), just="bottom")
                        grid.roundrect(gp=gpar(fill=fill))
                }
        }
        grid.clip()
        grid.text(labels[1], x=unit(2,"mm"), y=unit(nlabel - 0.5, "lines"), just="left", gp=gpar(fontface="bold"))
        grid.text(labels[2:nlabel], x=unit(2,"mm"), y=unit((nlabel-1):1 - 0.5, "lines"), just="left")
        popViewport()
}
boxGrob <- function(labels, x=.5, y=.5) { 
  grob(labels=labels, x=x, y=y, cl="box") 
} 
drawDetails.box <- function(x, ...) { 
  tableBox(x$labels, x$x, x$y) 
} 
xDetails.box <- function(x, theta) { 
  nlines <- length(x$labels) 
  height <- unit(nlines, "lines") 
  width <- unit(4, "mm") + max(stringWidth(x$labels)) 
  grobX(roundrectGrob(x=x$x, y=x$y, width=width, height=height), 
        theta) 
} 
yDetails.box <- function(x, theta) { 
  nlines <- length(x$labels) 
  height <- unit(nlines, "lines") 
  width <- unit(4, "mm") + max(stringWidth(x$labels)) 
  grobY(rectGrob(x=x$x, y=x$y, width=width, height=height), 
        theta) 
}
# draws the connections
connect<-function(box1, box2, lab1, lab2, arr=arrow(type="closed",angle=15,length=unit(2,"mm")), dir="up"){
        #lab1 is an index of item in the box1, headding has index 0
        xa<-grobX(box1, "east")
        ya<-grobY(box1, "north")-unit(lab1+0.5,"lines")
        xb<-grobX(box2, "west")
        yb<-grobY(box2, "north")-unit(lab2+0.5,"lines")
        c<-1
        if (dir!="up"){c<--1}
        grid.curve(xa,ya,xb,yb,curvature=c,inflect=TRUE,arrow=arr,gp=gpar(fill="black"))
}
