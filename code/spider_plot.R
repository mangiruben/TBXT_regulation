#version2; creating spider plots around the viewpoint(Single chromosome)
#corrects splines drwaing to +-1MB to the viewpoint#
#Adopted from Splinter 2012 with adjustments by Reuben 

#reuben.yaa@kcl.ac.uk ; Wardle lab 2019
#Runs on the terminal >>> Rscript spider_plot.R
#########################################
library(optparse)


option_list <- list(
make_option(c("--dom"),help="path to matrix or data.frame must have 2 columns showing, start and end position of the interactions", type="character",metavar="dom"),
make_option(c("--vp.loc"),help="coordinate of the viewpoint",type="numeric",metavar="vp.loc"),
make_option(c("--outdir"),help="path for output file, default writing in the working directory",type="character",metavar="OUTDIR",default= "./spider_plots"),
make_option(c("--name"),help="name of the condition or test",type="character",metavar="NAME")

)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$dom))
{
  stop("You haven't specified the data frame")
}

dir.create(file.path(opt$outdir,opt$name), recursive=T, showWarnings=F)

drawChrom <- function(chrom.len =50e6, max.wid = 100e6, hei=5, chrom.wid=0, y.loc=5)
{
  r <- seq(0,2*pi, len=1000)
  chrom.wid = chrom.wid/2

  dim.val <- par("din")
  left  <- r[501:1000]
  right <- r[1:500]
  
  correction <- ((max.wid*chrom.wid)/(2*hei)) / (dim.val[1]/dim.val[2])

  x <- sin(left)*correction + correction
  y <- cos(left)*chrom.wid+y.loc

  x1 <- x
  y1 <- y

  x <- c(x, sin(right)*correction+chrom.len -correction)
  y <- c(y, cos(right)*chrom.wid+y.loc )

  x2 <- sin(right)*correction+chrom.len -correction
  y2 <- cos(right)*chrom.wid+y.loc
  polygon(x,y, col='white')
  
  col.seq <- seq(0.5,1,len=250)
  col.seq <- c(col.seq,rev(col.seq))
  cols <- rgb(col.seq,col.seq,col.seq)
  segments(x1,y1,x2,rev(y2), col = cols)
  polygon(x,y, lwd=2)


}

drawLocalChrom <- function(labels=c("Xa","Xi"), yloc1 = 3, yloc2 = 5, wid = 2, num.chrom=1, vp.loc, outdir, OUTPUT)
{
  chrom.start <-(opt$vp.loc-1000000)
  chrom.len <-  (opt$vp.loc+1000000)
  plot(c(chrom.start, chrom.len), c(yloc1-wid,yloc2+wid), type='n', axes=F, xlab="", ylab="", cex.lab=1)
  
  mb <- floor(chrom.len/1e6)
  mb10 <- floor(mb/10)*10
  #draw an axis
  label <- c(seq(0, mb10, by=10),mb)
  #segments(label*1e6, -1e9, label*1e6, 1e9, lwd=2, lty=2, col='grey90')
  mid.y <- (yloc1+yloc2)/2
  #segments(label*1e6, mid.y-0.15, label*1e6, mid.y+0.15, lwd=2)
  #segments(0, mid.y, mb*1e6, mid.y, lwd=2)
  
  
  #active X
  if(num.chrom==2){
    drawChrom(chrom.len=chrom.len, max.wid=chrom.len, y.loc=yloc1)
    text(-1e6,yloc1, labels[2],cex=1)
  }	
  #inactive X 
  drawChrom(chrom.len=chrom.len, max.wid=chrom.len, y.loc=yloc2)
  
  text(-1e6,yloc2, labels[1],cex=1)
  at <- 0:mb
  #axis(1, at=at*1e6, labels=NA, cex.axis=1, lwd=1,las=1)
  #axis(1, at=label*1e6, labels=label, cex.axis=1, lwd=2,las=1)
  #axis(3, at=at*1e6, labels=NA, cex.axis=1, lwd=1,las=1)
  #axis(3, at=label*1e6, labels=label, cex.axis=1, lwd=2)

}
drawSplines.domain <- function( dom, vp.loc, y.base=5.5, y.arc=8, plot=F, col='blue', chrom.size=166e6, relative=F, OUTPUT)
{
  dom<-read.delim(opt$dom, header=FALSE)
  if(nrow(dom) == 0)
    return
  for(i in 1:nrow(dom))
  {
    #start <- min(dom[i,1], dom[i,2])but this depends on which columns your interactions are
    #end <- max(dom[i,1], dom[i,2]) but this depends on which columns your interactions are
    start <- dom[i,2]
    end <- dom[i,3]
    xspline(c(opt$vp.loc,(end + opt$vp.loc)/2, start, end, (end + opt$vp.loc)/2, opt$vp.loc), c(y.base,y.arc,y.base,y.base,y.arc,y.base), open=F, shape=c(0,1,0,0,1,0), col=col, border=col)
  }
}
Spider<- function(name, outdir){
a=file.path(opt$outdir,opt$name, paste0(opt$name," spline.png"))
png(file =a,2000,800)
plot.new()
drawChrom()
par(new=TRUE)
drawLocalChrom()
drawSplines.domain()

}
Spider()


