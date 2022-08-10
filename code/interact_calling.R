##### Boundary caller from 4Cseq experiments that uses in bulk deWitLab/peakC package (Geeven et al. 2018, https://doi.org/10.1093/nar/gky443) with some minors adjustements #####
##### Rafael Domínguez Acemel (José Luis Gómez Skarmeta lab); rdacemel@gmail.com ####


library(isotone)
library(peakC)
library(stringr)
library(optparse)
library(tools)

single_peakC <- function(frag1,frag2,chrom,coord,name,outdir,swSize,lwSize,zwSize,pval,ithr,zthr){
  message(paste0("Processing ",name))
  input1 <- write_input(frag1,chrom,name,"rep1",outdir)
  input2 <- write_input(frag2,chrom,name,"rep2",outdir)
  input <- c(input1,input2)
  PEAKC_OUTPUT <- run_peakC(input,name,coord,swSize,lwSize)
  PEAKC_OUTPUT$zeros1 <- zerosmooth(PEAKC_OUTPUT$score1,window=zwSize)
  PEAKC_OUTPUT$zeros2 <- zerosmooth(PEAKC_OUTPUT$score2,window=zwSize)
  PEAKC_OUTPUT$zeros <- PEAKC_OUTPUT$zeros1 + PEAKC_OUTPUT$zeros2
  PEAKC_OUTPUT <- selectPeaks(PEAKC_OUTPUT,pval,ithr,zthr)
  outfname <- file.path(outdir,"outputs",paste0(name,".Coutput"))
  write.table(PEAKC_OUTPUT,
              file=outfname,
              sep="\t",
              row.names = F)
  boundaries <- getBoundaries(PEAKC_OUTPUT,chrom,name)
  plotfile <- file.path(outdir,"plots",paste0(name,".png"))
  ub <- as.integer(boundaries[2])
  db <- as.integer(boundaries[3])
  plot_peakC(PEAKC_OUTPUT,plotfile,ithr,zthr,zwSize,ub,db)
  return(boundaries)
}

write_input <- function(frag,chrom,name,rep,outdir){
  FRAG <- read.table(frag,
                     col.names = c("chromosome","start","end","score"),
                     colClasses = c("character","integer","integer","integer"))
  FRAG$middle <- as.integer(floor((FRAG$start + FRAG$end)/2))
  FRAG <- subset(FRAG,
                 subset = chromosome == chrom,
                 select = c("middle","score"))
  
  input <- file.path(outdir,"inputs",paste0(name,"_",rep,".Cinput"))
  write.table(FRAG,
              file = input,
              sep = "\t",
              quote = F,
              row.names = F,
              col.names = F)
  return(input)
}

run_peakC <- function(input,name,coord,swSize,lwSize){
  EXPS <- readMultiple(input,vp.pos = coord,window = lwSize)
  R <- combined.analysis(EXPS,num.exp = 2,vp.pos =coord,wSize = swSize)
  REP1 <- as.data.frame(EXPS[[1]])
  REP2 <- as.data.frame(EXPS[[2]])
  colnames(REP1) <- c("start","score1")
  colnames(REP2) <- c("start","score2")
  PEAKC_OUTPUT <- merge(REP1,REP2,by ="start")
  SM <- R$dbR
  SM$pval <- R$p.value
  colnames(SM) <- c("start","sm_score1","sm_score2","iso1","iso2","pval")
  PEAKC_OUTPUT <- merge(PEAKC_OUTPUT,SM,by="start")
  return(PEAKC_OUTPUT)
}

zerosmooth <- function(vec,window=25){
  zeros <- c()
  for(i in 1:length(vec)){
    if ((i-window)<=0){
      ups <- vec[0:i]
      downs <- vec[i:(i+window)]
    }else if((i+window)>length(vec)){
      ups <- vec[(i-window):i]
      downs <- vec[i:length(vec)]
    }else{
      ups <- vec[(i-window):i]
      downs <- vec[i:(i+window)]
    }
    val <- max(c(sum(ups>0),sum(downs>0)))
    zeros <- c(zeros,val)
  }
  return(zeros)
}

selectPeaks <- function(OUTPUT,pval,ithr,zthr){
  miniso1 <- ithr*max(OUTPUT$iso1)
  miniso2 <- ithr*max(OUTPUT$iso2)
  zthr <- quantile(OUTPUT$zeros,probs = c(zthr))
  OUTPUT$peak <- (OUTPUT$pval < 0.1) & (OUTPUT$zeros > zthr) & (OUTPUT$iso1 > miniso1) & (OUTPUT$iso2 > miniso2)
  return(OUTPUT)
}

getBoundaries <- function(OUTPUT,chrom,name){
  pstart <- unlist(subset(OUTPUT,subset = peak,select=start))
  boundaries <- c("chrom"=chrom,"lstart"=pstart[1],"lend"=pstart[length(pstart)],"name"=name)
  return(boundaries)
}

plot_peakC <- function(OUTPUT,fname,ithr,zthr,zwSize,ub,db){
  miniso1 <- max(OUTPUT$iso1)*ithr
  miniso2 <- max(OUTPUT$iso2)*ithr
  zthr <- quantile(OUTPUT$zeros,probs = c(zthr))
  png(fname,800,800)
  par(mfrow=c(3,1))
  plot(sm_score1~start,data=OUTPUT,type="h",col=ifelse(OUTPUT$peak,"red","black"))
  lines(iso1~start,data=OUTPUT,col="blue",lwd=5)
  abline(h=miniso1,col="blue",lwd=5)
  abline(v=c(ub,db),lwd=2)
  plot(sm_score2~start,data=OUTPUT,type="h",col=ifelse(OUTPUT$peak,"red","black"))
  lines(iso2~start,data=OUTPUT,col="blue",lwd=5)
  abline(h=miniso2,col="blue",lwd=5)
  abline(v=c(ub,db),lwd=2)
  plot(zeros~start,data=OUTPUT,type="h",ylim=c(0,zwSize*2),col=ifelse(OUTPUT$peak,"red","black"))
  abline(h=zthr,col="blue",lwd=5)
  abline(v=c(ub,db),lwd=2)
  dev.off()
}


option_list <- list(
  make_option(c("--frag1"),help="bedgraph of rep1, not needed if --input is used for bulk",type="character",metavar="FRAG1"),
  make_option(c("--frag2"),help="bedgraph of rep2, not needed if --input is used for bulk",type="character",metavar="FRAG2"),
  make_option(c("--chrom"),help="Chromosome of the viewpoint, not needed if --input is used for bulk",type="character",metavar="CHROM"),
  make_option(c("--coord"),help="Coordinate of the viewpoint, not needed if --input is used for bulk",type="numeric",metavar="COORD"),
  make_option(c("--name"),help="Name of the experiment, not needed if --input is used for bulk",type="character",metavar="NAME",default="unknown"),
  make_option(c("--input"),help="TSV file for bulk analysis with one row per experiment with frag1-frag2-chrom-coord-name",type="character",metavar="INPUT"),
  make_option(c("--outdir"),help="Output folder path, will be created (default ./peakC)",type="character",metavar="OUTDIR",default="./peakC"),
  make_option(c("--swSize"),help="Smoothing window size in fragments (default=51)",type="integer",metavar="SWSIZE",default=51),
  make_option(c("--lwSize"),help="Window to look for interactions in bp (default=+-1000000)",type="integer",metavar="LWSIZE",default=1000000),
  make_option(c("--zwSize"),help="Sparsity smoothing window size, in fragments (default=25)",type="integer",metavar="ZWSIZE",default=25),
  make_option(c("--pval"),help="Threshold for peakC to call significant interactions (default=0.1)",type="integer",metavar="PVAL",default=1),
  make_option(c("--ithr"),help="Height threshold for the isotonic (default=0.025)",type="numeric",metavar="ITHR",default=0.025),
  make_option(c("--zthr"),help="Sparsity threshold, expressed in percentile (default=0.75)",type="numeric",metavar="ZTHR",default=0.75)
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

checkargs1 <- (length(opt$frag1) > 0) && (length(opt$frag2) > 0) && (length(opt$chrom) > 0) && (length(opt$coord) > 0)
checkargs2 <- length(opt$input) > 0

if (!xor(checkargs1,checkargs2)){
  stop("Not proper arguments specified")
} 

dir.create(opt$outdir,recursive=T,showWarnings = F)
dir.create(file.path(opt$outdir,"inputs"),recursive=T,showWarnings = F)
dir.create(file.path(opt$outdir,"plots"),recursive=T,showWarnings = F)
dir.create(file.path(opt$outdir,"outputs"),recursive=T,showWarnings = F)

if (checkargs1){
  boundaries <- single_peakC(opt$frag1,opt$frag2,opt$chrom,opt$coord,opt$name,opt$outdir,opt$swSize,opt$lwSize,opt$zwSize,opt$pval,opt$ithr,opt$zthr)
  BOUND <- data.frame("chrom"=boundaries[1],
                      "lstart"=boundaries[2],
                      "lend"=boundaries[3],
                      "name"=boundaries[4])
  bfile <- file.path(opt$outdir,paste0(opt$name,".boundaries"))
}


if(checkargs2){
  INPUT <- read.table(opt$input,sep="\t",
                      col.names = c("frag1","frag2","chrom","coord","name"),
                      colClasses = c("character","character","character","numeric","character"))
  bchroms <- c()
  bstarts <- c()
  bends <- c()
  bnames <- c()
  for (i in 1:nrow(INPUT)){
    frag1 <- INPUT$frag1[i]
    frag2 <- INPUT$frag2[i]
    chrom <- INPUT$chrom[i]
    coord <- INPUT$coord[i]
    name <- INPUT$name[i]
    boundaries <- single_peakC(frag1,frag2,chrom,coord,name,opt$outdir,opt$swSize,opt$lwSize,opt$zwSize,opt$pval,opt$ithr,opt$zthr)
    bchroms <- c(bchroms,boundaries[1])
    bstarts <-  c(bstarts,boundaries[2])
    bends <- c(bends,boundaries[3])
    bnames <- c(bnames,boundaries[4])
  }
  
  BOUND <- data.frame("chrom"=bchroms,
                      "lstart"=bstarts,
                      "lends"=bends,
                      "name"=bnames)
  bfile = file.path(opt$outdir,paste0(basename(file_path_sans_ext(opt$input)),".boundaries"))
}

write.table(BOUND,
            file = bfile,
            sep = "\t",
            quote = F,
            row.names = F)


