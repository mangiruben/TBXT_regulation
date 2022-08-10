#version2; creating bed files from significant fragments/reads to be used for merging, and visualisation of spider plots
#written by Reuben June,2019 
#condition path to Coutput file.. name=cell_line)
#reuben.yaa@kcl.ac.uk ; Wardle lab
#########################################


library(optparse)

option_list <- list(
make_option(c("--condition"),help="put the path of the Coutput file", type="character",metavar="condition"),
make_option(c("--name"),help="name of the condition",type="character",metavar="name"),
make_option(c("--chrom"),help="chrom number of the viewpoint",type="character",metavar="chrom"),
make_option(c("--outdir"),help="Output file, will be created in the working directory",type="character",metavar="OUTDIR",default= "./")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$condition))
{
	stop("You haven't specified the condition file path")
}

True_Int <-function (condition, name, chrom, outdir) 
{ 
	peak_table <-read.table(opt$condition, header=T)
	peak_table2<- peak_table[ 1:(nrow(peak_table)-1),]
	peak_table2$end <-peak_table$start[2:nrow(peak_table)]
	peak_table2$chrom <-opt$chrom
	Final<- peak_table2[peak_table2$peak == "TRUE",c("chrom","start","end", "peak")]
	True_peak <- file.path(opt$outdir, paste0(opt$name, "-True_peak.bed")) 
	write.table(Final, file=True_peak, sep ="\t" , col.names= FALSE, row.names = FALSE, quote =FALSE)

}

True_Int()
