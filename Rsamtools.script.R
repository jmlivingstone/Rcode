library("Rsamtools")
require(multicore)

# fill in information or parse from .gff file
gene.info <- NA
gene.info$gene.name <- "ERBB2"
gene.info$chr <- "chr17"
gene.info$gene.start <- 37844167
gene.info$gene.end <- 37886679

files <- list.files()
files.2 <- dir("", full.names=T)

# index all files in directory
mclapply(1:length(files.2),function(i){
        print(i)
        output.bam.file <- files.2[i]
        index.bam(output.bam.file)
})

# Get the read counts for each file
# each patient will get an arbritary number - need to remember to match to file vector later
for(i in 1:length(files.2)){
        print(i)
        output.file <- paste("nreads/PAT_",i,".nreads",sep="")
        x <- samtools.read.count(files.2[i])
        write.table(x,output.file,row.names=FALSE,col.names=FALSE,sep=",")
}

# Prepare smaller BAM files for the region given
# make sure subbam folder exists
mclapply(1:length(files.2),function(i){
        output.bam.file <- paste("subbam/PAT_",i,".bam",sep="")
        print(output.bam.file)
        subBam(files.2[i],output.bam.file, gene.info$chr, gene.info$gene.start,gene.info$gene.end)
})

# Create the index files, and depth files for smaller region
# make sure depth and subbam folder exists
mclapply(1:length(files.2),function(i){
        print(i)
        output.bam.file <- paste("subbam/PAT_",i,".bam",sep="")
        output.depth.file  <- paste("depth/PAT_",i,".depth",sep="")
        index.bam(output.bam.file)
        genomeCoverage.depth(output.bam.file,output.depth.file)
})

N <- get.nreads(files.2)
depth.matrix <- get.depth.matrix(gene.info)

for(i in 1:length(files.2)){
        depth.matrix[,i] <-     depth.matrix[,i] / (N[i] / 100000000)
}

# only want to look at the exprs across specific bases that are isoform specific
region <- depth.matrix[4:371,]

# calculate the average expression across the region
avg.exprs <- apply(region, 2, mean)

#order the patients based on the average exprssion
o <- order(avg.exprs)

#### functions ###

# calls installed samtools program and uses already formatted hg19 genome file
genomeCoverage.depth <- function(bam.file,output.depth.file) {
        system.c <- paste("./bin/genomeCoverageBed -ibam ",bam.file," -g hg19.genome -bga -split > ",output.depth.file,sep="")
        system(system.c)

}

# create a matrix will the read counts
get.nreads <- function(files.2){
        N <- sapply(1:length(files.2),function(i) {
                nreads.file <- paste("nreads/PAT_",i,".nreads",sep="")
                read.delim(nreads.file,header=FALSE)[1,1]
        })
	return(N)
}
