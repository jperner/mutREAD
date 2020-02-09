##
##
##

#########################################
## LIBRARIES AND PARAMETERS
#########################################

# libraries
library(VariantAnnotation)
library(SomaticSignatures)
library(BSgenome.Hsapiens.UCSC.hg19)
library(nnls)
library(parallel)

# example vcf file
vcffile <- "data/test.vcf"
fileID <- "test"

# read predefined signatures
MUTSIGN <- read.delim("data/secrier_signatures129.txt")

# filter mutations on these chromosomes
filterchrom <- c("chrM", "hs37d5", "NC_007605")

# output 
rrs_files <- "data/regions/RestrictionEnzymeDigest_v45_newEnzSelection/"

# regions for exomes
exome_file <- "data/regions/nexterarapidcapture_exome_targetedregions_v1.2.bed"

# expanded exome data
expexome_file <- "data/regions/nexterarapidcapture_expandedexome_targetedregions.bed"

#########################################
## READ AND FILTER MUTATIONS
#########################################

# read VCF
mut <- readVcf(vcffile, "hg19")

# filter some chromosomes
mut <- mut[!seqnames(mut)%in%filterchrom]
seqlevels(mut) <- setdiff(seqlevels(mut), filterchrom)

# remove multi-alleles (e.g. G,T)
filtermuts <- elementNROWS(rowRanges(mut)$ALT)>1
mut <- mut[!filtermuts,]

# make simplified GenomicRange
snv <- as(mut, "VRanges")
snv$FileID <- rep(fileID, length(snv))

#########################################
## CALCULATE EXPOSURES BASED ON ALL SNVS
#########################################

# function to calculate exposure from snvs
fitExposures <- function(cursnv){
	 
	# calculate frequencies of trinucleotide mutations
	mutcontext <- mutationContext(cursnv, BSgenome.Hsapiens.UCSC.hg19)
	motifmat <- motifMatrix(mutcontext, group = "FileID", normalize = TRUE)

	# match trinucleotide names
	row.names(motifmat) <- sapply(strsplit(row.names(motifmat), ""), function(x) paste0(x[4],"[",x[1],">",x[2],"]",x[6]))

	# fit exposures and scale to sum 1
	if(ncol(motifmat)>0){
		exposures <- matrix(nnls(as.matrix(MUTSIGN), motifmat[row.names(MUTSIGN),])$x, nrow=1, ncol=ncol(MUTSIGN))
		if( rowSums(exposures)>0 ) exposures[1,] <- exposures[1,]/sum(exposures[1,])
	}
	row.names(exposures) <- colnames(motifmat)
	colnames(exposures) <- colnames(MUTSIGN)
	return(exposures)
}

exposures_wgs <- fitExposures(snv)

#########################################
## SIMULATION NR. OF MUTATIONS VS COSSIM
#########################################

# cosine similarity function
cosineSim <- function(a,b) crossprod(a,b)/sqrt(crossprod(a)*crossprod(b))

# for each number of mutations
cossim_list <- list()
for(n in c(10, 25)){ #, 50, 75, 100, 250, 500, 750, 1000, 2500, 5000, 7500, 10000, 25000, 50000)){
	if(length(snv)>n){
		# make a 100 repititions
		sample_idx <- sapply(1:100, function(i) sample(1:length(snv), n, replace=FALSE))
                sigmat <- apply(sample_idx, 2, function(sampind) fitExposures(snv[sampind]))
                sample_rep_list <- apply(sigmat, 2, function(sigs) cosineSim(exposures_wgs[1,], sigs))
	}else{
		sample_rep_list <- NA
	}
        cossim_list[[as.character(n)]] <- sample_rep_list
}

#########################################
## READ SIMULATED REGIONS FOR EACH SEQUENCING METHOD
#########################################

# load exome regions
readExomes <- function(file){
        # bed-like file
        exome <- read.table(file)
        names(exome) <- c("chr", "start", "end")

        # make GRange object for easy overlapping
        exome <- makeGRangesFromDataFrame(exome)

        return(exome)
}

exome <- readExomes(exome_file)
expexome <- readExomes(expexome_file)

# read the insilico digestion fasta
readIDRanges <- function(file){

        # read fasta file
        tmp <- readLines(file)
        tmp <- grep(">", tmp, value=TRUE)
        tmp <- strsplit(tmp, " \\| |: ")

        if(length(tmp)<1) return(NULL)

        # format to bed
        bed <- data.frame(      chr=sapply(tmp, "[", 14),
                                start=as.numeric(sapply(tmp, "[", 10)),
                                end=as.numeric(sapply(tmp, "[", 12)),
                                strand=sapply(tmp, "[", 8))

        # switch for fragments on negative strand; GRange prerequisite start < end
        s <- bed$start
        s[bed$strand=="-"] <- bed$end[bed$strand=="-"]
        bed$end[bed$strand=="-"] <- bed$start[bed$strand=="-"]
        bed$start <- s
        bed$strand <- rep("*", nrow(bed))

        # make GRange for easy overlapping
        bed <- makeGRangesFromDataFrame(bed)

        # simulate reads from each side (let's say 150bp)
        lread <- resize(bed, width=150, fix="start")
        rread <- resize(bed, width=150, fix="end")
        bed <- c(lread, rread)

        return(bed)
}

# read and convert all fragment files
allrrsfiles <- list.files(rrs_files, "*.fasta", full.name=TRUE, recursive=TRUE)
library(parallel)
beds <- mclapply(allrrsfiles, readIDRanges, mc.cores=3)
names(beds) <- gsub("hg19_|_400-500_v45|\\.fasta", "", basename(allrrsfiles))

# remove combinations without fragments
beds <- beds[sapply(beds, length)>0]

#########################################
## SIMULATE SEQUENCING METHODS VS COSSIM
#########################################

# function to subselect mutations that overlap the regions
fitExposuresSubRegion <- function(cursnv, regions){
	# select overlapping regions
	subsnv <- cursnv[overlapsAny(cursnv, regions),]

	# tmp
	exposures <- matrix(0, nrow=1, ncol=ncol(MUTSIGN), dimnames= list(unique(cursnv$FileID), colnames(MUTSIGN)))
	
	# fit exposures with subset
	if(length(subsnv)){
		etmp <- fitExposures(subsnv)
		exposures[1, colnames(MUTSIGN)] <- etmp[1, colnames(MUTSIGN)]
	}
	return(exposures)
}

# fit
exposures_exome <- fitExposuresSubRegion(snv, exome)
exposures_expexome <- fitExposuresSubRegion(snv, expexome)
rrs_exposures <- mclapply(beds, function(x) fitExposuresSubRegion(snv, x), mc.cores=3)

# cosine similarity to original
exome_cossim <- cosineSim(exposures_exome[1,], exposures_wgs[1,])
expexome_cossim <- cosineSim(exposures_expexome[1,], exposures_wgs[1,])
rrs_cossim <- mclapply(rrs_exposures, function(emat) cosineSim(emat[1,], exposures_wgs[1,]), mc.cores=3)

#########################################
## BASE-PAIRS COVERED SIMULATION
#########################################

# get the number of base pairs sequenced
sumbps <- sapply(beds, function(x) sum(width(x)))
bps_exome <- sum(width(exome))
bps_expexome <- sum(width(expexome))

#########################################
## BASE-PAIRS COVERED SIMULATION
#########################################

# cut the genome in 1MB bins
bins <- tileGenome(seqinfo(Hsapiens), tilewidth=1e6, cut.last.tile.in.chrom=TRUE)

# count the overlapping loci per bin
exome_locpermb <- table(factor(queryHits(findOverlaps(bins, exome)), levels=1:length(bins)))
expexome_locpermb <- table(factor(queryHits(findOverlaps(bins, expexome)), levels=1:length(bins)))
rrs_locpermb <- lapply(beds, function(x) table(factor(queryHits(findOverlaps(bins, x)), levels=1:length(bins))))


