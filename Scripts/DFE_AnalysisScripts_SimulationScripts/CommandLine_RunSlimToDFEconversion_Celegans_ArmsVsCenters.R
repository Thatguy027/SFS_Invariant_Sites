# the function:

make.est_dfe.input <- function(poly.dat, genome.dat, fixed.dat, generation, num.inds.sampled, genome.size, filename, fold=FALSE, use.manual.sample=FALSE){

	# because diploid:
	sample.size <- 2 * num.inds.sampled
	
	# make data frames of all possible neutral and deleterious mutations at all time points recorded
	neut.muts <- NULL
	seln.muts <- NULL
			
	## WHERE poly.dat IS A FILE OF ROWS OF MUTATIONS OCCURRING
	#	m1 = neutral site in coding
	#	m2, 3 = deleterious selected site in coding
	#	m4 = beneficial selected site in coding
		
# tack on fixed data and then can include counts for 0's
	fixed.mut.dat <- fixed.dat[fixed.dat$gen.fixed <= as.numeric(generation) ,]
		# this gives only mutations that have fixed PRIOR to and INCLUDING WITHIN the current generation time point sampled
	fixed.neut.muts <- c(which(fixed.mut.dat$mut.type == "m1"))
	fixed.seln.mut.IDs <- fixed.mut.dat$mut.ID[-fixed.neut.muts]
	fixed.neut.mut.IDs <- fixed.mut.dat$mut.ID[fixed.neut.muts]
	
	num.neut.muts.fixed <- length(fixed.neut.mut.IDs)
	num.seln.muts.fixed <- length(fixed.seln.mut.IDs)
	
	
	neut.muts <- poly.dat[poly.dat$mut.type == "m1" ,]
	seln.muts <- poly.dat[poly.dat$mut.type != "m1" ,]
			
	if(use.manual.sample == FALSE){
		sfs.total <- table(poly.dat$mut.prev)
		sfs.neut <- table(neut.muts$mut.prev)
		sfs.seln <- table(seln.muts$mut.prev)
		
		# add on fixed things to the fixed section of the table
		if(is.na(sfs.total[as.character(sample.size)])){
			sfs.total[as.character(sample.size)] <- num.neut.muts.fixed + num.seln.muts.fixed
		}else{sfs.total[as.character(sample.size)] <- sfs.total[as.character(sample.size)] + num.neut.muts.fixed + num.seln.muts.fixed}
		if(is.na(sfs.neut[as.character(sample.size)])){
			sfs.neut[as.character(sample.size)] <- num.neut.muts.fixed
		}else{sfs.neut[as.character(sample.size)] <- sfs.neut[as.character(sample.size)] + num.neut.muts.fixed}
		if(is.na(sfs.seln[as.character(sample.size)])){
			sfs.seln[as.character(sample.size)] <- num.seln.muts.fixed
		}else{sfs.seln[as.character(sample.size)] <- sfs.seln[as.character(sample.size)] + num.seln.muts.fixed}
		
		sfs.total["0"] <- genome.size - sum(sfs.total)

		genome.size.neut <- 0.25*genome.size
		genome.size.seln <- 0.75*genome.size

		sfs.neut["0"] <- genome.size.neut - sum(sfs.neut)
		sfs.seln["0"] <- genome.size.seln - sum(sfs.seln)
	}
	if(use.manual.sample == TRUE){
		## would have to use this section of code when manually subsampling a full sample output, because don't have the allele frequencies in the sample...
		
		# all possible mut IDs of any type (need for later calcs):
		
		neutral.mut.IDs <- c(neut.muts$mut.ID)
		selected.mut.IDs <- c(seln.muts$mut.ID)
		
		
		all.mutations <- unlist(lapply(as.character(genome.dat[,2]), FUN=strsplit, split=" "))
		just.neut.muts <- all.mutations[which(all.mutations %in% neutral.mut.IDs)]
		just.seln.muts <- all.mutations[which(all.mutations %in% selected.mut.IDs)]
		
		# get the allele frequencies:
		freqs.total <- table(all.mutations)
		freqs.neut <- table(just.neut.muts)
		freqs.seln <- table(just.seln.muts)

		# get the frequency spectra
		sfs.total <- table(freqs.total)
		sfs.neut <- table(freqs.neut)
		sfs.seln <- table(freqs.seln)


		# add on fixed things to the fixed section of the table
		if(is.na(sfs.total[as.character(sample.size)])){
			sfs.total[as.character(sample.size)] <- num.neut.muts.fixed + num.seln.muts.fixed
		}else{sfs.total[as.character(sample.size)] <- sfs.total[as.character(sample.size)] + num.neut.muts.fixed + num.seln.muts.fixed}
		if(is.na(sfs.neut[as.character(sample.size)])){
			sfs.neut[as.character(sample.size)] <- num.neut.muts.fixed
		}else{sfs.neut[as.character(sample.size)] <- sfs.neut[as.character(sample.size)] + num.neut.muts.fixed}
		if(is.na(sfs.seln[as.character(sample.size)])){
			sfs.seln[as.character(sample.size)] <- num.seln.muts.fixed
		}else{sfs.seln[as.character(sample.size)] <- sfs.seln[as.character(sample.size)] + num.seln.muts.fixed}
		
		sfs.total["0"] <- genome.size - sum(sfs.total)
		
		genome.size.neut <- 0.25*genome.size
		genome.size.seln <- 0.75*genome.size

		sfs.neut["0"] <- genome.size.neut - sum(sfs.neut)
		sfs.seln["0"] <- genome.size.seln - sum(sfs.seln)
	}
	
	# each SFS (for DFE) must have the counts for 0 and for fixation
	# if 100 inds sampled, that makes lenght of sfs 201 (0 - 200)
	full.list <- as.character(0:sample.size)
	
	# have to fill in any missing ones with zero to make the neutral match the selected SFS
	neut.missing <- setdiff(full.list, names(sfs.neut))
	sfs.neut[neut.missing] <- 0
	sfs.neut.contents <- as.numeric(paste(sfs.neut))
	sfs.neut.labels <- as.numeric(names(sfs.neut))
	temp.neut <- data.frame(cbind(sfs.neut.labels, sfs.neut.contents))
	ordered.neut <- temp.neut[order(temp.neut[,1]), c(1,2)]
	final.sfs.neut <- ordered.neut[,2]

	seln.missing <- setdiff(full.list, names(sfs.seln))
	sfs.seln[seln.missing] <- 0
	sfs.seln.contents <- as.numeric(paste(sfs.seln))
	sfs.seln.labels <- as.numeric(names(sfs.seln))
	temp.seln <- data.frame(cbind(sfs.seln.labels, sfs.seln.contents))
	ordered.seln <- temp.seln[order(temp.seln[,1]), c(1,2)]
	final.sfs.seln <- ordered.seln[,2]
	
	if(fold == TRUE){
		# fold the site frequency table back on itself
		# take freqs 0-99 and bin with freqs 200-101, then 100 stays on its own at the end (but it's actually 101 because R starts counting at 1, not 0)
		final.sfs.seln <- c(final.sfs.seln[1:num.inds.sampled] + final.sfs.seln[(sample.size+1):(num.inds.sampled+2)], final.sfs.seln[(num.inds.sampled + 1)], rep(0, num.inds.sampled))
		final.sfs.neut <- c(final.sfs.neut[1:num.inds.sampled] + final.sfs.neut[(sample.size+1):(num.inds.sampled+2)], final.sfs.neut[(num.inds.sampled + 1)], rep(0, num.inds.sampled))
	}

	dfe.input <- paste(c(
	"1
", sample.size,"
", paste(c(final.sfs.seln), collapse=" "),"
", paste(c(final.sfs.neut), collapse=" ")
	), collapse="")
	
	if(fold == FALSE){
		write(dfe.input, file=paste(c("Unfolded", filename), collapse=""))
	}
	if(fold == TRUE){
		write(dfe.input, file=paste(c("Folded", filename), collapse=""))	
	}
}
#____________________________________________________________________________________________________#







inds.sampled <- 100
#pop.size <- 10000


args <- commandArgs(trailingOnly=TRUE)

gen <- as.numeric(unlist(strsplit(as.character(args[1]), split="_"))[1])

gsize <- as.numeric(args[4])
gsize.arms <- 0.72*gsize
gsize.centers <- 0.28*gsize

right.arm1 <- 1439999
left.arm2 <- 2560000
right.arm2 <- 5439999
left.arm3 <- 6560000
right.arm3 <- 9439999
left.arm4 <- 10560000
right.arm4 <- 13439999
left.arm5 <- 14560000
right.arm5 <- 17439999
left.arm6 <- 18560000
right.arm6 <- 21439999
left.arm7 <- 22560000
right.arm7 <- 23999999

arms <- c(0:right.arm1, left.arm2:right.arm2, left.arm3:right.arm3, left.arm4:right.arm4, left.arm5:right.arm5, left.arm6:right.arm6, left.arm7:right.arm7)
centers <- setdiff(0:right.arm7, arms)

pop.size <- as.numeric(args[5])

setwd(as.character(args[3]))

#____________________________________________________________________________________________________#


## full data output
full.file <- paste(c("ModifiedSampleOutput_", as.character(args[1])), collapse="")

full.samp.muts.start <- as.numeric(unlist(strsplit(system(paste(c("grep -n Mutations ", full.file), collapse=""), intern=TRUE), split=":"))[1])
full.samp.genomes.start <- as.numeric(unlist(strsplit(system(paste(c("grep -n Genomes ", full.file), collapse=""), intern=TRUE), split=":"))[1])
full.samp.file.end <- as.numeric(head(tail(unlist(strsplit(system(paste(c("wc -l ", full.file), collapse=""), intern=TRUE), split=" ")), n=2), n=1))	

pdat <- read.table(full.file, skip=full.samp.muts.start, nrow=((full.samp.genomes.start-1) - full.samp.muts.start), sep=" ")
names(pdat) <- c("mut.ID", "unique.mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "generation_arose", "mut.prev")
pdat.arms <- pdat[which(pdat$base_position %in% arms) ,]
pdat.centers <- pdat[which(pdat$base_position %in% centers) ,]

pmut.arm.IDs <- pdat.arms$mut.ID
pmut.center.IDs <- pdat.centers$mut.ID
		
gdat <- read.table(full.file, skip=full.samp.genomes.start, nrow=(full.samp.file.end - full.samp.genomes.start), sep="A")
gdat.arms <- cbind(gdat[,1],  rep(NA, 2*inds.sampled))
gdat.centers <- cbind(gdat[,1],  rep(NA, 2*inds.sampled))
for(k in 1:(inds.sampled*2)){
    temp.genome <- as.numeric(unlist(strsplit(as.character(gdat[k,2]), split=" ")))
    temp.arms <- paste(temp.genome[which(temp.genome %in% pmut.arm.IDs)], collapse=" ")
    temp.centers <- paste(temp.genome[which(temp.genome %in% pmut.center.IDs)], collapse=" ")

    gdat.arms[k,2] <- temp.arms
    gdat.centers[k,2] <- temp.centers
}
gdat.arms <- data.frame(gdat.arms)
gdat.centers <- data.frame(gdat.centers)

		

## fixed data output
fixed.mut.id.start <- 2
fdat <- read.table(paste(c("FixedOutput_", paste(unlist(strsplit(as.character(args[1]), split="_"))[-1], collapse="_")), collapse=""), skip=fixed.mut.id.start)
names(fdat) <- c("mut.ID", "unique.mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "gen_arose", "gen.fixed")
fdat.arms <- fdat[which(fdat$base_position %in% arms) ,]
fdat.centers <- fdat[which(fdat$base_position %in% centers) ,]

#____________________________________________________________________________________________________#

if(args[2] == "subsample"){
    outfile.arms <- paste(c("DFE_SFS_subsamp_ARMS_gen", as.character(args[1])), collapse="")
    make.est_dfe.input(poly.dat=pdat.arms, genome.dat=gdat.arms, fixed.dat=fdat.arms,
    generation=gen, num.inds.sampled=inds.sampled, genome.size=gsize.arms,
    filename=outfile.arms, fold=TRUE, use.manual.sample=TRUE)
    outfile.centers <- paste(c("DFE_SFS_subsamp_CENTERS_gen", as.character(args[1])), collapse="")
    make.est_dfe.input(poly.dat=pdat.centers, genome.dat=gdat.centers, fixed.dat=fdat.centers,
    generation=gen, num.inds.sampled=inds.sampled, genome.size=gsize.centers,
    filename=outfile.centers, fold=TRUE, use.manual.sample=TRUE)
}else{
    outfile.arms <- paste(c("DFE_SFS_full_ARMS_gen", as.character(args[1])), collapse="")
    make.est_dfe.input(poly.dat=pdat.arms, genome.dat=gdat.arms, fixed.dat=fdat.arms,
    generation=gen, num.inds.sampled=pop.size, genome.size=gsize.arms,
    filename=outfile.arms, fold=TRUE, use.manual.sample=FALSE)
    outfile.centers <- paste(c("DFE_SFS_full_CENTERS_gen", as.character(args[1])), collapse="")
    make.est_dfe.input(poly.dat=pdat.centers, genome.dat=gdat.centers, fixed.dat=fdat.centers,
    generation=gen, num.inds.sampled=pop.size, genome.size=gsize.centers,
    filename=outfile.centers, fold=TRUE, use.manual.sample=FALSE)
}

	
