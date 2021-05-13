# the function:

make.est_alpha_omega.input <- function(poly.dat, genome.dat, fixed.dat, generation, num.inds.sampled, genome.size, filename, use.manual.sample=FALSE){

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
	
	# take just the number of fixations (last number in the sfs) as the "selected" and "neutral" "differences" because it's supposed to be a species level thing
	seln.fixed.sites <- final.sfs.seln[(sample.size + 1)]
	neut.fixed.sites <- final.sfs.neut[(sample.size + 1)]


	seln.line <- paste(c("1", as.character(format(0.75*genome.size, scientific=FALSE)), seln.fixed.sites), collapse=" ")
	neut.line <- paste(c("0", as.character(format(0.25*genome.size, scientific=FALSE)), neut.fixed.sites), collapse=" ")
	
	alpha_omega.input <- paste(c(seln.line, neut.line), collapse="\n")

	write(alpha_omega.input, file=paste(c("Divergence", filename), collapse=""))	
}
#____________________________________________________________________________________________________#







inds.sampled <- 100
#pop.size <- 10000


args <- commandArgs(trailingOnly=TRUE)

gen <- as.numeric(unlist(strsplit(as.character(args[1]), split="_"))[1])

gsize <- as.numeric(args[4])

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
		
gdat <- read.table(full.file, skip=full.samp.genomes.start, nrow=(full.samp.file.end - full.samp.genomes.start), sep="A")



## fixed data output
fixed.mut.id.start <- 2
fdat <- read.table(paste(c("FixedOutput_", paste(unlist(strsplit(as.character(args[1]), split="_"))[-1], collapse="_")), collapse=""), skip=fixed.mut.id.start)
names(fdat) <- c("mut.ID", "unique.mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "gen_arose", "gen.fixed")

#____________________________________________________________________________________________________#

if(args[2] == "subsample"){
	outfile <- paste(c("_alphaOmega_subsamp_gen", as.character(args[1])), collapse="")
	make.est_alpha_omega.input(poly.dat=pdat, genome.dat=gdat, fixed.dat=fdat, 
	generation=gen, num.inds.sampled=inds.sampled, genome.size=gsize, 
	filename=outfile, use.manual.sample=TRUE)
}else{
	outfile <- paste(c("_alphaOmega_full_gen", as.character(args[1])), collapse="")
	make.est_alpha_omega.input(poly.dat=pdat, genome.dat=gdat, fixed.dat=fdat, 
	generation=gen, num.inds.sampled=pop.size, genome.size=gsize, 
	filename=outfile, use.manual.sample=FALSE)
}

	
