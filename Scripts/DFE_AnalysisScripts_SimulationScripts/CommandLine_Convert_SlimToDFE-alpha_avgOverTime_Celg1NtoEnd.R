
## MAKE MODIFIED DATA FILES FROM OUTPUT FOR EACH GEN

# ARG 1 = Sample File base name
# Arg 2 = directory of file

args <- commandArgs(trailingOnly=TRUE)






num.inds.sampled <- 100

# args 1 is basename of file
# args 2 is working directory (inputs dir)
# args 3 is genosize
# args 4 is pop size

setwd(as.character(args[2]))

gsize <- as.numeric(args[3])
pop.size <- as.numeric(args[4])

inds.sampled <- 100


taking.in.file <- paste(c("SampleOutput_", as.character(args[1])), collapse="")

gens <- matrix(unlist(strsplit(system(paste(c("grep -n '#OUT: ' ", taking.in.file), collapse=""), intern=TRUE), split=" ")), byrow=TRUE, ncol=4)[,2]
#final.gen <- format(as.numeric(tail(gens, n=1))+pop.size, scientific=FALSE)


# remove 1N, 2N, 3N - 9N:
#gens <- gens[-c(1:9)]

# add last gen
#gens <- c(gens, as.character(final.gen))

for(i in gens){
	spitting.out.file.name <- paste(c("ModifiedSampleOutput_", i, "_", as.character(args[1])), collapse="")

	# subset the file to include the mutations section and the individuals section (eventually become polydat and genodat)

	taking.in.file <- paste(c("SampleOutput_", as.character(args[1])), collapse="")
	start.line <- as.numeric(unlist(strsplit(system(paste(c("grep -n '#OUT: ", i, " ' ", taking.in.file), collapse=""), intern=TRUE), split=":"))[1])
	if(i == tail(gens, n=1)){	# this is the last gen
		end.line <- as.numeric(unlist(strsplit(system(paste(c("wc -l ", taking.in.file), collapse=""), intern=TRUE), split=" "))[1])
                system(paste(c("sed -n '", start.line, ",", end.line, " p' ", taking.in.file, " > ", spitting.out.file.name), collapse=""))
	}else{					# this is gens in the middle
		temp.pop.size <- format(as.numeric(pop.size), scientific=FALSE)
		temp.i <- format(as.numeric(i), scientific=FALSE)
		end.line <- as.numeric(unlist(strsplit(system(paste(c("grep -n '#OUT: ", format(as.numeric(temp.i)+(as.numeric(temp.pop.size)/10), scientific=FALSE), " ' ", taking.in.file), collapse=""), intern=TRUE), split=":"))[1])
		system(paste(c("sed -n '", start.line, ",", (end.line-1), " p' ", taking.in.file, " > ", spitting.out.file.name), collapse=""))
	}
}


# spits out SFS at each gen into a "Modified" sample file





# the function for calculating SFS:

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







	

# get the SFS for each gen from the Modified files, output either Folded or unfolded SFS

for(j in gens){
	gen <- format(j, scientific=FALSE)
	full.file <- paste(c("ModifiedSampleOutput_", gen, "_", as.character(args[1])), collapse="")
	
	full.samp.muts.start <- as.numeric(unlist(strsplit(system(paste(c("grep -n Mutations ", full.file), collapse=""), intern=TRUE), split=":"))[1])
	full.samp.genomes.start <- as.numeric(unlist(strsplit(system(paste(c("grep -n Genomes ", full.file), collapse=""), intern=TRUE), split=":"))[1])
	full.samp.file.end <- as.numeric(head(tail(unlist(strsplit(system(paste(c("wc -l ", full.file), collapse=""), intern=TRUE), split=" ")), n=2), n=1))	
	
	pdat <- read.table(full.file, skip=full.samp.muts.start, nrow=((full.samp.genomes.start-1) - full.samp.muts.start), sep=" ")
	names(pdat) <- c("mut.ID", "unique.mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "generation_arose", "mut.prev")
			
	gdat <- read.table(full.file, skip=full.samp.genomes.start, nrow=(full.samp.file.end - full.samp.genomes.start), sep="A")
	
	## fixed data output
	fixed.mut.id.start <- 2
	fdat <- read.table(paste(c("FixedOutput_", as.character(args[1])), collapse=""), skip=fixed.mut.id.start)
	names(fdat) <- c("mut.ID", "unique.mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "gen_arose", "gen.fixed")
	
	#____________________________________________________________________________________________________#
	
	outfile <- paste(c("SFS_gen", gen, "_", as.character(args[1])), collapse="")
	make.est_dfe.input(poly.dat=pdat, genome.dat=gdat, fixed.dat=fdat, 
	generation=gen, num.inds.sampled=inds.sampled, genome.size=gsize, 
	filename=outfile, fold=TRUE, use.manual.sample=TRUE)

}





#____________________________________________________________________________________________________#

# MAKE THE AVERAGE DFE

sample.size <- num.inds.sampled*2

outfile <- paste(c("Avg_SFS_10NtoEnd_", as.character(args[1])), collapse="")
files <- system(paste(c("ls Folded*", as.character(args[1])), collapse=""), intern=TRUE)
sfs.dat <- vector("list", length(gens))
for(j in 1:length(gens)){
	sfs.dat[[j]] <- as.matrix(read.table(files[j], skip=2))
}
avged.sfs.dat <- apply(simplify2array(sfs.dat), MARGIN=1:2, FUN=mean)
	
dfe.input <- paste(c(
'1
', sample.size,'
', paste(c(avged.sfs.dat[1,]), collapse=" "),'
', paste(c(avged.sfs.dat[2,]), collapse=" ")
), collapse="")

write(dfe.input, file=outfile)




#____________________________________________________________________________________________________#
# make the divergence file from number of fixed sites in the fixed file only, this goes from 0N to 10N

fixed.mut.id.start <- 2

outfile <- paste(c("Divergence_", as.character(args[1])), collapse="")
fixed.file <- read.table(paste("FixedOutput_", as.character(args[1]), sep=""), skip=fixed.mut.id.start)
names(fixed.file) <- c("mut.ID", "unique.mut.ID", "mut.type", "base_position", "seln_coeff", "dom_coeff", "subpop_ID", "gen_arose", "gen.fixed")

fixed.neut.muts <- c(which(fixed.file$mut.type == "m1"))
fixed.seln.mut.IDs <- fixed.file$mut.ID[-fixed.neut.muts]
fixed.neut.mut.IDs <- fixed.file$mut.ID[fixed.neut.muts]

num.neut.muts.fixed <- length(fixed.neut.mut.IDs)
num.seln.muts.fixed <- length(fixed.seln.mut.IDs)

genome.size.neut <- 0.25*gsize
genome.size.seln <- 0.75*gsize

seln.line <- paste(c("1", as.character(format(genome.size.seln, scientific=FALSE)), num.seln.muts.fixed), collapse=" ")
neut.line <- paste(c("0", as.character(format(genome.size.neut, scientific=FALSE)), num.neut.muts.fixed), collapse=" ")

alpha_omega.input <- paste(c(seln.line, neut.line), collapse="\n")

write(alpha_omega.input, file=outfile)	



