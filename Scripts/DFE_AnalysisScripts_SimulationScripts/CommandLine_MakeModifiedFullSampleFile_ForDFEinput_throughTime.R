
## MAKE A FULL DATA FILE, MODIFIED FROM SAMPLE OUTPUT
##	AT ALL GENERATIONS SAMPLED


args <- commandArgs(trailingOnly=TRUE)

# args 1 is basename of file
# args 2 is working directory (inputs dir)

setwd(as.character(args[2]))


taking.in.file <- paste(c("SampleOutput_", as.character(args[1])), collapse="")

gens <- matrix(unlist(strsplit(system(paste(c("grep -n '#OUT: ' ", taking.in.file), collapse=""), intern=TRUE), split=" ")), byrow=TRUE, ncol=4)[,2]

for(i in 1:length(gens)){
	spitting.out.file.name <- paste(c("ModifiedSampleOutput_", gens[i], "_", as.character(args[1])), collapse="")
	
	# subset the file to include the mutations section and the individuals section (eventually become polydat and genodat)
	
	start.line <- as.numeric(unlist(strsplit(system(paste(c("grep -n '#OUT: ", gens[i], " ' ", taking.in.file), collapse=""), intern=TRUE), split=":"))[1])

	if(i == length(gens)){	# then it's the last gen
		system(paste(c("tail -n+", start.line, " ", taking.in.file, " > ", spitting.out.file.name), collapse=""))
	}else{					# otherwise we're still in the middle of the file
		end.line <- as.numeric(unlist(strsplit(system(paste(c("grep -n '#OUT: ", gens[i+1], " ' ", taking.in.file), collapse=""), intern=TRUE), split=":"))[1])
		system(paste(c("sed -n '", start.line, ",", (end.line-1), " p' ", taking.in.file, " > ", spitting.out.file.name), collapse=""))
	}
}

