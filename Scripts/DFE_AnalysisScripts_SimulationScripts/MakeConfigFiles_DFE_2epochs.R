
args <- commandArgs(trailingOnly=TRUE)

# args[1]: T/F -- False = delet only, True = with beneficials, so make an alpha_omega config file
# args[2]: string -- the directory containng the sfs file to be analyzed
# args[3]: numeric -- the sfs file to be analyzed
# args[4]: numeric -- mean s value	from DFE example, = -0.1
# args[5]: numeric -- mean beta		from DFE example, = 0.5



# make class 0 config file:
line1 <- "data_path_1 /cap1/kgilbert/New_DFE_alpha/ProgramData
"
sfs.input <- paste(c("sfs_input_file ", args[2], args[3],"
"), collapse="")
class0.output <- paste(c("est_dfe_results_dir Outputs/OutputClass0_", unlist(strsplit(as.character(args[3]), split=".txt"))[1],"
"), collapse="")
line4.end <- "site_class 0
fold 1
epochs 2
search_n2 1
t2_variable 1
t2 1000"

file.text <- paste(c(line1, sfs.input, class0.output, line4.end), collapse="")
write(file.text, file=paste(c(args[2],"config_class0.txt"), collapse=""))



# make class 1 config file:
class1.output <- paste(c("est_dfe_results_dir Outputs/OutputClass1_", unlist(strsplit(as.character(args[3]), split=".txt"))[1], "
"), collapse="")
class0output.class1input <- paste(c("est_dfe_demography_results_file Outputs/OutputClass0_", unlist(strsplit(as.character(args[3]), split=".txt"))[1], "/est_dfe.out
"), collapse="")
line5.end <- paste(c("site_class 1
fold 1
epochs 2
mean_s_variable 1
mean_s ", args[4], "
beta_variable 1
beta ", args[5]), collapse="")

file.text <- paste(c(line1, sfs.input, class1.output, class0output.class1input, line5.end), collapse="")
write(file.text, file=paste(c(args[2],"config_class1.txt"), collapse=""))



# make an omega_alpha config file
##	if(args[1] == TRUE){
##	always make alpha omega file for now - just see what it estimates in the ones with delet muts only
line1 <- "data_path_1 /cap1/kgilbert/DFE_alpha/ProgramData/
"
	name.of.div.file <- unlist(strsplit(as.character(args[3]), split="SFS"))[2]
	divergence.data.line <-	paste(c("divergence_file ", args[2], "Divergence_alphaOmega", name.of.div.file,"
"), collapse="")
	new.results.line <- paste(c("est_alpha_omega_results_file output_est_alpha_omega", name.of.div.file, ".out
"), collapse="")
	prev.dfe.results <- paste(c("est_dfe_results_file Outputs/OutputClass1_", unlist(strsplit(as.character(args[3]), split=".txt"))[1], "/est_dfe.out
"), collapse="")
	prev.neut.results <- paste(c("neut_egf_file Outputs/OutputClass0_", unlist(strsplit(as.character(args[3]), split=".txt"))[1], "/neut_egf.out
"), collapse="")
	prev.sel.results <- paste(c("sel_egf_file Outputs/OutputClass1_", unlist(strsplit(as.character(args[3]), split=".txt"))[1], "/sel_egf.out
"), collapse="")
	options <- "do_jukes_cantor         0
remove_poly             0"

	file.text <- paste(c(line1, divergence.data.line, new.results.line, prev.dfe.results, prev.neut.results, prev.sel.results, options), collapse="")
	write(file.text, file=paste(c(args[2],"config_alpha_omega.txt"), collapse=""))

##	}
