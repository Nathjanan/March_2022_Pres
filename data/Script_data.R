#################################### make individual smiles
df = read.csv('vEGFR.smi',sep = "" , header = FALSE, na.strings ="", stringsAsFactors= F)
for (i in 1:nrow(df)) {
  each_smi = df[i, ]
  each_smi_sel = each_smi[, c(1, 2)]
  name_id = as.character(each_smi$V2)
  write.table(each_smi_sel, sep = "\t", file = paste0('ligands/smi/',
                                                      name_id, ".smi"),
              col.names = FALSE, row.names = FALSE, quote = FALSE)
  
}

############################## make sdf files 

library(parallel)
library(doSNOW)
library(benchmarkme)
smi_files = list.files(path= "ligands/smi/",
                       pattern = ".smi", full.names = FALSE, recursive = FALSE)

cl <- makeCluster(benchmarkme::get_cpu()$no_of_cores)
registerDoSNOW(cl)

foreach(i = 1:length(smi_files)) %dopar% {
  smi_file = smi_files[i]
  input_smi = paste0("ligands/smi/",
                     smi_file)
  output_sdf = paste0("ligands/sdf/", 
                      gsub('.{3}$', "", smi_file), "sdf")
  command = paste0("obabel ", input_smi, " -O ", output_sdf, " --gen3D -p")
  system(command, intern = FALSE)
  
} 
stopCluster(cl)



##### docked sdf molecules
library(parallel)
library(doSNOW)
library(benchmarkme)
sdf_files = list.files(path= "ligands/sdf/",
                       pattern = ".sdf", full.names = FALSE, recursive = FALSE)

cl <- makeCluster(benchmarkme::get_cpu()$no_of_cores)
registerDoSNOW(cl)

foreach(i = 1:length(sdf_files)) %dopar% {
  sdf_file = sdf_files[i]
  input_sdf = paste0("ligands/sdf/",
                     sdf_file)
  output_sdf = paste0("ligands/dock/", 
                      gsub('.{3}$', "", sdf_file), "sdf")
  command = paste0("./smina.static ", 
                   ' -r prepared_4asd_protein.pdbqt --autobox_ligand ligand_4ASD_BAX.pdb', ' -l ',
                   input_sdf, " -o ", output_sdf,  
                   ' --num_modes 1 --seed 0 --cpu 1 --size_x 25 --size_y 25 --size_z 25 --exhaustiveness 1')
  system(command)
  
} 

stopCluster(cl)

#### dock the over head molecles 
empty_files = system("find . -name '*.sdf' -size 0", intern = TRUE)
want_files = c()
for (i in 1:length(empty_files)) {
  sdf_files = strsplit(empty_files, "/")
  sdf_file = sdf_files[[i]][4]
  want_files = c(want_files, sdf_file)
}



for (i in 1:length(want_files)) {
  sdf_file = want_files[i]
  input_sdf = paste0("ligands/sdf/",
                     sdf_file)
  output_sdf = paste0("ligands/dock/", 
                      gsub('.{3}$', "", sdf_file), "sdf")
  command = paste0("./smina.static ", 
                   ' -r prepared_4asd_protein.pdbqt --autobox_ligand ligand_4ASD_BAX.pdb', ' -l ',
                   input_sdf, " -o ", output_sdf,  
                   ' --num_modes 1 --seed 0 --cpu 1 --size_x 25 --size_y 25 --size_z 25 --exhaustiveness 1')
  system(command)
  
} 
