
fao <- format_afterOrtho("../app_tagrun/groupsASAHP1.txt", "groups")

data <- read.table("../app_tagrun/TAG_update2.txt", header = T)


tag_matrix <- analyze_OrthoMCL(fao, data, model = "lmeR2ind", species_name = "Treatment", resp = "tgl_fly", rndm1 = "Treatment", rndm2 = "Experiment", sig_digits = 3)
manhat_grp(fao_groups, tag_matrix,tree = file)
manhat_grp(fao, tag_matrix)
pdgplot(data, tag_matrix, COG = "ASAHP01261", species_colname = "Treatment", data_colname = "tgl_fly", xlab = "Taxa", ylab = "TAG Content", tree = file)

file <- system.file('extdata', 'muscle_tree2.dnd', package='MAGNAMWAR')
plot(read.tree(file))

write_mcl(tag_matrix_3, "tagmatrix_short.csv")

phydataerror(phy = file, data = data, mcl_matrix = tag_matrix, species_colname = "Treatment", data_colname = "tgl_fly", xlabel = "TAG Content", COG = "ASAHP03027")

finalkey[] <- lapply(finalkey, as.character)
finalkey[finalkey=="efOG"]<-"efog"

write_mcl(tag_matrix, "tag_4_17.csv")
