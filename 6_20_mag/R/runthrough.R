# 
# # RASTtogbk(inpuut_fasta, input_reference, out_path)
# 
# 
# # formatted_file <- format_MCLfastas("../Documents/COLLEGE/ChastonLab/RUN/FASTAS/")
# 
# # RUN OrthoMCL SOFTWARE #
# 
# stl_afterortho<- format_afterOrtho("../stl2/groupsstv_tag_lifs.txt", format="groups")
# 
# matrix_TAG_old <- analyze_OrthoMCL(after_ortho_format1, pheno_data, "Treatment", "lm", resp="RespVar")
# 
# pheno_data2 <- read.table("../TAG_data.txt", header =T)
# starv_data <- read.csv("../stl2/starv_data.txt", header =T)
# 
# 
# pheno_data2[1,]
# 
# TAG_mixmtrx <- analyze_OrthoMCL(stl_afterortho, pheno_data2, "Treatment", "lmeR2ind", resp="tgl_fly", rndm1="Treatment", rndm2="Experiment")
# 
# 
# 
# install.packages("survival")
# install.packages("parallel")
# install.packages("foreach")
# install.packages("doParallel")
# 
# starv_mtrx <- analyze_OrthoMCL(stl_afterortho, starv_data, "abbrev", "survmulti", time="time", event="event", rndm1="exp", rndm2="vial", multi=19)
# 
# 
# 
# 
# 
# ### john running lifespan
# 
# life_data <- read.csv('/home/sblanch/Downloads/allF_RLcensored.csv', header=T)
# table(list(life_data$trt2))
# life_mtrx <- analyze_OrthoMCL(haplo_data = stl_afterortho, var_file = life_data,species_name =  "trt2", model = "survmulticensor",fix2 = "baclost",rndm1 = "exp",rndm2 = "vial2",multi = 19,time = "t1",event = "event",time2 = "t2")
# 
# 
# 
# 
# 
# ### lindy and corinne running PERSISTENCE
# 
# pers_data <- read.csv("/home/sblanch/Downloads/allF_RLcensored_persist_FINAL.csv", header =T)
# 
# pers2 <- aggregate(Days.post.Picking ~ trt2 + rep +  exp + rep + vial2, pers_data, FUN="mean")
# 
# pers2[1,]
# 
# 
# #mcl_mtrx1 <- analyze_OrthoMCL(stl_afterortho, pers_data, "trt2", "lmeR2nest", resp="Days.post.Picking", rndm1="trt2", rndm2="exp")
# 
# persistence_mtrx <- analyze_OrthoMCL(stl_afterortho, pers_data, "trt2", "lmeR2ind", resp="Days.post.Picking", rndm1="trt2", rndm2="exp")
# 
# library(seqinr)
# 
# pick_file <- pick_repseq(stl_afterortho, "compliantFasta-original/")
# joined_pers_mtrx <- join_repset(pick_file, persistence_mtrx)
# qqplotter(persistence_mtrx)
# 
# write_mcl(joined_pers_mtrx, "")
# 
# ### Alec's STARVATION
# 
# starv_data <- read.csv("../stl2/starv_data.txt", header =T)
# 
# starv_mtrx <- analyze_OrthoMCL(stl_afterortho, starv_data, "abbrev", "survmulti", time="time", event="event", rndm1="exp", rndm2="vial", multi=18)
# 
# pick_file <- pick_repseq(stl_afterortho, "compliantFasta-original/")
# 
# dim(starv_mtrx)
# 
# dim(persistence_mtrx)
# 
# setwd('~/Desktop/mycode/drchaston_lab/stl2/')
# 
# starv_mtrx6 <- analyze_OrthoMCL(stl_afterortho, starv_data, "abbrev", "survmulti", time="time", event="event", rndm1="exp", rndm2="vial", multi=18)
# 
# starv_mtrx7 <- read.csv('outstarve.2csv', header=T)
# colnames(starv_mtrx7) <- c("COG","p-val","corrected_p-val","mean_COGContain","mean_COGLack","taxa_contain","taxa_miss")
# 
# 
# 
# ## format file by unlisting it - this will be tricky. i might actually need to write a perl script to parse it
# 
# library(seqinr)
# getwd()
# setwd("..")
# pick_file <- pick_repseq(stl_afterortho, "compliantFasta-original/")
# joined_starv_mtrx <- join_repset(pick_file, starv_mtrx7)
# 
# write_mcl(joined_starv_mtrx, "join.csv")
# 
# qqplotter(joined_starv_mtrx)
# 
# 
# 
# 
