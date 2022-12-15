rm(list=ls())

#options(max.print=1000000)
library(data.table)
library(readr)
library(haven)
library(plyr)
library(dplyr)
library(vegan)
library(scales)
library(RColorBrewer)
library(grid)
library(pheatmap)
library(lme4)
library(ggplot2)
library(cowplot)
library(ggsignif)
library(car)
library(Maaslin2) 
library(gridExtra)
library(knitr)
library(readxl)
library(ggpubr)
library(nlme)

### Pre-processing steps 
#### (1) read in channing metadata diet data 

z_mlvs_exposure <- read.csv('/n/home02/rmehta/balskus/mlvs_exposure_wenjie.csv', header=T, row.names = 1)

# some of your nutrients have w1avg and w2avg in the name
colnames(z_mlvs_exposure) <- gsub("w1avg","w1", colnames(z_mlvs_exposure))
colnames(z_mlvs_exposure) <- gsub("w2avg","w2", colnames(z_mlvs_exposure))

colnames(z_mlvs_exposure) <- gsub("plasma1","w1", colnames(z_mlvs_exposure))
colnames(z_mlvs_exposure) <- gsub("plasma2","w2", colnames(z_mlvs_exposure))

z_mlvs_agebmi <- z_mlvs_exposure [ , colnames(z_mlvs_exposure) %in% c('agemlvs', 'bmi12','act10')] 

keep_metadata_cum <- c('aofib10v','calor10v','frtaf10v','ceraf10v','vegaf10v')
z_metadata_cum <- z_mlvs_exposure[, keep_metadata_cum]

keep_metadata_avg <- c(
  #ffq questionnaire
  'acid_avg', 'alc_avg', 'abx_avg', 'prep_avg',
  'selfdiet_avg', 'probx_avg', 'bristolcat_avg', 'bristol_avg', 'oral_avg', 'yogurt_avg',
  #ffq energy/fiber
  'calor_avg', 'aofib_avg', 
  #ffq patterns 
  'west_avg', 'westq_avg', 'prud_avg', 'prudq_avg',
  #ffq nutrients
  'sulfur_avg', 'prot_sulf_avg', 'nonprot_sulf_avg', 'aprot_avg', 'dprot_avg',
  'prot_avg', 'vprot_avg', 'afat_avg', 'dfat_avg', 'tfat_avg', 'vfat_avg',
  'carbo_avg', 'sucr_avg', 'fruct_avg', 'lact_avg', 'st_avg', 'glu_avg',
  'calc_wo_avg', 'calc_avg', 'iron_avg', 'iron_wo_avg', 'heme_avg', 'fdfol_avg',
  'fol98_avg', 'fol98_wo_avg', 'folic_avg', 'dfe_avg', 'b12_avg', 'b12_wo_avg',
  'choline_avg', 'choline_wo_avg', 'dvitk_avg', 'vitk_avg', 'vitk2_m4_avg',
  'vitk_wo_avg', 'satfat_avg', 'sft07_avg', 'mft07_avg', 'monfat_avg',
  'ply07_avg', 'poly_avg', 'trn07_avg', 'ttrn07_avg', 'alco_avg', 'vitd_avg',
  'vitd_wo_avg', 'dvitd_avg', 'glcsin_avg', 
  #ffq food groups
  'promeatavg', 'redmeatavg',
   'fishavg', 'poultavg', 'eggsavg', 'butteravg', 'margavg',
  'lowdaiavg', 'highdaiavg', 'wineavg', 'liquoravg', 'beeravg', 'teaavg',
  'coffeeavg', 'fruitavg', 'frujuavg', 'cruvegavg', 'yelvegavg', 'tomatoavg',
  'leafvegavg', 'legumeavg', 'othvegavg', 'potatoavg', 'frenchavg', 'wholegavg',
  'refingavg', 'pizzaavg', 'sugdrkavg', 'lowdrkavg', 'snackavg', 'nutsavg',
  'mayoavg', 'crmsoupavg', 'sweetsavg')

z_metadata_avg <- z_mlvs_exposure[, keep_metadata_avg]

keep_metadata_w1 <- c(
  #questionnaire
  'acid_w1', 'alc_w1', 'abx_w1', 'prep_w1',
  'selfdiet_w1', 'probx_w1', 'bristolcat_w1', 'bristolcat5', 'bristolcat6', 'oral_w1', 'yogurt_w1',
  'bristol_w1', 'bristol5', 'bristol6',
  #ddr1
  'calor_fs_dr_w1', 'a_aofib_fs_dr_w1', 'a_aofib_fo_dr_w1', 'a_pect_fo_dr_w1', 'dchocffq1','mchocffq1',
  #ddr1 nutrients
  'a_pro_fs_dr_w1', 'a_vpro_fo_dr_w1', 'a_aprot_fo_dr_w1',
  'a_fat_fs_dr_w1', 'a_carbo_fs_dr_w1', 'a_sucr_fs_dr_w1',
  'a_fruct_fs_dr_w1', 'a_lact_fo_dr_w1', 'a_star_fo_dr_w1',
  'a_glu_fs_dr_w1', 'a_calc_fs_dr_w1', 'a_iron_fs_dr_w1',
  'a_dfe_fs_dr_w1', 'a_fdfol_fo_dr_w1', 'a_fol98_fo_dr_w1',
  'a_folic_fs_dr_w1', 'a_b12_fs_dr_w1', 'a_choline_fs_dr_w1',
  'a_vk_fo_dr_w1', 'a_sfa_fs_dr_w1', 'a_mfa_fs_dr_w1',
  'a_pfa_fs_dr_w1', 'a_ttfa_fo_dr_w1', 'a_alco_fo_dr_w1',
  'a_vd_fs_dr_w1',
  'crp_w1','tchdl_w1','logcrp_w1','logtchdl_w1', 'logtc_w1' , 'loghdl_w1', 'logtg_w1'
)

keep_metadata_w2 <- c(
  #questionnaire
  'acid_w2', 'alc_w2', 'abx_w2', 'prep_w2',
  'selfdiet_w2', 'probx_w2', 'bristolcat_w2', 'bristolcat7', 'bristolcat8', 'oral_w2', 'yogurt_w2',
  'bristol_w2', 'bristol7', 'bristol8',
  #ddr2 energy/fiber
  'calor_fs_dr_w2', 'a_aofib_fs_dr_w2', 'a_aofib_fo_dr_w2', 'a_pect_fo_dr_w2', 'dchocffq2','mchocffq2',
  #ddr2 nutrients
  'a_pro_fs_dr_w2', 'a_vpro_fo_dr_w2', 'a_aprot_fo_dr_w2',
  'a_fat_fs_dr_w2', 'a_carbo_fs_dr_w2', 'a_sucr_fs_dr_w2',
  'a_fruct_fs_dr_w2', 'a_lact_fo_dr_w2', 'a_star_fo_dr_w2',
  'a_glu_fs_dr_w2', 'a_calc_fs_dr_w2', 'a_iron_fs_dr_w2',
  'a_dfe_fs_dr_w2', 'a_fdfol_fo_dr_w2', 'a_fol98_fo_dr_w2',
  'a_folic_fs_dr_w2', 'a_b12_fs_dr_w2', 'a_choline_fs_dr_w2',
  'a_vk_fo_dr_w2', 'a_sfa_fs_dr_w2', 'a_mfa_fs_dr_w2',
  'a_pfa_fs_dr_w2', 'a_ttfa_fo_dr_w2', 'a_alco_fo_dr_w2',
  'a_vd_fs_dr_w2',
  'crp_w2','tchdl_w2','logcrp_w2','logtchdl_w2', 'logtc_w2' , 'loghdl_w2', 'logtg_w2'
)

z_metadata_w1 <- z_mlvs_exposure[, keep_metadata_w1]
z_metadata_w2 <- z_mlvs_exposure[, keep_metadata_w2]
z_metadata_all <- cbind(z_mlvs_agebmi, z_metadata_cum, z_metadata_avg, z_metadata_w1, z_metadata_w2)


z_metadata_all$choc_avg <- (z_metadata_all$dchocffq1+z_metadata_all$mchocffq1+z_metadata_all$dchocffq2+z_metadata_all$mchocffq2)/2
z_metadata_all$choc_avg[is.na(z_metadata_all$choc_avg)] <- 0

# metadata_all doesn't have a participants variable 
z_metadata_all$id <- rownames(z_metadata_all)
z_metadata_all$id <- as.factor(z_metadata_all$id) 

####### read in Min's MGX files and average these and then asinsqrt ###############

gphcdh<- read.table("/n/home02/rmehta/balskus/mindata/Dehydroxylase_MGX_aug.txt",header=T,sep="\t")
gphcdh$ID1 <- substr(gphcdh$Sample.ID,1,8)
gphcdh$ID1 <- as.factor(gphcdh$ID1)

###### ID key ####

idkey <- read.table("/n/home02/rmehta/balskus/MLVSidKey.txt",header=T,sep="\t") 
idkey$ID1 <- as.factor(idkey$ID1) 
matched1 <- inner_join(idkey,gphcdh,by="ID1")
matched1$id <- as.factor(matched1$id) 
gphmeta <- inner_join(matched1,z_metadata_all,by="id")

dim(gphmeta)
#[1] 909 219

############################################################################


#### MAASLIN2 time ##############

##### start with patterns, exploratory, adjust for age, bmi, and the other dietary pattern which hones in on gph meta2 ###

library(tidyr)
gphmeta2 <- gphmeta %>% filter(west_avg > -2) %>% drop_na(eCadh_mgx)
df_input_data <- gphmeta2 %>% select(c(6,7,8))
df_input_data <- asin(sqrt(df_input_data))
df_input_metadata_noID <- gphmeta2 %>% select(2,29,31,9,11,10)

dim(df_input_data) #899
dim(df_input_metadata_noID) #899 

str(df_input_data)
str(df_input_metadata_noID)

fit_data_pattern = Maaslin2(
    input_data = df_input_data, 
    input_metadata = df_input_metadata_noID, 
    output = "maaslin_output_pattern_gphSeptMgx", 
    transform="none",
    heatmap_first_n = 100,
    max_significance = 0.25,
    normalization="none",
    standardize=TRUE,
    fixed_effects = c("prud_avg" ,"west_avg","agemlvs" , "bmi12" ),random_effects=c("ID1"))

#west_EAgenes<-ggplot(gphmeta2 ,aes(x=west_avg,y=asin(sqrt(Gedh2)))) + geom_point() + geom_smooth(method = "lm") + theme_cowplot() 

library(ggpmisc)
my.formula <- y~x 
prud_gph<-ggplot(gphmeta2 , aes(x = prud_avg, y = asin(sqrt(GpHcdh_mgx)))) +
    geom_point(col = "#3C5488FF",alpha=0.8) + theme_bw()+ geom_smooth(method = "lm",
        col = "#C42126", formula = my.formula,
        #se = FALSE,
        size = 1) +theme_cowplot(20) +
theme(strip.background = element_blank(),axis.text.x = element_text(angle=45,vjust = 1, hjust = 1))+xlab("Prudent pattern score")+ ylab("Rel. abund. (asin sqrt)")+
scale_x_continuous(labels = function(x) format(x, scientific = TRUE))# + annotate("text", x = 3, y = 0.13, size = 5, label = "coef: 0.0034\nFDR: 2e-4\nn=909", fontface = "italic",hjust = 1)
ggsave(filename='./gphprud2_aug22.png', plot=prud_gph, width = 6, height = 6, dpi = 600) 

### food groups -- null for components (of western pattern) 

df_input_data <- gphmeta2 %>% select(c(6,7))
df_input_data <- asin(sqrt(df_input_data))
df_input_metadata_noID <- gphmeta2 %>% select(c(2,9,11,29,31,81:114,79))

dim(df_input_data) #909
dim(df_input_metadata_noID) #909

str(df_input_data)
str(df_input_metadata_noID)
 
#this is foods 
fit_data_food = Maaslin2(
    input_data = df_input_data, 
    input_metadata = df_input_metadata_noID, 
    output = "maaslin_output_foods_aug22", 
    transform="none",
    heatmap_first_n = 100,
    max_significance = 0.10,
    normalization="none",
    standardize=TRUE,
    fixed_effects = c("promeatavg" ,"redmeatavg", "fishavg" ,   "poultavg",
"eggsavg"    ,"butteravg" , "margavg"    ,"lowdaiavg" , "highdaiavg",
 "wineavg"   , "liquoravg" , "beeravg"   , "teaavg"   ,  "coffeeavg",
"fruitavg"   ,"frujuavg"   ,"cruvegavg"  ,"yelvegavg"  ,"tomatoavg",
 "leafvegavg", "legumeavg"  ,"othvegavg" , "potatoavg" , "frenchavg",
"wholegavg"  ,"refingavg" , "pizzaavg"  , "sugdrkavg" , "lowdrkavg",
 "snackavg"  , "nutsavg"  , "choc_avg","prud_avg" ,"west_avg","agemlvs","bmi12"),random_effects=c("ID1"))

library(nlme)
b <- lme(asin(sqrt(GpHcdh_mgx))~coffeeavg + cruvegavg + agemlvs + bmi12 + abx_avg+bristol_avg, random=~1|ID1, data=gphmeta ,na.action=na.omit,control=lmeControl(returnObject=TRUE))
summary(b)
#coffeeavg    0.00162268 0.000581818 295  2.788978  0.0056
#cruvegavg    0.00670172 0.002132856 295  3.142136  0.0018
#agemlvs     -0.00013649 0.000197699 295 -0.690418  0.4905
#bmi12       -0.00035132 0.000229199 295 -1.532838  0.1264
#abx_avg     -0.00145135 0.001796646 295 -0.807812  0.4198
#bristol_avg -0.00281976 0.000911754 295 -3.092677  0.0022

cruv_gph <-ggplot(gphmeta2 , aes(x = cruvegavg, y = asin(sqrt(GpHcdh_mgx)))) +
    geom_point(col = "#3C5488FF",alpha=0.8) + theme_bw()+ geom_smooth(method = "lm",
        col = "#C42126", formula = my.formula,
        #se = FALSE,
        size = 1) +theme_cowplot(20) +
theme(strip.background = element_blank(),axis.text.x = element_text(angle=45,vjust = 1, hjust = 1))+xlab("Cruciferous veg.")+ ylab("Rel. abund. (asin sqrt)")+
scale_x_continuous(labels = function(x) format(x, scientific = TRUE)) #+ annotate("text", x = 3, y = 0.13, size = 5, label = "coef: 0.0021\nFDR: 0.07\nn=909", fontface = "italic",hjust = 1)
ggsave(filename='./gphcruv2_Aug22.png', plot=cruv_gph, width = 6, height = 6, dpi = 600)

coff_gph <-ggplot(gphmeta2 , aes(x = coffeeavg, y = asin(sqrt(GpHcdh_mgx)))) +
    geom_point(col = "#3C5488FF",alpha=0.8) + theme_bw()+ geom_smooth(method = "lm",
        col = "#C42126", formula = my.formula,
        #se = FALSE,
        size = 1) +theme_cowplot(20) +
theme(strip.background = element_blank(),axis.text.x = element_text(angle=45,vjust = 1, hjust = 1))+xlab("Coffee")+ ylab("Rel. abund. (asin sqrt)")+
scale_x_continuous(labels = function(x) format(x, scientific = TRUE))#+ annotate("text", x = 6, y = 0.13, size = 5, label = "coef: 0.0015\nFDR: 0.07\nn=909", fontface = "italic",hjust = 1) 
ggsave(filename='./coff_gph2_Aug22.png', plot=coff_gph , width = 6, height = 6, dpi = 600)


############ HEATMAP - patterns ##############


a<-read.table("/n/home02/rmehta/balskus/mindata/heatmap_pattern_min.txt",header=T,row.names=1,sep="\t",fill = TRUE, comment.char = "" , check.names = FALSE)
test_labels <- read.table("/n/home02/rmehta/balskus/mindata/heatmapqval_patt_min.txt",header=T,row.names=1,sep="\t",fill=TRUE, comment.char="",check.names=FALSE)
test_labels[is.na(test_labels)] <- ""

color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))
breaksList = seq(-.004, 0.004, by = 0.0002)

  
p<-pheatmap(
  mat               = a,
#  color             = inferno(10),
  border_color      = "black",
  show_colnames     = T,
  show_rownames     = T,
  drop_levels       = TRUE,
  fontsize          = 14,
 cluster_rows = FALSE, cluster_cols = FALSE,cellwidth = 14, cellheight = 14,kmeans_k = NA,clustering_method = "complete",display_numbers = test_labels,fontsize_number=14,
color=colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,legend_breaks = c(-.004, 0, 0.004), 
main = "", legend_labels = c("-0.004", "0", "0.004"))

ggsave(filename='./patternmin_aug22.png', plot=p, width = 5, height = 12, scale=2, dpi=600)

############ HEATMAP - foods ##############


a<-read.table("/n/home02/rmehta/balskus/mindata/foodtable_min_beta.txt",header=T,row.names=1,sep="\t",fill = TRUE, comment.char = "" , check.names = FALSE)
test_labels <- read.table("/n/home02/rmehta/balskus/mindata/foodtable_min_q.txt",row.names=1,header=T,sep="\t",fill=TRUE, comment.char="",check.names=FALSE)

color = colorRampPalette(rev(brewer.pal(n = 5, name = "RdBu")))
breaksList = seq(-.004, 0.004, by = 0.0002)

  
p<-pheatmap(
  mat               = a,
#  color             = inferno(10),
  border_color      = "black",
  show_colnames     = T,
  show_rownames     = T,
  drop_levels       = TRUE,
  fontsize          = 14,
 cluster_rows = FALSE, cluster_cols = FALSE,cellwidth = 14, cellheight = 14,kmeans_k = NA,clustering_method = "complete",display_numbers = test_labels,fontsize_number=14,
color=colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(length(breaksList)), # Defines the vector of colors for the legend (it has to be of the same lenght of breaksList)
         breaks = breaksList,legend_breaks = c(-.006, 0, 0.006), 
main = "", legend_labels = c("-0.006", "0", "0.006"))

ggsave(filename='./foodsmin_aug22.png', plot=p, width = 5, height = 12, scale=2, dpi=600)



############################################################################## checking for relationship with phenotypic data -- interaction ##################################################


####### checking for relationship with phenotypic data -- first take average of biomarker data ######

#read in plasma markers 

plasma <- read.csv('/n/home02/rmehta/balskus/mindata/metadata.csv', header=T) %>% select(c(idf,idfs,tc_plasma, hdlc_plasma,hba1cp, tg_plasma, crp_plasma, ldlc_plasma, METscore)) #### NEED TO REMEMBER WHERE THIS FILE CAME FROM 
idkey2 <- read.csv('/n/home02/rmehta/balskus/mindata/bugsdata.csv', header=T) %>% select(c(idf,id,idfs))
plasmaid <- inner_join(idkey2,plasma,by="idfs")

plasmaid$id <- as.factor(plasmaid$id)

plasmameta <- inner_join(plasmaid,gphmeta,by="id")
plasmameta$numdiff <- (as.numeric(substr(plasmameta$Sample.ID,13,13))-as.numeric(substr(plasmameta$idfs,7,7)))

plasmameta2 <-  plasmameta %>% filter(numdiff==4)
dim(plasmameta2)

plasmameta2$time <- substr(plasmameta2$Sample.ID,13,13)
plasmameta2<- plasmameta2%>% mutate(group=ifelse(as.numeric(time)<7,1,2))
plasmameta2$idgroup <- paste0(plasmameta2$id,"_",plasmameta2$group)
plasmameta2$time <- as.numeric(plasmameta2$time)
plasmameta2$id <- as.numeric(plasmameta2$id)

plasmameta3 <-plasmameta2%>% select(-c(ID1,samples,Sample.ID,run_accession,idfs,idf.y)) %>% group_by(as.factor(idgroup)) %>% summarise_at(vars(id:choc_avg), mean, na.rm = TRUE)

plasmameta3 <- data.frame(plasmameta3)
plasmameta3$id <- as.factor(plasmameta3$id)

median(plasmameta3$GpHcdh_mgx)
plasmameta3 <- plasmameta3 %>% mutate(gphmed = ifelse(GpHcdh_mgx>0.0008390005,1,0))
plasmameta3$gphmed <- as.factor(plasmameta3$gphmed)
plasmameta3 <- plasmameta3 %>% mutate(gphmed2 = ifelse(GpHcdh_mgx>0.0,1,0))
plasmameta3$gphmed2 <- as.factor(plasmameta3$gphmed2)

median(plasmameta3$eCadh_mgx)
plasmameta3 <- plasmameta3 %>% mutate(ecadmed = ifelse(eCadh_mgx>0.0005071445,1,0))
plasmameta3$ecadmed<- as.factor(plasmameta3$ecadmed)

median(plasmameta3$Hcdh_mgx)
plasmameta3 <- plasmameta3 %>% mutate(hcdhmed= ifelse(Hcdh_mgx>0.001378689,1,0))
plasmameta3$hcdhmed<- as.factor(plasmameta3$hcdhmed)


#continuous 
ccc1 <- lme(log(crp_plasma)~coffeeavg +asin(sqrt(GpHcdh_mgx))+ asin(sqrt(GpHcdh_mgx))*coffeeavg , random=~1|id, data=plasmameta3 ,na.action=na.omit,control=lmeControl(returnObject=TRUE))

#binary univariate 
ccc2 <- lme(log(crp_plasma)~coffeeavg +gphmed+ gphmed*coffeeavg , random=~1|id, data=plasmameta3 ,na.action=na.omit,control=lmeControl(returnObject=TRUE))

 
###### cru veg 

#continuous (p=0.06)
ddd1<- lme(log(crp_plasma)~cruvegavg+asin(sqrt(GpHcdh_mgx))+ asin(sqrt(GpHcdh_mgx))*cruvegavg+ agemlvs + bmi12 + abx_avg+bristol_avg+coffeeavg , random=~1|id, data=plasmameta3 ,na.action=na.omit,control=lmeControl(returnObject=TRUE))

#binary 
ddd2<- lme(log(crp_plasma)~cruvegavg+gphmed+ gphmed*cruvegavg+ agemlvs + bmi12 + abx_avg+bristol_avg+coffeeavg , random=~1|id, data=plasmameta3 ,na.action=na.omit,control=lmeControl(returnObject=TRUE))
#                      Value Std.Error  DF   t-value p-value
#(Intercept)       -5.807279 1.2505364 291 -4.643831  0.0000
#cruvegavg         -0.737216 0.2341448 291 -3.148549  0.0018
#gphmed1           -0.325484 0.1908371 158 -1.705562  0.0901
#agemlvs            0.055399 0.0156966 291  3.529337  0.0005
#bmi12              0.083453 0.0180933 291  4.612358  0.0000
#abx_avg            0.096256 0.1406141 291  0.684539  0.4942
#bristol_avg       -0.027835 0.0725444 291 -0.383701  0.7015
#coffeeavg         -0.048313 0.0465059 291 -1.038851  0.2997
#cruvegavg:gphmed1  0.839039 0.2966652 158  2.828235  0.0053



pcrp <-ggplot(plasmameta3 ) +
  aes(x = cruvegavg, y = log(crp_plasma), color = gphmed) +
  geom_point(alpha=0.4,size=2) +
  geom_smooth(method = "lm",se=FALSE,size=1.5)+theme_bw(base_size=24)+ scale_color_manual(breaks = c("0", "1"),
                        values=c("red", "dodgerblue2"))
#annotate("text", x = 2, y = (-2), size = 7, label = "pint = 0.005", fontface = "italic",hjust = 1) 


ggsave(filename='./pcrp_aug22.png', plot=pcrp, width = 5, height = 4, scale=2, dpi=600)


#specificity tests (for other genes, null) 

#binary #P > 0.12, 0.98
ddd3<- lme(log(crp_plasma)~cruvegavg+ecadmed+ ecadmed*cruvegavg+ agemlvs + bmi12 + abx_avg+bristol_avg+coffeeavg , random=~1|id, data=plasmameta3 ,na.action=na.omit,control=lmeControl(returnObject=TRUE)) #borderline sig 
ddd4<- lme(log(crp_plasma)~cruvegavg+hcdhmed+ hcdhmed*cruvegavg+ agemlvs + bmi12 + abx_avg+bristol_avg+coffeeavg , random=~1|id, data=plasmameta3 ,na.action=na.omit,control=lmeControl(returnObject=TRUE))


#specificity (other markers)

ddd16<- lme(log(hba1cp)~cruvegavg+gphmed+ gphmed*cruvegavg+ agemlvs + bmi12 + abx_avg+bristol_avg+coffeeavg , random=~1|id, data=plasmameta3 ,na.action=na.omit,control=lmeControl(returnObject=TRUE))

ddd8<- lme(log(tc_plasma)~cruvegavg+gphmed+ gphmed*cruvegavg+ agemlvs + bmi12 + abx_avg+bristol_avg+coffeeavg , random=~1|id, data=plasmameta3 ,na.action=na.omit,control=lmeControl(returnObject=TRUE))
ddd9<- lme(log(hdlc_plasma)~cruvegavg+gphmed+ gphmed*cruvegavg+ agemlvs + bmi12 + abx_avg+bristol_avg+coffeeavg , random=~1|id, data=plasmameta3 ,na.action=na.omit,control=lmeControl(returnObject=TRUE))
ddd11<- lme(log(tg_plasma)~cruvegavg+gphmed+ gphmed*cruvegavg+ agemlvs + bmi12 + abx_avg+bristol_avg+coffeeavg , random=~1|id, data=plasmameta3 ,na.action=na.omit,control=lmeControl(returnObject=TRUE))
ddd13<- lme(log(ldlc_plasma)~cruvegavg+gphmed+ gphmed*cruvegavg+ agemlvs + bmi12 + abx_avg+bristol_avg+coffeeavg , random=~1|id, data=plasmameta3 ,na.action=na.omit,control=lmeControl(returnObject=TRUE))

summary(ddd8) #P=0.08
summary(ddd9) #p=0.93
summary(ddd11) #p=0.30
summary(ddd13) #p=0.09
summary(ddd16) #p=0.28

#### make quintiles 

plasmameta3$cruquart<- cut(
  plasmameta3$cruveg,
  breaks = quantile(plasmameta3$cruveg, c(0, 0.25, 0.5, 0.75, 1)),
  labels = c("QQ1", "QQ2", "QQ3", "QQ4"),
  right  = FALSE,
  include.lowest = TRUE
)

plasmameta_low <- plasmameta3 %>% filter(gphmed==0)
plasmameta_hi <- plasmameta3 %>% filter(gphmed==1)

plasmameta_low %>% group_by(cruquart) %>% summarize(meancrp=mean(crp_plasma),sdcrp=sd(crp_plasma),na.rm=TRUE)
# A tibble: 4 x 4a),na.rm=TRUE)
#  cruquart meancrp sdcrp na.rm
#  <fct>      <dbl> <dbl> <lgl>
#1 QQ1        2.32   3.88 TRUE
#2 QQ2        2.29   3.90 TRUE
#3 QQ3        0.938  1.04 TRUE
#4 QQ4        0.887  1.18 TRUE


plasmameta_hi %>% group_by(cruquart) %>% summarize(meancrp=mean(crp_plasma),sdcrp=sd(crp_plasma),na.rm=TRUE)
# A tibble: 4 x 4),na.rm=TRUE)
#  cruquart meancrp sdcrp na.rm
#  <fct>      <dbl> <dbl> <lgl>
#1 QQ1         1.31  1.58 TRUE
#2 QQ2         2.24  4.58 TRUE
#3 QQ3         1.81  3.04 TRUE
#4 QQ4         2.19  3.70 TRUE

#ptrend

aaa2<- lme(log(crp_plasma)~cruvegavg+ agemlvs + bmi12 + abx_avg+bristol_avg+coffeeavg , random=~1|id, data=plasmameta_hi,na.action=na.omit,control=lmeControl(returnObject=TRUE))#0.54 hi
bbb2<- lme(log(crp_plasma)~cruvegavg+ agemlvs + bmi12 + abx_avg+bristol_avg+coffeeavg , random=~1|id, data=plasmameta_low,na.action=na.omit,control=lmeControl(returnObject=TRUE))#0.003 lo


##### test gordonibacteria specificity 

gordoni <-read.table("/n/home02/rmehta/balskus/mindata/MinSpecies_Sept_MGX.txt",header=T,sep="\t")%>% rename(Sample.ID=ID)
plasmametag<-inner_join(gordoni,plasmameta2,by="Sample.ID") %>% select(-c(ID1,ID_full,samples,Sample.ID,run_accession,idfs,idf.y)) %>% group_by(as.factor(idgroup)) %>% summarise_at(vars(Gordonibacter:choc_avg), mean, na.rm = TRUE)

median(plasmametag$Gordonibacter)
plasmametag <- plasmametag %>% mutate(gordmed = ifelse(Gordonibacter>0.0229425,1,0))
plasmametag$gordmed <- as.factor(plasmametag$gordmed)


ddd10<- lme(log(crp_plasma)~cruvegavg+gordmed + gordmed *cruvegavg+ agemlvs + bmi12 + abx_avg+bristol_avg+coffeeavg , random=~1|id, data=plasmametag,na.action=na.omit,control=lmeControl(returnObject=TRUE)) #p=0.23

