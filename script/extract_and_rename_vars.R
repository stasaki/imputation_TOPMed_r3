


library(data.table)
library(tidyverse)



array=fread("tmp_plink/combined_cleaned.bim")
divs=fread("tmp_plink/DIVERSE_cleaned.bim")
wgs=fread("AMPAD_WGS/tmp/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08.recalibrated_variants_hg38_dn.QCed1_filtered-updated.bim")

head(array)
head(divs)

array%>%
  dplyr::select(-V3)%>%
  inner_join(.,divs%>%
               dplyr::select(-V3),by=c("V1","V4","V5","V6")) -> ds


ds%>%
  inner_join(.,wgs%>%
               dplyr::select(-V3),by=c("V1","V4","V5","V6")) -> ds_all

head(ds_all)

ds_all%>%
  dplyr::select(V2.x)%>%
  fwrite(.,file = "tmp_plink/Array_keep_vars.txt",
         append = F,quote = F,sep = "\t",row.names = F,col.names = F)

ds_all%>%
  dplyr::select(V2.y)%>%
  fwrite(.,file = "tmp_plink/Diverse_keep_vars.txt",
         append = F,quote = F,sep = "\t",row.names = F,col.names = F)

ds_all%>%
  dplyr::select(V2)%>%
  fwrite(.,file = "tmp_plink/AMPAD_keep_vars.txt",
         append = F,quote = F,sep = "\t",row.names = F,col.names = F)

rm(ds)
rm(array)
rm(divs)
rm(wgs)

ampad = fread("tmp_plink/AMPAD_cleaned.bim")
array = fread("tmp_plink/combined_cleaned.bim")
diverse = fread("tmp_plink/DIVERSE_cleaned.bim")

indx=match(array$V2,ds_all$V2.x)
head(indx)
array$V2 = ds_all$V2[indx]

indx=match(diverse$V2,ds_all$V2.y)
head(indx)
diverse$V2 = ds_all$V2[indx]

fwrite(array,file = "tmp_plink/combined_cleaned.bim",append = F,quote = F,sep = "\t",row.names = F,col.names = F)
fwrite(diverse,file = "tmp_plink/DIVERSE_cleaned.bim",append = F,quote = F,sep = "\t",row.names = F,col.names = F)
