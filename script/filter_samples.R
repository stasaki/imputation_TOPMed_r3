
setwd("/data/shinya/shinya/Project/imputation_TOPMed_r3")
library(tidyverse)
library(data.table)
fread("AMPAD_WGS/tmp/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08.recalibrated_variants_hg38_dn.QCed1_filtered-updated-chr1.fam")%>%
  mutate(type="AMPAD_WGS")%>%
  bind_rows(.,
            read.table("Diverse_Joint/tmp/DIVERSE.joint.callset_filtered-updated-chr1.fam")%>%
              mutate(type="Diverse"))%>%
  bind_rows(.,read.table("combined_plink/TOPMed_r3_array_hg38.fam")%>%
              mutate(type="Array"))%>%
  group_by(V1)%>%
  mutate(flag=duplicated(V1))%>%
  filter(flag)%>%
  dplyr::select(V1,V2,type)%>%
  mutate(V1=as.character(V1)%>%str_pad(.,8,pad = "0"),
         V2=as.character(V2)%>%str_pad(.,8,pad = "0"))%>%
  group_by(type)%>%
  nest() -> ds

ds$data[[1]]%>%
  write.table(.,file="tmp_plink/Diverse_remove.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)

ds$data[[2]] %>%
  write.table(.,file="tmp_plink/Array_remove.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)



ds$data[[1]]%>%mutate(V3=paste0(V1,"_",V2))%>%dplyr::select(V3)%>%
  write.table(.,file="tmp/Diverse_remove.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = F)

ds$data[[2]]%>%mutate(V3=paste0(V1,"_",V2))%>%dplyr::select(V3)%>%
  write.table(.,file="tmp/Array_remove.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = F)




fread("AMPAD_WGS/tmp/NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08.recalibrated_variants_hg38_dn.QCed1_filtered-updated-chr1.fam")%>%
  mutate(type="AMPAD_WGS")%>%
  bind_rows(.,
            read.table("Diverse_Joint/tmp/DIVERSE.joint.callset_filtered-updated-chr1.fam")%>%
              mutate(type="Diverse"))%>%
  bind_rows(.,
            read.table("BU_genotyping/tmp/BU_measured_uniq_pos_QCed_hg38_filtered-updated-chr1.fam")%>%
              mutate(type="BU"))%>%
  bind_rows(.,
            read.table("San1_Broad_1709/tmp/Merge.phaseIV.final_uniq_pos_QCed_hg38_filtered-updated-chr1.fam")%>%
              mutate(type="Broad"))%>%
  bind_rows(.,
            read.table("San1_CHOP_382/tmp/chop.rosmap.euam.vFinal_uniq_pos_QCed_hg38_filtered-updated-chr1.fam")%>%
              mutate(type="CHOP"))%>%
  bind_rows(.,
            read.table("San1_CHOP_RMM_AA/tmp/chop.rosmap.afam.vFinal_uniq_pos_QCed_hg38_filtered-updated-chr1.fam")%>%
              mutate(type="CHOP_RMM_AA"))%>%
  dplyr::select(V1,type)%>%
  mutate(V1=as.character(V1)%>%str_pad(.,8,pad="0")) -> array_type

ds

array_type%>%
  group_by(type)%>%
  nest() -> array_type
  

  array_type%>%
    filter(!type%in%c("AMPAD_WGS","Diverse"))%>%
    unnest()%>%
    filter(!V1%in%(  ds%>%
                      filter(type=="Array")%>%
                      unnest()%>%pull(V1)))%>%
    bind_rows(.,  array_type%>%
                filter(type%in%c("Diverse"))%>%
                unnest()%>%
                filter(!V1%in%(  ds%>%
                                   filter(type=="Diverse")%>%
                                   unnest()%>%pull(V1))))%>%
    bind_rows(.,  array_type%>%
                filter(type%in%c("AMPAD_WGS"))%>%
                unnest()) -> array_type
  

  

fread("combined_plink/TOPMed_r3_array_and_WGS_pca.eigenvec")%>%
  mutate(V1=as.character(V1)%>%str_pad(.,8,pad="0"))%>%
  mutate(V2=as.character(V2)%>%str_pad(.,8,pad="0"))%>%
  left_join(.,array_type,by="V1") -> pca_eigenvec

colnames(pca_eigenvec) = c("IID","FID",paste0("PC",1:20),"Type")
  
fwrite(pca_eigenvec,file="combined_plink/TOPMed_r3_array_and_WGS_pca.eigenvec.txt.gz",append = F,quote = F,sep = "\t",row.names = F,col.names = T)
  
  
fread("combined_plink/TOPMed_r3_all_hg38_pca.eigenvec")%>%
  mutate(V1=as.character(V1)%>%str_pad(.,8,pad="0"))%>%
  mutate(V2=as.character(V2)%>%str_pad(.,8,pad="0"))%>%
  left_join(.,array_type,by="V1") -> pca_eigenvec

colnames(pca_eigenvec) = c("IID","FID",paste0("PC",1:20),"Type")

fwrite(pca_eigenvec,file="combined_plink/TOPMed_r3_all_hg38_pca.eigenvec.txt.gz",append = F,quote = F,sep = "\t",row.names = F,col.names = T)


pca_eigenvec%>%
  ggplot(.,aes(PC1,PC2))+geom_point()

