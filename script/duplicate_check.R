


ibd_result = fread("tmp_duplicate/ibd_results.genome")
ibd_result%>%
  dplyr::select(IID1,IID2,PI_HAT)%>%
  arrange(desc(PI_HAT))%>%
  mutate(SAMPLE.ID1=IID1,
         SAMPLE.ID2=IID2,
         IID1=gsub(".+-","",IID1),
         IID2=gsub(".+-","",IID2)) -> ibd_result


# check IID consistency across data
ibd_result%>%
  filter(IID1==IID2)%>%
  filter(PI_HAT<0.3)

# check duplicated genotype
ibd_result%>%
  filter(IID1!=IID2)%>%#pull(PI_HAT)%>%hist
  filter(PI_HAT>0.9)

ibd_result%>%
  filter(IID1!=IID2)%>%#pull(PI_HAT)%>%hist
  filter(PI_HAT>0.9)%>%
  dplyr::select(SAMPLE.ID1:SAMPLE.ID2)%>%
  bind_rows(.,
            ibd_result%>%
              filter(IID1!=IID2)%>%#pull(PI_HAT)%>%hist
              filter(PI_HAT>0.9)%>%
              dplyr::select(SAMPLE.ID1:SAMPLE.ID2)%>%
              dplyr::mutate(SAMPLE.ID1_=SAMPLE.ID2,
                            SAMPLE.ID2_=SAMPLE.ID1)%>%
              mutate(SAMPLE.ID1=SAMPLE.ID1_,
                     SAMPLE.ID2=SAMPLE.ID2_)%>%
              dplyr::select(SAMPLE.ID1:SAMPLE.ID2))%>%
  mutate(Current_IID = gsub(".+-","",SAMPLE.ID1),
         Potential_Swap_IID = gsub(".+-","",SAMPLE.ID2))%>%
  dplyr::select(Current_IID:Potential_Swap_IID) -> res
  
res%>%fwrite(file = "combined_plink/problematic_samples.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)
res%>%fwrite(file = "combined_vcf/problematic_samples.txt",append = F,quote = F,sep = "\t",row.names = F,col.names = T)

