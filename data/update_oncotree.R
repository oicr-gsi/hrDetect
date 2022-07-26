install.packages("mskcc.oncotree")
library(mskcc.oncotree)

tumor_types <- as.data.frame(get_tumor_types(oncotree_version = "oncotree_2021_11_02"))[,c("oncotree_version","oncotree_code","tissue")]
tumor_types <- tumor_types[order(tumor_types$tissue),]

write.table(x = tumor_types,file = '~/oncotree_2021_11_02.csv',quote = F,row.names = F,sep = ",")
