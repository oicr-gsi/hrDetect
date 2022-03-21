library(data.table)
library(tidyverse)

setwd('/Volumes/')

test1 <- fread('cgi/scratch/fbeaudry/sigTools_test/sigTools.test1.text')


ggplot(test1 %>% filter(V4 == "HRD_point"),aes(x=V2,y=as.numeric(V5))) + 
  geom_point() + theme_bw() + labs(x="Sample",y="Probability of HRD")
