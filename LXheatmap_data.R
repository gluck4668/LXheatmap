
setwd("D:/R-lin study/R-packages/LXheatmap")
library(openxlsx)
data_heatmap01 <- read.xlsx("data_heatmap01.xlsx")
data_heatmap02 <- read.xlsx("data_heatmap02.xlsx")

usethis::use_data(data_heatmap01,overwrite = T)
usethis::use_data(data_heatmap02,overwrite = T)

rm(list=ls())

data(data_heatmap01)
data(data_heatmap02)
