\name{LXheatmap}
\alias{LXheatmap}
\title{LXheatmap:for the analysis of heatmap}
\usage{
LXheatmap(data_file,cutree_rows_n)
}
\description{
LXheatmap:for the analysis of heatmap
}
\examples{
install.packages("devtools")

library(devtools)

install_github("gluck4668/LXheatmap")

library(LXheatmap)

#---------------
data(data_heatmap01)
data(data_heatmap02)
#---------------

rm(list=ls())

setwd("D:/LXheatmap_example")

data_file <- c("data_heatmap02.xlsx")

cutree_rows_n=1

LXheatmap(data_file,cutree_rows_n)

}
