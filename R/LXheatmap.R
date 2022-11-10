
LXheatmap <- function(data_file,cutree_rows_n){

 all_packages <- data.frame(installed.packages())

 pack <- data.frame(c("plyr","ggplot2","pheatmap","openxlsx","stringr","psych","dplyr","ggplotify","conflicted"))

  pack$type <- pack[,1] %in% all_packages$Package

for (i in 1:nrow(pack)){if (!requireNamespace(pack[i,1], quietly=TRUE))
  install.packages(pack[i,1],update = F,ask = F)
 }
 rm(i)


 packages <- as.character(pack[,1])

 for(i in packages){
   library(i, character.only = T)
 }
 rm(i)

 conflict_scout()

 conflict_prefer("%+%", "psych")
 conflict_prefer("distance", "psych")

 conflict_prefer("filter", "dplyr")
 conflict_prefer("lag", "dplyr")
 conflict_prefer("arrange", "dplyr")
 conflict_prefer("select", "dplyr")
 conflict_prefer("summarise", "dplyr")
 conflict_prefer("summarize", "dplyr")
 conflict_prefer("rename", "dplyr")
 conflict_prefer("mutate", "dplyr")
 conflict_prefer("count", "dplyr")
 conflict_prefer("failwith", "dplyr")
 conflict_prefer("id", "dplyr")
 conflict_prefer("combine", "dplyr")
 conflict_prefer("desc", "dplyr")
 conflict_prefer("collapse", "dplyr")
 conflict_prefer("slice", "dplyr")

 conflict_prefer("Position", "ggplot2")

 #--------------------------------------------------------
df <- read.xlsx(data_file)

colnames(df)[1]="ID"

df <- data.frame(distinct(df, ID, .keep_all = TRUE))

rownames(df) <- df[,1]

df <- df[,-1]


groups_names <- data.frame(colnames(df))
groups_n <- nrow(groups_names)

if(tolower(trimws(groups_names[groups_n,]))==tolower(trimws("type"))){
    df_data <- df[,-groups_n]
    type <- data.frame(df[,groups_n])
    rownames(type) <- rownames(df)
    colnames(type) <- c("types") } else
    {df_data <- df
    type <- NULL}

df_data$sum <- rowSums(df_data) # sum  all the data of the rows

df_data <- df_data[df_data[,10]>0,] # screen out the rows without 0.

df_data <- df_data[,-10]

genes_names <- data.frame(rownames(df_data))

groups <- data.frame(colnames(df_data))

groups$names <- gsub("\\d+$","", groups[,1])

rownames(groups) <- groups[,1]

groups <- dplyr::select(groups,names)

if (nrow(df_data)<60)
  {tree_col=30
  tree_row=60} else
  {tree_col=0
  tree_row=0}

if (nrow(df_data)<50)
  show_row=TRUE else
    show_row=FALSE

fontsize_r <- dplyr::case_when(nrow(df_data)<20 ~10,
                               nrow(df_data)<30 ~9,
                               nrow(df_data)<50 ~8,
                               TRUE ~8
                               )
print("---------------------------------------------------------------------")
print("The analysis results can be found in the folder of <analysis result>")

heat_map <- pheatmap::pheatmap(df_data,scale ="row",
              cluster_rows = T,cluster_cols = F,
             # clustering_distance_rows = "euclidean",
              clustering_method = "ward.D2",
              main="Heatmap graphics", fontsize=14,
              show_rownames = show_row,show_colnames = T,
              angle_col=45,
              annotation_names_row = T,annotation_names_col = T,
              fontsize_row = fontsize_r,fontsize_col = 12,
              annotation_col = groups,
              annotation_row = type,
              treeheight_col = tree_col,treeheight_row = tree_row,
              cutree_rows = cutree_rows_n, cutree_cols =1,
              border_color = NA,
              border=FALSE,
              color=colorRampPalette(c("deepskyblue","white","red"))(100))

if(dir.exists("analysis result")==FALSE)
  dir.create("analysis result")
#ggsave("analysis result/Heatmap graphics 01.png",heat_map,width=1200, height =1000, dpi=180,units = "px")


as_ggplot <- as.ggplot(heat_map)
ggplot_heat <- as_ggplot+ theme(plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
ggplot_heat

ggsave("analysis result/Heatmap graphics.png",ggplot_heat,width=1200, height =1000, dpi=180,units = "px")


#----------------------------------------------------------------
row_cluster <- cutree(heat_map$tree_row,k=cutree_rows_n)
row_cluster

neworder <- heat_map$tree_row$order
neworder

neworder_df=df[neworder,]

neworder_df[,ncol(neworder_df)+1]=row_cluster[match(rownames(neworder_df),names(row_cluster))]
colnames(neworder_df)[ncol(neworder_df)]="Cluster"
neworder_df

write.xlsx(neworder_df,"analysis result/new_order_data.xlsx",rowNames=T,colNames=T)

}

