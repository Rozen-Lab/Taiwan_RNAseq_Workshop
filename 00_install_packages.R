# install all the required packages
install.packages("data.table")
install.packages('Matrix')
BiocManager::install("ggplot2")
install.packages('patchwork')
install.packages('gridExtra')
install.packages('RColorBrewer')
BiocManager::install("Seurat")
install.packages('clustree')
install.packages('pheatmap')
BiocManager::install("clusterProfiler")
install.packages('msigdbr')
install.packages('hdf5r')

# create the folder to store downloaded h5 files
if(!dir.exists("0data/")){dir.create("0data/")}
system("unzip bmRawCounts.zip")