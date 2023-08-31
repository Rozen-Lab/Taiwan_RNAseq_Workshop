# install all the required packages
install.packages("data.table")
install.packages('Matrix')
# new
install.packages("BiocManager")

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

# For chapter 4
install.packages("pdist")
install.packages("phate")

if(!dir.exists("0data/")){
  # presumably we need to download the data
  # create the folder to store downloaded h5 files
  dir.create("0data/")
  # get the data and put the unzipped data in 0data
  # The following calls to system do not work on Windows
  system("wget -O 0data/bmRawCounts.zip https://www.dropbox.com/s/7cs0mhdlbu1e437/bmRawCounts.zip?dl=0")  
  system("cd 0data; unzip bmRawCounts.zip")
}
