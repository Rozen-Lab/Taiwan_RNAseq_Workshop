###### Section 1
install.packages("devtools")
library(devtools)
install_github("kendomaniac/rCASC")

# downloading the required docker containers
library(rCASC)
downloadContainers(group="docker", containers.file="full")

###### Section 2.2 10X Genomics
home <- getwd()
system(paste("mkdir ", home, "/data", sep=""))
system(paste("mkdir ", home, "/data/scratch", sep=""))
system(paste("mkdir ", home, "/data/genomes", sep=""))
system(paste("mkdir ", home, "/data/genomes/cellranger_hg19mm10", sep=""))
system(paste("mkdir ", home, "/data/fastqs", sep=""))

#download refrence genome
setwd(paste(home, "/data/genomes/cellranger_hg19mm10", sep=""))
#getting the human and mouse cellranger index
system("wget http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-hg19-and-mm10-2.1.0.tar.gz")
# must untar this
system("tar -xvf refdata-cellranger-hg19-and-mm10-2.1.0.tar.gz")

# download sample fastq files
setwd(paste(home, "/data/fastqs", sep=""))
# downloading 100 cells 1:1 Mixture of Fresh Frozen Human (HEK293T) and Mouse (NIH3T3) Cells
system("wget http://cf.10xgenomics.com/samples/cell-exp/2.1.0/hgmm_100/hgmm_100_fastqs.tar")
system("tar xvf hgmm_100_fastqs.tar")
# The cellranger analysis is run without the generation of the secondary analysis

setwd(home)
cellrangerCount(group="docker",  transcriptome.folder=paste(home, "data/genomes/cellranger_hg19mm10", sep=""),  
                fastq.folder=paste(home, "data/fastqs", sep=""),  expect.cells=100, 
                nosecondary=TRUE, scratch.folder=paste(home, "data/scratch", sep=""))

###### Section 2.3 Spatial Transcriptomics
library(rCASC)
home <- getwd()
system(paste("mkdir ", home, "/data/V1_Mouse_Kidney_fastqs", sep=""))
system(paste("mkdir ", home, "/data/genomes/cellranger_mm10", sep=""))
system(paste("mkdir ", home, "/data/spatial_scratches", sep=""))
system(paste("mkdir ", home, "/data/spatial_outputs", sep=""))

setwd(paste(home, "/data/V1_Mouse_Kidney_fastqs", sep=""))
system("wget http://s3-us-west-2.amazonaws.com/10x.files/samples/spatial-exp/1.0.0/V1_Mouse_Kidney/V1_Mouse_Kidney_fastqs.tar")
system("tar xvf V1_Mouse_Kidney_fastqs.tar")
# DatasetImage
system("wget http://cf.10xgenomics.com/samples/spatial-exp/1.0.0/V1_Mouse_Kidney/V1_Mouse_Kidney_image.tif")

# referenceGenomeMM10
setwd(paste(home, "/data/genomes/cellranger_mm10", sep=""))
system("wget http://cf.10xgenomics.com/supp/spatial-exp/refdata-cellranger-mm10-3.0.0.tar.gz")
system("tar xvf refdata-cellranger-mm10-3.0.0.tar.gz")

# Create count matrix from spatial transcriptomics fasta
setwd(home)
stpipeline(group="docker", scratch.folder=paste(home, "/data/spatial_scratches", sep=""), 
           data.folder=paste(home, "/data/spatial_outputs", sep=""), #result folder
           genome.folder=paste(home, "/data/genomes/cellranger_mm10/refdata-cellranger-mm10-3.0.0", sep=""), 
           fastqPathFolder=paste(home, "/data/V1_Mouse_Kidney_fastqs/V1_Mouse_Kidney_fastqs", sep=""), 
           ID="hey",imgNameAndPath=paste(home, "/data/V1_Mouse_Kidney_fastqs/V1_Mouse_Kidney_image.tif", sep=""),
           slide="V19L29-096",area="B1")

# analyzes the data that came up from permutationClustering script with spatial transcriptomics data
spatialAnalysis2( group = "docker", scratch.folder=paste(home, "/data/spatial_scratches", sep=""), 
                  file= paste(home, "/data/spatial_outputs", "/HBC_BAS1_expr-var-ann_matrix.csv", sep=""),
                  nCluster=9, separator=",", 
                  tissuePosition= paste(home, "/data/spatial_outputs", "spatial/tissue_positions_list.csv", sep=""),
                  Sp = 0.8, percentageIncrease = 10
)

###### Section 3 Counts matrix manipulation
# Subset of a mouse single cell experiment made to quickly test Lorenz filter.
setwd(home)
system("wget http://130.192.119.59/public/testSCumi_mm10.csv.zip")
unzip("testSCumi_mm10.csv.zip")
tmp <- read.table("testSCumi_mm10.csv", sep=",", header=T, row.names=1)
dim(tmp)
#27998   806
write.table(tmp, "testSCumi_mm10.txt", sep="\t", col.names=NA)
filterZeros(file=paste(getwd(),"testSCumi_mm10.txt",sep="/"), threshold=0, sep="\t")
#Out of 27998 genes 11255 are left after removing genes with no counts
#output is filtered_testSCumi_mm10.txt
