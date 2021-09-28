##### Installing riboview #####
#install.packages("reticulate")
#install.packages("Rtsne")
#install.packages("~/Desktop/RiboVIEWreticulate_2.1.tar.gz", repos = NULL, type ="source")

##### Loading up packages #####
library("devtools")
library("reticulate")
library(RiboVIEWreticulate)
library(rminiconda)

##### Configure reticulate to work with python 2 #####
#remotes::install_github("hafen/rminiconda")
#rminiconda::install_miniconda(version = 2, name = "py2")
py2 = rminiconda::find_miniconda_python("py2")
reticulate::use_python(py2, required = TRUE)
rminiconda_pip_install("pysam","py2")
rminiconda_pip_install("biopython==1.76","py2")

##### Set Up For Quality Control #####
# Defining path to reference sequence, GTF file, and CDS table
refFASTA = "/Users/keeganflanagan/Desktop/Khanh_position/genomes_and_samples/Fly/BDGP6.32.cdna.fa"
refGTF = "/Users/keeganflanagan/Desktop/Khanh_position/genomes_and_samples/Fly/Drosophila_melanogaster.BDGP6.32.103.gtf"
refCDS = "/Users/keeganflanagan/Desktop/Khanh_position/genomes_and_samples/Fly/dmel_CDStable.tsv"
pathout = "/Users/keeganflanagan/Desktop/Khanh_position/RiboView_test/resu/"

# Creating CDS table
gtf2table(refGTF, refCDS)

# Define the address and name of aligned sequences (in BAM format) of each sample,
reads_cl_2 = "/Users/keeganflanagan/Desktop/Khanh_position/genomes_and_samples/Fly/dmel_control_RPF_2_Aligned.sortedByCoord.out.bam"

# Define the list of files
list.bam = list(reads_c1_2)

# Names for the conditions of the experiment
XP.conditions = c("Control")
XP.conditions.i = c(1)
XP.names = c("Ctrl.R1")

##### QC Preliminary Calculations #####

####!!!!! Important! In order to run this command properly you need to go into the 
#!!!!! Source code for periodicity_line 47 and change all of the int(element) calls
# !!!!! to int(float(element)). You will need to reload the library after you make this change. 

# !!!! You also are going to have to fix the plotting function somehow.
# !!!! I think the plastid tutorial on plotting might have some suggestions. 
periodicity(list.bam, refCDS, refFASTA, pathout, XP.names, 
            versionStrip = FALSE, python.messages = TRUE, mitochondrion = FALSE)


