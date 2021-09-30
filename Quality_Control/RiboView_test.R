##### Installing riboview #####
#install.packages("reticulate")
#install.packages("Rtsne")
#install.packages("~/Desktop/RiboVIEWreticulate_2.1.tar.gz", repos = NULL, type ="source")

##### Loading up packages #####
library("devtools")
library("reticulate")
library(RiboVIEWreticulate)
library(rminiconda)

# Used to install pysam in the miniconda environment used downloaded by reticulate,
# Not the py2 version used by rminiconda.
#py_install("pysam")
#py_install("biopython==1.76")


##### Configure reticulate to work with python 2 #####
#remotes::install_github("hafen/rminiconda")
#rminiconda::install_miniconda(version = 2, name = "py2")
py2 = rminiconda::find_miniconda_python("py2")
reticulate::use_python(py2, required = TRUE)
use_condaenv(py2)
#rminiconda_pip_install("pysam","py2")
#rminiconda_pip_install("biopython==1.76","py2")

##### Set Up For Quality Control #####
# Defining path to reference sequence, GTF file, and CDS table
refFASTA = "/Users/keeganflanagan/Desktop/Khanh_position/genomes_and_samples/Fly/BDGP6.32.103.cdna.all.fa"
refGTF = "/Users/keeganflanagan/Desktop/Khanh_position/genomes_and_samples/Fly/Drosophila_melanogaster.BDGP6.32.103.gtf"
refCDS = "/Users/keeganflanagan/Desktop/Khanh_position/genomes_and_samples/Fly/dmel_CDStable.tsv"
pathout = "/Users/keeganflanagan/Desktop/Khanh_position/RiboView_test/resu/"

# Creating CDS table
gtf2table(refGTF, refCDS)

# Define the address and name of aligned sequences (in BAM format) of each sample,
reads_c1_r1 = "/Users/keeganflanagan/Desktop/Khanh_position/genomes_and_samples/Transcriptome_dmel/control_RPF_2_Aligned.toTranscriptome.out.bam"
reads_c2_r1 = "/Users/keeganflanagan/Desktop/Khanh_position/genomes_and_samples/Transcriptome_dmel/Fmr1_RPF_2_Aligned.toTranscriptome.out.bam"

# Define the list of files
list.bam = list(reads_c1_r1, reads_c2_r1)

# Names for the conditions of the experiment
XP.conditions = c("Control","Fmr1")
XP.conditions.i = c(1,2)
XP.names = c("Ctrl.R1","Fmr1.R1")

##### QC Preliminary Calculations #####

####!!!!! Important! In order to run this command properly you need to go into the 
#!!!!! Source code for periodicity_line 47 and change all of the int(element) calls
# !!!!! to int(float(element)). You will need to reload the library after you make this change. 
# !!!! Hmmm, this txt file ends up totally empty. 

# !!!!!! Note that you must use an alignment to the Transcriptome, NOT an 
# Alignment to the genome. 
#periodicity(list.bam, refCDS, refFASTA, pathout, XP.names, 
#            versionStrip = FALSE, python.messages = TRUE, mitochondrion = FALSE)

# What on earth is montage and why is it so evil? 

# Review and select adequate footprint length via the following command:
attach(listminmax <- select.FPlen(list.bam, pathout, XP.names))

# The previous command simply does not work, so I looked at the graphs manually and
# chose mini and maxi myself. 
mini=25
maxi=32

# Well this does just does not work at all. I need to figure out why this is 
# not calling to the right version of python. Maybe it is still trying to use
# That darn thing that reticulate downloaded for itself. must learn how to 
# install pysam with reticulate.

# See if you can figure out how to add a new conda environment to your conda
# environment path, might make it easier to configure the one installed 
# by reticulate.

# AHHHAHAHAHAHAHAHAHA. It also has the float int problem. 

# I had to replace all // with / to avoid generating floats from division. 

# Slightly altered the location of peak list in the hopes that it might improve something.

# Note that setting python.messages to true is a terrible idea for this one.
# The output will go on forever. 

# must switch all print to print()
enrichmentNoccupancy(list.bam, refCDS, refFASTA, mini, maxi, XP.names,
                     pathout, versionStrip = FALSE, python.messages=FALSE,
                     mitochondrion = FALSE)

# Note to sef: This trace and untrace thing is really useful for altering packages.
# I should remember to learn how to use this, especially as I get better with R. 

# Perhaps this whole rminiconda thing was a mistake, I am sure reticulate has 
# built in functions that can do this. Let me try to downgrade the python version in
# there, then I might be able to make this work. 

