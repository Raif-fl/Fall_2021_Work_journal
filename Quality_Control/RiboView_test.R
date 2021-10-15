##### Installing riboview #####
#install.packages("reticulate")
#install.packages("Rtsne")
#install.packages("~/Downloads/RiboVIEWreticulate_2.1.tar.gz", repos = NULL, type ="source")

##### Loading up packages #####
library("devtools")
library("reticulate")
library(RiboVIEWreticulate)
library('tseriesChaos')

##### Creating a proper conda environemnt #####
#conda_create(envname = "reticulate_py2", 
#             packages= c("pysam", "biopython==1.76"), python_version = "2.7")

# Used to install pysam in the miniconda environment used downloaded by reticulate,
# Not the py2 version used by rminiconda.
#py_install("pysam")
#py_install("biopython==1.76")


##### Configure reticulate to work with python 2 #####
use_condaenv('reticulate_py2', required = TRUE)

##### Set Up For Quality Control #####
# Defining path to reference sequence, GTF file, and CDS table
refFASTA = "/home/keeganfl/Desktop/Work_Fall_2021/genomes_&_samples/dmel/Drosophila_melanogaster.BDGP6.32.cds.all.fa"
refGTF = "/home/keeganfl/Desktop/Work_Fall_2021/genomes_&_samples/dmel/Drosophila_melanogaster.BDGP6.32.103.gtf"
refCDS = "/home/keeganfl/Desktop/Work_Fall_2021/genomes_&_samples/dmel/dmel_CDStable.tsv"
pathout = "/home/keeganfl/Desktop/Work_Fall_2021/RiboView_test/resu/"

# Creating CDS table
gtf2table(refGTF, refCDS)

# Define the address and name of aligned sequences (in BAM format) of each sample,
reads_c1_r1 = '/home/keeganfl/Desktop/Work_Fall_2021/genomes_&_samples/tra_dmel/control_RPF_2_Aligned.toTranscriptome.out.bam'
reads_c2_r1 = '/home/keeganfl/Desktop/Work_Fall_2021/genomes_&_samples/tra_dmel/Fmr1_RPF_2_Aligned.toTranscriptome.out.bam'

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
periodicity(list.bam, refCDS, refFASTA, pathout, XP.names, 
            versionStrip = FALSE, python.messages = TRUE, mitochondrion = FALSE)

# What on earth is montage and why is it so evil? 

# Review and select adequate footprint length via the following command:
attach(listminmax <- select.FPlen(list.bam, pathout, XP.names))

# The previous command simply does not work, so I looked at the graphs manually and
# chose mini and maxi myself. 
mini=25
maxi=34

# Well this does just does not work at all. I need to figure out why this is 
# not calling to the right version of python. 

# I just commented out line 96 cause it does nothing but create an error. 

# It also has the float int problem.

# I had to replace all // with / to avoid generating floats from division. 

# Slightly altered 1the location of peak list in the hopes that it might improve something.
# and by this I mean I altered lines 688-700

# changed lines 1150 to 1160 to fix referenced before assignment glitch.

# Note that setting python.messages to true is a terrible idea for this one.
# The output will go on forever. 

# must switch all print to print()
enrichmentNoccupancy(list.bam, refCDS, refFASTA, mini, maxi, XP.names,
                     pathout, versionStrip = FALSE, python.messages=FALSE,
                     mitochondrion = FALSE)

generate.m.s(XP.conditions,XP.names,pathout,B=1000)

# Does not seem to be working, must investigate
visu.m.s.enrichmnt.res <- visu.m.s.enrichmnt(XP.conditions, XP.names, pathout)

visu.tracks.res <- visu.tracks(XP.conditions,
              XP.names,
              pathout,
              refCDS, # Why and how did they forget this input?
              mRNA="random", codon.labels=FALSE, codon.col="darkslateblue")

Venn.all.res <- Venn.all(XP.names, XP.conditions, pathout)

enricht.aroundA.res <- enricht.aroundA(XP.conditions, XP.names, pathout)

repl.correl.counts.Venn.res <- repl.correl.counts.Venn(XP.conditions, XP.names,
                                                         pathout)

#Thew an error
repl.correl.gene.res <- repl.correl.gene(XP.conditions, XP.names, pathout)

#Again, getting an error, something about width
repl.correl.codon.res <- repl.correl.codon(list.bam, refCDS, refFASTA,
                                             mini, maxi,
                                             XP.names, XP.conditions, pathout)

repl.correl.heatmap.res <- repl.correl.heatmap(XP.conditions.i, XP.names, pathout)

chx.artefacts.res <- chx.artefacts(XP.conditions, XP.names, pathout)

ntcodon.freq.nt.res <- ntcodon.freq.nt(XP.conditions, XP.names, pathout)

ntcodon.freq.cod.res <- ntcodon.freq.cod(XP.conditions, XP.names, pathout)

batch.effects.lm.e.res <- batch.effects.lm.e(XP.conditions, XP.names, pathout)

batch.effects.pca.res <- batch.effects.pca(XP.conditions, XP.names, pathout)

metagene.res <- metagene.all(XP.conditions, XP.names, pathout)

outputQc(pathout, XP.conditions)

outputMine(pathout, XP.conditions)


