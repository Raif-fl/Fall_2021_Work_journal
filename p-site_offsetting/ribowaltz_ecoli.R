library(GenomicFeatures)

gfftx = makeTxDbFromGFF(file = "/Users/keeganflanagan/Desktop/Khanh_position/genomes_and_samples/Ecoli/GCF_015291845.1_ASM1529184v1_genomic.gff",
                format = "gff3",
                dataSource = "gff3 annotation file for E.coli assembly GCA_000005845.2",
                organism = "Escherichia coli")
# Load the ribowaltz package
library(riboWaltz)

# Load in the fly gff3 annotation file
df = create_annotation(txdb = gfftx,
                       dataSource = "gff3 annotation file for E.coli assembly GCA_000005845.2",
                       organism = "Escherichia coli")

# Ribowalts requires a file that has been aligned to the transcriptome rather than 
# a file that has been aligned to the genome. 

# Load up bam_list from transcriptome
bam_list = bamtolist("/Users/keeganflanagan/Desktop/Khanh_position/genomes_and_samples/Transcriptome_ecol",
                     annotation = df)

# Calculate P_site offset
offsets = psite(data = bam_list, start = FALSE, extremity = "3end",
                plot = TRUE, plot_dir = "/Users/keeganflanagan/Desktop/Khanh_position/Data/P-site_offsets")

# Filter the offsets dataframe to only include the information needed by plastid. 
plastid_offsets = data.frame(length=offsets$length,
                             p_offset=offsets$corrected_offset_from_3)

write.table(plastid_offsets, "/Users/keeganflanagan/Desktop/Khanh_position/Data/P-site_offsets/Fmr1_RPF_2_P_Offsets",sep = "\t" , row.names = FALSE)

