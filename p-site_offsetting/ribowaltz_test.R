# Load the ribowaltz package
library(riboWaltz)

# Load in the fly GTF annotation file
df = create_annotation(gtfpath = "/Users/keeganflanagan/Desktop/Khanh_position/genomes_and_samples/Fly/Drosophila_melanogaster.BDGP6.32.103.gtf")

# Ribowalts requires a file that has been aligned to the transcriptome rather than 
# a file that has been aligned to the genome. 

# Load up bam_list from transcriptome
bam_list = bamtolist("/Users/keeganflanagan/Desktop/Khanh_position/genomes_and_samples/Transcriptome_dmel",
                     annotation = df)

# Calculate P_site offset
offsets = psite(data = bam_list, start = FALSE, extremity = "3end",
                plot = TRUE, plot_dir = "/Users/keeganflanagan/Desktop/Khanh_position")

# Filter the offsets dataframe to only include the information needed by plastid. 
plastid_offsets = data.frame(length=offsets$length,
                             p_offset=offsets$corrected_offset_from_3)

write.table(plastid_offsets, "/Users/keeganflanagan/Desktop/Khanh_position/Fmr1_RPF_2_P_Offsets",sep = "\t" , row.names = FALSE)

