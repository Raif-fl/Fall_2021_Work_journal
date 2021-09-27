# import plastid
# data structure for mapping read alignments to genomic positions
from plastid import BAMGenomeArray, VariableFivePrimeMapFactory, \
                        GTF2_TranscriptAssembler, Transcript, ThreePrimeMapFactory
import numpy as np
import numpy
import pandas as pd
import warnings

# Define the path to our Bam files
PATH = "/Users/keeganflanagan/Desktop/Khanh_position/genomes_and_samples/Fly/"

# Create a new mapping function in order to use the three prime
# mappings created by RiboWaltz.
def variable_threeprime_map_function(alignments,segment,p_offsets):
        reads_out = []
        count_array = numpy.zeros(len(segment))
        for read in alignments:
            for length, offset in zip(p_offsets["length"],p_offsets["p_offset"]):
                if length != len(read.positions):
                    continue # skip read if offset is outside read boundaries

             # count offset 3' to 5' if the `segment` is on the plus-strand
             # or is unstranded
                if segment.strand == "+":
                    p_site = read.positions[-offset - 1]
                elif segment.strand == ".":
                    p_site = read.positions[-offset - 1]
             # count offset from other end if `segment` is on the minus-strand
                elif segment.strand == "-":
                    p_site = read.positions[offset]

                if p_site >= segment.start and p_site < segment.end:
                    reads_out.append(read)
                    count_array[p_site - segment.start] += 1
        return reads_out, count_array

def VariableThreePrimeMapFactory(p_offsets):
    def new_func(alignments,segment):
        return variable_threeprime_map_function(alignments,segment,p_offsets=p_offsets)

        return new_func

# Load in the table of P-site offsets.
p_offsets=pd.read_csv("control_RPF_2_p_offsets.txt", sep="\t")

# load the transcript annotations from the GTF file.
# GTF2_TranscriptAssembler returns an iterator, so here we convert it to a list.
transcripts = list(GTF2_TranscriptAssembler(open(PATH + "Drosophila_melanogaster.BDGP6.32.103.gtf"),return_type=Transcript))

# Remove non-protein coding transcripts from transcripts list.
protein_coding = []
for transcript in transcripts:
    if transcript.attr['transcript_biotype'] == 'protein_coding':
        protein_coding.append(transcript)
