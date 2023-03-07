import re
import pysam

def get_transcriptome(file_path='../../data/chr11_transcriptome.fasta'):
    transcripts = {}
    with open(file_path, 'r') as f:
        cur_transcript = None
        for line in f:
            if line[0] == '>':
                cur_transcript = line[1:].rstrip()
                transcripts[cur_transcript] = ""
            else:
                transcripts[cur_transcript] += line.rstrip()
    return transcripts

