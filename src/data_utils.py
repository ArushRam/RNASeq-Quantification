import re
import pysam

class data_utils:
    def __init__(self, transcriptome_path, alignment_path):
        self.transcriptome = None
        self.transcriptome_path = transcriptome_path
        self.alignment_path = alignment_path
        self.alignment_file = pysam.AlignmentFile(self.alignment_path, 'rb')
    
    def get_transcriptome(self):
        '''
        Build transcriptome dictionary (name:sequence) if not already built, and return it.
        '''
        if self.transcriptome is not None:
            return self.transcriptome
        
        # Build Transcript Dictionary from FASTA file
        self.transcriptome = {}
        with open(self.transcriptome_path, 'r') as f:
            cur_transcript = None
            for line in f:
                if line[0] == '>':
                    cur_transcript = line[1:].rstrip()
                    self.transcriptome[cur_transcript] = ""
                else:
                    self.transcriptome[cur_transcript] += line.rstrip()
        return self.transcriptome
    
    def get_transcript_lengths(self):
        lengths = {}
        if not self.transcriptome:
            self.get_transcriptome()
        for transcript in self.transcriptome:
            lengths[transcript] = len(self.transcriptome[transcript])
        return lengths
    
    def get_read_count(self):
        '''
        Count total number of reads in BAM file
        '''
        return self.alignment_file.count(until_eof=True)
    
    def get_alignment_iterator(self):
        return self.alignment_file.fetch(until_eof=True)
    
    def bam_to_fastq(self, fname):
        '''
        Convert reads from BAM alignment to FASTQ (single-end read data) for use in Kallisto comparison.
        One entry per READ, not alignment
        '''
        reads = self.get_alignment_iterator()
        seen_reads = set()
        read_re = 'read([0-9]*)'
        with open(fname, 'w') as fastq:
            for read in reads:
                name = re.search(read_re, read.query_name).group(0)
                if name in seen_reads:
                    continue
                seen_reads.add(name)
                seq = read.get_forward_sequence()
                qual = read.qual
                if read.is_reverse:
                    qual = qual[::-1]
                fastq.write(f'@{name}\n{seq}\n+\n{qual}\n')
        print('Total Reads: ', len(seen_reads))
