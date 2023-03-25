from src.data_utils import data_utils
import numpy as np
import pandas as pd
import pickle
import os
import re
pd.options.display.float_format = lambda x : '{:.0f}'.format(x) if int(x) == x else '{:,.2f}'.format(x)

class IsoEM:
    def __init__(self, transcriptome_path, alignment_path, mean_fragment_length=200, max_iter=500, convergence_threshold=0.01):
        self.data = data_utils(transcriptome_path, alignment_path)
        transcriptome = self.data.get_transcriptome()
        self.n_transcripts = len(transcriptome)
        self.mu = mean_fragment_length
        self.convergence_threshold = convergence_threshold
        self.max_iter = max_iter

        # Build numpy array and dictionary indexing transcripts
        self.transcript_ids = np.fromiter(transcriptome.keys(), dtype='U15')
        self.transcript_indices = {tid:i for i, tid in enumerate(transcriptome.keys())}

        # Transcript Lengths
        self.transcript_lengths = np.asarray([len(transcriptome[tid]) for tid in transcriptome.keys()])

        # Result Arrays
        self.theta = np.ones(self.n_transcripts)/self.n_transcripts
        self.expected_counts = None
        self.median_relative_differences = []
        self.mean_relative_differences = []

        # Get Weights
        self.get_weights()

    def get_weights(self):
        if os.path.isfile('data/weights.pickle'):
            with open('data/weights.pickle', 'rb') as handle:
                self.weights = pickle.load(handle)
            return
        
        alignments = self.data.get_alignment_iterator()
        self.weights = {}
        read_re = 'read([0-9]*)'
        for alignment in alignments:
            read_id, ref = int(re.search(read_re, alignment.query_name).group(1)), alignment.reference_name
            refidx = self.transcript_indices[ref]
            if read_id not in self.weights:
                self.weights[read_id] = {}
            if ref in self.weights[read_id]:
                self.weights[read_id][refidx] += 1
            else:
                self.weights[read_id][refidx] = 1

        with open('data/weights.pickle', 'wb') as handle:
            pickle.dump(self.weights, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
    def get_median_relative_difference(self, theta):
        eps = 1e-9
        relative_error = np.abs(self.theta-theta)/(self.theta + eps)
        median = np.median(relative_error)
        mean = np.mean(relative_error)
        return median, mean
    
    def normalize(self, read_counts, transcript_lengths, mu):
        return read_counts/(transcript_lengths - mu + 1) * (transcript_lengths >= mu)
    
    def estep(self):
        read_counts = np.zeros(self.n_transcripts)
        for read in self.weights:
            weight_sum = sum([self.theta[j] * self.weights[read][j] for j in self.weights[read]])
            for j in self.weights[read]:
                read_counts[j] += (self.theta[j] * self.weights[read][j])/weight_sum
        self.expected_counts = read_counts
        return read_counts
    
    def mstep(self):
        length_normalized_counts = self.normalize(self.expected_counts, self.transcript_lengths, self.mu)
        new_theta = length_normalized_counts/np.sum(length_normalized_counts)
        return new_theta
    
    def EM(self):
        med_rel_diff, epoch = 1, 0
        while med_rel_diff > self.convergence_threshold and epoch < self.max_iter:
            self.estep()
            new_theta = self.mstep()
            med_rel_diff, mean_rel_diff = self.get_median_relative_difference(new_theta)
            self.median_relative_differences.append(med_rel_diff)
            self.mean_relative_differences.append(mean_rel_diff)
            self.theta = new_theta
            epoch += 1
            print(f"Iteration: {epoch}, Median Relative Difference: {med_rel_diff}")
        return self.theta
    
    def get_results(self):
        tpms = self.theta * 1e6
        tpms[tpms==0.] = 0.
        effective_lengths = (self.transcript_lengths - self.mu + 1) * (self.transcript_lengths >= 200)
        df = pd.DataFrame({
            'transcript': self.transcript_ids,
            'length': self.transcript_lengths,
            'effective_length': effective_lengths,
            'expected_counts': self.expected_counts,
            'tpm': tpms
        })
        df['expected_counts'] = df['expected_counts'].astype(int)
        return df, self.median_relative_differences, self.mean_relative_differences
        

