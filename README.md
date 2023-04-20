# RNASeq-Quantification

[IsoEM](https://almob.biomedcentral.com/articles/10.1186/1748-7188-6-9) expectation maximization implementation for quantifying relative transcript expression using RNA-Seq data.

Achieved comparable median performance to [Kallisto](https://www.nature.com/articles/nbt.3519).

### References

1. Nicolae, M., Mangul, S., Măndoiu, I.I. et al. Estimation of alternative splicing isoform frequencies from RNA-Seq data. Algorithms Mol Biol 6, 9 (2011). https://doi.org/10.1186/1748-7188-6-9
2. Bray, N., Pimentel, H., Melsted, P. et al. Near-optimal probabilistic RNA-seq quantification. Nat Biotechnol 34, 525–527 (2016). https://doi.org/10.1038/nbt.3519

### Implementation Note

Data is in a separate folder. It contains a BAM file with all alignments and a FASTA file containing the transcriptome.
