# Some ideas about ChIP-Seq

> Chromatin-immunoprecipitation sequencing

## Supervised Peak Calling for ChIP-Seq

* Many peak detection algorithms have been proposed for CHIP-seq data analysis, **but it is not obvious which algorithm and what parameters are optimal for any given data set.**
* Conducted by Toby Hocking, their group propose a supervised approach that uses labels produced from visual inspection of the aligned read counts.
* The main idea is to **manually label** a small subset of the genome, and then learn a model that makes consistent peak predictions on the rest of the genome.

## Peak Calling for ChIP-seq

* Read Shifting: the aligned  reads are form fragments of 150-300bp in length and, as most ChIP-seq data from single-end sequencing, only one end of a fragment is read.

## Reference

* [Analysis of ChIP-seq data](https://galaxyproject.org/tutorials/chip/)
* [Peak Calling for ChIP-seq](https://www.sigmaaldrich.com/technical-documents/articles/biology/peak-calling-for-chip-seq.html)
* [Supervised Peak Calling for ChIP-Seq](http://www.computationalgenomics.ca/BourqueLab/project/supervised-peak-calling-for-chip-seq/)