# hmmer_compare
compare hmmer3 profiles, generate profile alignments, and generate profile trees.

It's kind of slow, so not suitable for comparing large numbers of profiles. More than 100 is ambitious.

Note: This package is a small slice of a more comprehensive tools suite that I hope to make available before the end of 2023. So apart from bugfixes, I may not maintain it much (unless I decide to use it as a dependency for the broader suite).


# INSTALLATION

## Download and extract the hmmer_compare package

Clone the repo or download a zip file of the repo.

To download a zip of the repo, go to the [hmmer_compare github page](https://github.com/seanrjohnson/hmmer_compare) and download and extract the zip file for the latest release (see the right side of the github page).

Extract the zip file.  In a terminal, navigate to the newly extracted folder and proceed with either the pip or Anaconda install instructions


## pip
On a system where you have python3 installed

```
    pip install -e .
```


## Anaconda
### Install Anaconda (if you don't already have it)

[https://docs.conda.io/en/main/miniconda.html](https://docs.conda.io/en/main/miniconda.html)

or

[https://docs.anaconda.com/anaconda/install/](https://docs.anaconda.com/anaconda/install/)



### install hmmer_compare
```
conda env create -f conda_env.yml
```

This will create a new conda environment called "hmmer_compare".

### Activate the environment
```
conda activate hmmer_compare
```

# Examples

starting in the same directory as this readme file
## Creating a UPGMA tree of hmmer3 profiles

```
hmmer_compare.py -i test/data/pdonr_hmms.hmm -r test/data/pdonr_hmms.hmm -o scores.tsv
table_to_tree.py -i scores.tsv -o pdonr_hmms.newick
```

## Creating hmmer3 profile alignments
```
hmmer_compare.py -i test/data/pdonr_hmms.hmm -r test/data/pdonr_hmms.hmm -o scores_with_alignments.txt --alignments
```
Note that `scores_with_alignments.txt` can also be used as input to `table_to_tree.py` as it will only read lines with three tab-separated columns, and thereby skip the alignments.


# Contributing

I think this program could be thousands of times faster if the core algorithm was rewritten in C/cython (or maybe Rust). I haven't had time to do that, but I would be very grateful to anyone who does have time! Any other suggestions/contributions are also welcome.

# Dependencies
    - numpy
    - scipy
    - pyhmmer

# References

Adapted from pseudocode in:

- Steinegger, Martin, Markus Meier, Milot Mirdita, Harald Vöhringer, Stephan J. Haunsberger, and Johannes Söding. “HH-Suite3 for Fast Remote Homology Detection and Deep Protein Annotation.” BMC Bioinformatics 20, no. 1 (September 14, 2019): 473. [https://doi.org/10.1186/s12859-019-3019-7](https://doi.org/10.1186/s12859-019-3019-7).
- Söding, Johannes. “Protein Homology Detection by HMM–HMM Comparison.” Bioinformatics 21, no. 7 (April 1, 2005): 951–60. [https://doi.org/10.1093/bioinformatics/bti125](https://doi.org/10.1093/bioinformatics/bti125).
