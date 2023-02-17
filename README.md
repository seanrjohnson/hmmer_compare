# hmmer_compare
compare hmmer3 profiles and generate profile trees.

It's kind of slow, so not suitable for comparing large numbers of profiles. More than 100 is ambitious.

Note: This package is a small slice of a more comprehensive tools suite that I hope to make available before the end of 2023. So apart from bugfixes, I may not maintain it much (unless I decide to use it as a dependency for the broader suite).


# INSTALLATION

## Download and extract the hmmer_compare package

Go to the [hmmer_compare github page](https://github.com/seanrjohnson/hmmer_compare) and download and extract the zip file for the latest release (see the right side of the github page).

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

# Example

```

```

# Contributing

I think this program could be thousands of times faster if the core algorithm was rewritten in C/cython (or maybe Rust). I haven't had time to do that, but I would be very grateful to anyone who does have time!

