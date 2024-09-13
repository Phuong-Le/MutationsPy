# MutationsPy

A Python package to manipulate mutation data. Package is under development 

Features are roughly similar to [MutationsR](https://github.com/Phuong-Le/mutationsR) 

# Installations

```
pip install git+https://github.com/Phuong-Le/MutationsPy
```

# Usage highlights

The package supports command line interface to convert mutation matrix from hdp to sigprofiler/msighdp format and vice versa

```
# converting from sigprofiler to hdp (trinucleotide SNVs only)
mutation_matrix sigprofiler_to_hdp --sigprofiler_path <path/to/sigprofiler> --hdp_outpath <output/path/to/hdp>

# converting from hdp to sigprofiler (SNVs only)
mutation_matrix hdp_to_sigprofiler --hdp_path <path/to/hdp> --sigprofiler_outpath <output/path/to/sigprofiler>

# concatenate multiple mutation matrices
mutation_matrix concat_sigprofiler_mutmats --mutmats <path/to/matrix1> <path/to/matrix2> ... <path/to/matrixn> --outpath <output/path>
```