
# probepicker

A CLI for identifying the most variable probes in larger-than-memory methylation array datasets.

## Overview

The `probepicker` tool allows you to process large methylation datasets and identify the most variable probes. This can be useful in various genomic studies where variability in methylation is of interest.

## Installation

Download the precompiled binaries in the latest release

## Usage

probepicker expects the input data to be a gzipped TSV file where the first column contains probe names, and all other columns contain sample beta values.

```
# Identify the most variable probes (file-input)
./probepicker identify example.tsv.gz --delim "\t" --nprobes 10

# Identify the most variable probes (stdin)
gzip -dc example.tsv.gz | ./probepicker identify --delim "\t" --nprobes 10


# Select a subset of samples and identify the most variable probes on this subset
gzip -dc example.tsv.gz | ./probepicker select --samples samples.txt | ./probepicker identify > probes.txt
```


An example TSV file with methylation data is below
```
probes    sample1           sample2    sample3             sample4             sample5
cg00000029 0.098069738       0.909142796 0.8603162695913249 0.28886065277317896 0.089729344
cg00000165 0.8837011024652779 0.870551483 0.28391987868948804 0.831906752       0.16676024603610298
cg00000236 0.9206346711071329 0.9170243029210999 0.92350756 0.913890037         0.920367436
cg00000289 0.7360994844654579 0.886445714 0.8526066164315022 0.732923475        0.684158594
cg00000292 0.8537021401033309 0.885411564 0.9211661560807459 0.6791142545452471 0.8747384317623911
cg00000321 0.784947493        0.854457667 0.31709493256129806 0.7733845509958441 0.45572892772571005
cg00000363 0.5085552760021951 0.8419091764229979 0.796398249 0.179835181        0.7268788134611089

```

An example sample list is below

```
TCGA-OR-A5JP-01A
TCGA-OR-A5JG-01A
TCGA-PK-A5HB-01A
TCGA-OR-A5JE-01A
TCGA-OR-A5KU-01A
TCGA-OR-A5L9-01A
TCGA-OR-A5JQ-01A
TCGA-OR-A5K4-01A
TCGA-OR-A5JL-01A
TCGA-OR-A5LS-01A
```