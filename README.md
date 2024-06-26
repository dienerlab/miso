<img src="https://github.com/dienerlab/miso/raw/main/inst/extdata/logo.png" width="25%" align="right">

![R-CMD-check](https://github.com/dienerlab/miso/workflows/R-CMD-check/badge.svg)
[![codecov](https://codecov.io/gh/dienerlab/miso/branch/main/graph/badge.svg)](https://codecov.io/gh/dienerlab/miso)



`miso`, short for **mi**crobiome **so**ftware is a collection of helpers that we use to analyze microbiome
data. It makes it easier to run some common analyses and is pretty opinionated towards our own experiences.

<br><br>

## General philosophy

`miso` is simply a one stop location for smaller analysis helper and visualizations we
often use in the lab. It is supposed to remove some common pain points we encountered or
to implement custom (mostly genomic) analysis steps.

At this point complex workflows in the lab have been ported to nextflow and are no longer
included here. See our [pipelines](https://github.com/dienerlab/pipelines) for this.

## Analysis

For `misos` an analysis step is based on input data and a configuration,
thus having the function signature `step(object, config)`.
Most steps can be chained with the pipe operator.
For instance, the following is possible with `miso`:

```r
library(miso)

config <- list(
    demultiplex = config_demultiplex(barcodes = c("ACGTA", "AGCTT")),
    preprocess = config_preprocess(truncLen = 200),
    denoise = config_denoise()
)

output <- find_read_files("raw") %>%
          demultiplex(config$demultiplex) %>%
          quality_control() %>%
          preprocess(config$preprocess) %>%
          denoise(config$denoise)
```

This clearly logs the used workflow and the configuration. The configuration
can also be saved and read in many formats, for instance yaml.

config.yaml:
```
preprocess:
  threads: yes
  out_dir: preprocessed
  trimLeft: 10.0
  truncLen: 200.0
  maxEE: 2.0
denoise:
  threads: yes
  nbases: 2.5e+08
  pool: no
  bootstrap_confidence: 0.5
  taxa_db: https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz?download=1
  species_db: https://zenodo.org/record/1172783/files/silva_species_assignment_v132.fa.gz?download=1
  hash: yes
```

This can now be reused by someone else:

```r
config <- read_yaml("config.yml")

output <- find_read_files("raw") %>%
          quality_control() %>%
          preprocess(config$preprocess) %>%
          denoise(config$denoise)
```

## Other functions

All other functions are usually functions that are meant to be inside
more complex code or functions that produce plots and endpoints of
an analysis. Most of them act on phyloseq objects and some on tidy data
tables. Some are general lab helpers, for instance to make plate layouts.
