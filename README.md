# Welcome to **`SpatialE`**!

[![R >= 2.10](https://img.shields.io/badge/R-%3E%3D%202.10-brightgreen)](https://www.r-project.org/)

<p align="center">
  <h3 align="center">SpatialE identifies spatial enrichment of amyotrophic lateral sclerosis genes in spatial transcriptomics</h3>
</p>

**SpatialE** is a gene set spatial enrichment tool that estimates the entropy weighted differential gene expression matrix at spatial clusters, and calculates the significance of how a predefined gene set is enriched to a spatial cluster (region) when compared to randomly sampling gene sets. We packaged it into a user-friendly generic R package. Combined with [Seurat](https://satijalab.org/seurat/index.html), we can determine the spatial region where a gene set or cell type is enriched.

<img src="https://github.com/TJ-zhanglab/SpatialE/blob/main/Web Image/SpatialE Workflow.png" width="1000"/>

## Installation

```r
install.packages("devtools")
library(devtools)

options(timeout=9999999) # Set a lengthy timeout for the SpatialE download.
install_github("https://github.com/TJ-zhanglab/SpatialE.git")
```

## Repository structure

**`Example`**: A demonstration of how to use SpatialE, including **Data demo** and **User's guidance**.

**`R`**: The directory contains scripts(.R) for SpatialE

**`man`**: The directory contains scripts(.Rd) for SpatialE

**`Web Image`**: The directory contains README image files (only for Github website display)

## Usage

SpatialE includes built-in datasets for easy accessibility. After installation, refer to the **`Example`** for usage instructions.

## Issues

SpatialE is still undergoing development and enhancement. We welcome all feedback, error reports, and improvement suggestions. Please feel free to post a comment on our Github.

## References

Will add paper link after publishing.

## Copyright

This tool is developed in TJ-zhanglab.

All rights reserved.


