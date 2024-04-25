# Welcome to **`SpatialE`**!

<p align="center">
  <h3 align="center">Spatial enrichment and genomic analyses reveal the link of NOMO1 with amyotrophic lateral sclerosis</h3>
</p>

**SpatialE** is a gene set spatial enrichment tool that estimates the entropy weighted differential gene expression matrix at spatial clusters, and calculates the significance of how a predefined gene set is enriched to a spatial cluster (region) when compared to randomly sampling gene sets. We packaged it into a user-friendly generic R package (installation and usage tutorials are on this GitHub), and also provided a convenient-to-use website (https://www.aibraingenomics.com/spatiale/). By Combing [Seurat](https://satijalab.org/seurat/index.html), we can estimate the spatial region where a gene set or cell type is enriched to.

<img src="https://github.com/TJ-zhanglab/SpatialE/blob/master/Web%20Image/SpatialE%20workflow.png" width="1000"/>

## Installation

```r
install.packages("devtools")
library(devtools) ## Install dependencies

options(timeout=9999999) ## Set a lengthy timeout for the SpatialE download.
install_github("https://github.com/TJ-zhanglab/SpatialE.git")
```

## Repository structure

**`R`**: The directory contains scripts(.R) for SpatialE

**`User's guidance`**: A demonstration of how to use SpatialE

**`Web Image`**: The directory contains image files (only for Github website display)

**`data`**: Data used in the **`User's guidance`**

**`man`**: The directory contains scripts(.Rd) for SpatialE

## Usage

SpatialE includes built-in datasets for easy accessibility. After installation, refer to the **`User's guidance`** for usage instructions.

## Issues

SpatialE is still undergoing development and enhancement. We welcome all feedback, error reports, and improvement suggestions. Please feel free to post a comment on our Github.

## References

Will add paper link after publishing.

## Copyright

This tool is developed in TJ-zhanglab.

All rights reserved.


