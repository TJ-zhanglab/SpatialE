# Welcome to **`SpatialE`**!

Here, we demonstrate how to use SpatialE to analyze the enrichment of gene sets in Spatial Transcriptome (ST) data. The analysis process mainly includes the following parts.

- Loading the required R packages 
- Preparing the ST Seurat object 
- Running SpatialE
   - Processing the ST expression matrix by clusters
   - Filtering genes with very low expression
   - Calculating the delta matrix
   - Weighting each gene based on the Entropy Weight Method
   - Importing gene set(s) of interest and conducting spatial enrichment analysis
- Visualizing spatial enrichment results

## Step1: Loading the required R packages

Before starting the analysis, we will load SpatialE and the other necessary packages for this vignette.

[![Seurat-4.3.0](https://img.shields.io/badge/Seurat-4.3.0-green)](https://cran.r-project.org/web/packages/Seurat)
[![dplyr-1.0.10](https://img.shields.io/badge/dplyr-1.0.10-orange)](https://cran.r-project.org/web/packages/dplyr)

```r
library(Seurat)

library(dplyr)

library(SpatialE) 
```
## Step2: Preparing the ST Seurat object 

We will use the ST data of a posterior coronal slice of mouse brain that was generated with the Visium technology from [10x Genomics](https://www.10xgenomics.com/resources/datasets?query=&page=1&configure%5Bfacets%5D%5B0%5D=chemistryVersionAndThroughput&configure%5Bfacets%5D%5B1%5D=pipeline.version&configure%5BhitsPerPage%5D=500&menu%5Bproducts.name%5D=Spatial%20Gene%20Expression#:~:text=Search-,Datasets,-Products), which contains the expression of 32,285 genes at 2,702 spatial spots. We packaged it in the SpatialE, which can be imported directly:

```r
data(Mouse_Brain_ST_Demo)
```

It will return a Seurat object containing both the spot-level expression matrix and the associated image of the tissue slice. We referred the data preprocessing workflow of [Seurat](https://satijalab.org/seurat/index.html) to perform dimension reduction and clustering on the ST data, and saved the clustering results in the built-in dataset of SpatialE. We can use the `SpatialFeaturePlot()` function in [Seurat](https://satijalab.org/seurat/index.html) to visualize the preprocessed results.

```r
SpatialDimPlot(Mouse_Brain_ST_Demo, label = T, label.size = 7)
```

<img src="https://github.com/TJ-zhanglab/SpatialE/blob/master/Web%20Image/SpatialDimPlot.png" width="600"/>

Different colors represent different clusters (idents) and correspond to different histological regions. We suggest that before the next step analysis, the parameters of data preprocessing, especially the spatial dimension reduction parameters, should be adjusted appropriately according to the known anatomical positions. This process is to ensure that the spatial clusters (idents) can well distinguish the histological positions.

## Step3: Running SpatialE
SpatialE is a tool to analyze the enrichment of a target/predefined gene set (can be a set of cell type marker genes, representing specific cell type) to spatial regions of a tissue slice. It uses the **Entropy Weight Method** to assign different **weights** to different marker genes, which could be used to generate a specific spatial enrichment of specific cell types.

### Processing the ST expression matrix by clusters
Firstly, we will use the `getExpMatrix()` function in SpatialE to get the preprocessed ST expression matrix from the Seurat object, where each row represents a gene, and each column represents a spatial spot.

```r
exp <- getExpMatrix(Mouse_Brain_ST_Demo)
str(exp)
```

```
##  num [1:18768, 1:2702] 0 0 0 1 1 3 1 9 0 0 ...
##  - attr(*, "dimnames")=List of 2
##   ..$ : chr [1:18768] "Xkr4" "Gm19938" "Sox17" "Mrpl15" ...
##   ..$ : chr [1:2702] "AAACAAGTATCTCCCA-1" "AAACAATCTACTAGCA-1" "AAACACCAATAACTGC-1" "AAACAGAGCGACTCCT-1" ...
```

From the result of `str(exp)`, we can know that preprocessed Seurat object contains 18,768 genes and 2,702 spots. Secondly, we will divide the expression matrix by spatial clusters and calculate the average expression of each gene in different spatial clusters.


```r
cluster_spot_exp <- getClusteredExp(Mouse_Brain_ST_Demo, exp)
gene_clustered_mean <- Clustered_mean_Exp(cluster_spot_exp)
```

Below are the top 5 rows of `gene_clustered_mean` matrix, where each row represents a gene, each column represents a spatial cluster, and each matrix element represents the mean expression.

<table class="table table-responsive-{sm|md|lg|xl}" style="font-size: 1px; width: auto !important; ">

 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> cluster0 </th>
   <th style="text-align:right;"> cluster1 </th>
   <th style="text-align:right;"> cluster2 </th>
   <th style="text-align:right;"> cluster3 </th>
   <th style="text-align:right;"> cluster4 </th>
   <th style="text-align:right;"> cluster5 </th>
   <th style="text-align:right;"> cluster6 </th>
   <th style="text-align:right;"> cluster7 </th>
   <th style="text-align:right;"> cluster8 </th>
   <th style="text-align:right;"> cluster9 </th>
   <th style="text-align:right;"> cluster10 </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Xkr4 </td>
   <td style="text-align:right;"> 0.05319149 </td>
   <td style="text-align:right;"> 0.1066667 </td>
   <td style="text-align:right;"> 0.0921659 </td>
   <td style="text-align:right;"> 0.07978723 </td>
   <td style="text-align:right;"> 0.09917355 </td>
   <td style="text-align:right;"> 0.07586207 </td>
   <td style="text-align:right;"> 0.144796380 </td>
   <td style="text-align:right;"> 0.07361963 </td>
   <td style="text-align:right;"> 0.12941176 </td>
   <td style="text-align:right;"> 0.06666667 </td>
   <td style="text-align:right;"> 0.13026820 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Gm19938 </td>
   <td style="text-align:right;"> 0.09574468 </td>
   <td style="text-align:right;"> 0.1650000 </td>
   <td style="text-align:right;"> 0.0875576 </td>
   <td style="text-align:right;"> 0.10638298 </td>
   <td style="text-align:right;"> 0.13223140 </td>
   <td style="text-align:right;"> 0.09310345 </td>
   <td style="text-align:right;"> 0.144796380 </td>
   <td style="text-align:right;"> 0.09202454 </td>
   <td style="text-align:right;"> 0.14117647 </td>
   <td style="text-align:right;"> 0.09629630 </td>
   <td style="text-align:right;"> 0.14176245 </td> 
  </tr>
  <tr>
   <td style="text-align:left;"> Sox17 </td>
   <td style="text-align:right;"> 0.13829787 </td>
   <td style="text-align:right;"> 0.1650000 </td>
   <td style="text-align:right;"> 0.1474654 </td>
   <td style="text-align:right;"> 0.09574468 </td>
   <td style="text-align:right;"> 0.13774105 </td>
   <td style="text-align:right;"> 0.10689655 </td>
   <td style="text-align:right;"> 0.167420814 </td>
   <td style="text-align:right;"> 0.09202454 </td>
   <td style="text-align:right;"> 0.12941176 </td>
   <td style="text-align:right;"> 0.15555556 </td>
   <td style="text-align:right;"> 0.14942529 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Mrpl15 </td>
   <td style="text-align:right;"> 1.06382979 </td>
   <td style="text-align:right;"> 1.2150000 </td>
   <td style="text-align:right;"> 1.2949309 </td>
   <td style="text-align:right;"> 1.44148936 </td>
   <td style="text-align:right;"> 0.93388430 </td>
   <td style="text-align:right;"> 1.23448276 </td>
   <td style="text-align:right;"> 1.190045249 </td>
   <td style="text-align:right;"> 1.28834356 </td>
   <td style="text-align:right;"> 1.34705882 </td>
   <td style="text-align:right;"> 1.22222222 </td>
   <td style="text-align:right;"> 1.32950192 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Lypla1 </td>
   <td style="text-align:right;"> 0.59574468 </td>
   <td style="text-align:right;"> 0.8350000 </td>
   <td style="text-align:right;"> 0.5391705 </td>
   <td style="text-align:right;"> 0.91489362 </td>
   <td style="text-align:right;"> 0.60881543 </td>
   <td style="text-align:right;"> 0.59310345 </td>
   <td style="text-align:right;"> 0.533936652 </td>
   <td style="text-align:right;"> 0.71779141 </td>
   <td style="text-align:right;"> 0.85294118 </td>
   <td style="text-align:right;"> 0.83703704 </td>
   <td style="text-align:right;"> 0.60919540 </td>
  </tr>
</tbody>
</table>

### Filtering genes with very low expression

In consideration of the actual sequencing quality, we will consider filtering genes with extremely low expression in order to avoid extreme value bias in the results. Here, we will calculate the row maximum of the `gene_clustered_mean` matrix, which is the maximum of the average expression of each gene across clusters. And set the filter threshold by analyzing the distribution of these maximum values.

We recommend setting the filter threshold to a value near the lower quartile (Q1) of maximum values. More importantly, we hope that users could customize the filter threshold value according to the experiment and biological characteristics of the given dataset. The larger the threshold is, the more genes will be filtered out.

```r
gene_clustered_max <- mutate(gene_clustered_mean, max = apply(gene_clustered_mean,1,max)) # Calculate row maximum and generate a new column for storage
summary(gene_clustered_max$max) # View the distribution of maximum values
gene_clustered_mean_filtered <- genefilter(gene_clustered_mean, threshhold = 0.1) # Set filter threshold

nrow(gene_clustered_max) # Number of genes before filtering
```

```
## [1] 18768
```

```r
nrow(gene_clustered_mean_filtered) # Number of genes after filtering
```

```
## [1] 13471
```

There are 18,768 genes in total, 5,279 genes will be filtered and 13,471 genes will be left.<br>
<br>
<details>
  <summary> Why is the value near the lower quartile set as the filter threshold? </summary>

The quartile is the value at the 25% (lower quartile or Q1) or 75% (upper quartile or Q3) position after sorting data from small to large. By comparing Q1 and Q3, we can basically know the distribution of data variables, such as gene expression levels in different spatial clusters. Similar to the principle of **Three Sigma Guidelines** in Normal distribution, the quartile is suitable for any dataset. Whether or not its distribution is known. In most cases, the Normal distribution cannot be fitted precisely, and the quantile could be a good choice.<br>
</details>
<br>

<details>
  <summary> Why not filter the part above the upper quartile? </summary>

The upper and lower quartiles represent genes with very extreme expression. The lower quartile corresponds to genes with very low expression, or even almost no expression. The upper quartile corresponds to genes with very high expression, or even only expressed in specific cell types, which are consistent with the definition of most marker genes. Therefore, we recommend considering the lower filter limit.<br>
</details>
<br>

<details>
  <summary> Why not directly consider the distribution of the entire average expression matrix? </summary>

We have also considered alternative threshold-setting methods. The first is to consider the data distribution of the entire average expression matrix `gene_clustered_mean` without distinguishing between genes and spatial clusters. The second is to calculate the column maximum values (focus on the distribution of gene expression in each cluster). Both of these cases are more likely to filter spatial clusters than genes, which cannot achieve our goal, and we cannot conduct enrichment analysis on incomplete spatial clusters.<br>
</details>

### Calculating the delta matrix

In this step, we will calculate the differential gene expression matrix. The **difference** here refers to the difference between the average expression of the gene set in a spatial cluster and the average expression of the gene set in all spatial spots. Each row of the differential matrix represents a gene, and each column represents a spatial cluster. The matrix elements change from the expression to the the expression difference.<br>

We design three modes (**“high expression”**, **“two sides”** and **“low expression”**) to calculate the differential gene expression matrix, which corresponds to a high-expression related differential matrix (high expression) that only considers the deviations derived from highly expressed genes; a general differential matrix (two sides) that considers the deviations from both directions; and a low-expression related differential matrix (low expression) that only considers the deviations derived from lowly expressed genes. Biologically, the highly expressed genes are often important functional genes (eg. cell markers).<br>

</font></center>
<br>
In the "high expression" mode, the average expression of a gene in a spatial cluster is generally higher than the average expression of all spatial spots. The expression difference is specifically defined as $E_{g,c}-E_g$ , the $E_{g,c}$ represents the average expression value of the gene $g$ in the spatial cluster $c$, the $E_g$ represents the average expression value of the gene $g$ in all spots. The gene expression difference can be defined as $D_{g,c}$. For the part lower than the average expression of all spatial spots, in order to avoid negative value, we will define it as zero.<br>
<br>

$$D_{g,c}=
\begin{cases}
E_{g,c} - E_g, & (E_{g,c} > E_g)\\
0, & (E_{g,c} \leq E_g)
\end{cases}
$$

</font></center>
<br>
In the "two sides" mode, the average expression of a gene in a spatial cluster may be higher or lower than the average expression of all spatial spots. We calculate the absolute value of the difference to ensure that the result is not negative. The meaning of variables in the formula is the same as the "high expression" mode.<br>
<br>

$$D_{g,c} = \lvert E_{g,c} - E_g\vert$$

</font></center>
<br>
In the "low expression" mode, the average expression of a gene in a spatial cluster is generally lower than the average expression of all spatial spots. For the part higher than the average expression of all spatial clusters, in order to avoid negative value, we will define it as zero.The meaning of variables in the formula is the same as the "high expression" mode.<br>
<br>

$$D_{g,c}=
\begin{cases}
E_g - E_{g,c}, & (E_{g,c} < E_g) \\
0, & (E_{g,c} \geq E_g)
\end{cases}
$$

<br>
In this demo, we used the "high expression" mode, as highly expressed genes in a specific cell type or tissue region are typically more important from a biological standpoint (eg.cell marker genes).<br>

```r
# Define the type as "high expression","low expression" or "two sides"
delta <- getdelta(gene_clustered_mean_filtered, exp, type = 'high expression')
```

Through testing in different datasets, we found that the "high expression" mode is more recommended in most cases in order to have a correct enrichment of known cell types to expected tissue regions.<br>

### Weighting each gene based on the Entropy Weight Method

Next, we will use the Entropy Weight method to assign different weights to different genes by their variance.

```r
rescale_gene <- rescale(delta) 
percentage_gene <- percentage(rescale_gene) 
entropy_gene <- entropy(percentage_gene) 
weight_gene <- weight(entropy_gene) 
#Add the corresponding weight into the delta matrix
delta$gene <- rownames(delta)
delta <- left_join(delta, weight_gene, by = 'gene')
rownames(delta) <- delta$gene
```

You can type `?function name()` to learn more about functions. Check the column information of `delta` matrix, and you will find the corresponding weight of gene is added to the last column.

```
 [1] "cluster1"  "cluster2"  "cluster3"  "cluster4"  "cluster5"  "cluster6"  "cluster7"  "cluster8" 
 [9] "cluster9"  "cluster10" "cluster11" "gene"      "weight" 
```
### Importing gene set(s) of interest and conducting spatial enrichment analysis

In this step, we will perform spatial enrichment analysis on the gene set(s) of interest, which may be marker gene sets of various cell types, or any gene set of interest (for example, disease related gene sets). Please provide the gene set(s) in **vector** format.

Here, we will demonstrate the process using the SpatiaE's built-in data sets. As demonstrated below, you can directly load them using `data(cell_type_marker_genes)` and `data(ALS_genes)` functions.

<font size=3.5>**1 cell_type_marker_genes **</font>. This data object is a dataframe that contains marker genes of two cell types: *Layer6 Intra Telencephalic Isocortex(L6 IT CTX)* and *Oligodendrocyte*, the corresponding p values and adjusted p values. It is a single cell transcriptome data of mouse brain tissue published on the [Allen Brain Map website](http://portal.brain-map.org/atlases-and-data/rnaseq/mouse-whole-cortex-and-hippocampus-smart-seq). We obtainded these marker genes by using `FindAllMarkers()` function in [Seurat](https://satijalab.org/seurat/index.html).<br> 
<br>
<font size=3.5>**2 ALS_genes **</font>. This data object is a vector that contains 260 genes related to Amyotrophic lateral sclerosis. Please refer to Table. S7 of paper for detailed data sources.<br>

Loading the gene set.
```r
data(cell_type_marker_genes) #or data(ALS_genes)
```
The data structure of the **ell_type_marker_genes** is shown as below (only the top 5 rows).


<table class="table table-responsive-{sm|md|lg|xl" style="font-size: 5px; width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;text-align: center;"> p_val </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;text-align: center;"> p_val_adj </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;text-align: center;"> cluster </th>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;text-align: center;"> gene </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;text-align: center;"> 0 </td>
   <td style="text-align:right;text-align: center;"> 0 </td>
   <td style="text-align:left;text-align: center;"> L6 IT CTX </td>
   <td style="text-align:left;text-align: center;"> Dmkn </td>
  </tr>
  <tr>
   <td style="text-align:right;text-align: center;"> 0 </td>
   <td style="text-align:right;text-align: center;"> 0 </td>
   <td style="text-align:left;text-align: center;"> L6 IT CTX </td>
   <td style="text-align:left;text-align: center;"> Osr1 </td>
  </tr>
  <tr>
   <td style="text-align:right;text-align: center;"> 0 </td>
   <td style="text-align:right;text-align: center;"> 0 </td>
   <td style="text-align:left;text-align: center;"> L6 IT CTX </td>
   <td style="text-align:left;text-align: center;"> Blnk </td>
  </tr>
  <tr>
   <td style="text-align:right;text-align: center;"> 0 </td>
   <td style="text-align:right;text-align: center;"> 0 </td>
   <td style="text-align:left;text-align: center;"> L6 IT CTX </td>
   <td style="text-align:left;text-align: center;"> Bmp3 </td>
  </tr>
  <tr>
   <td style="text-align:right;text-align: center;"> 0 </td>
   <td style="text-align:right;text-align: center;"> 0 </td>
   <td style="text-align:left;text-align: center;"> L6 IT CTX </td>
   <td style="text-align:left;text-align: center;"> Col6a1 </td>
  </tr>
</tbody>
</table>  

Here we will take **L6 IT CTX** cell type in the **manual_76k** dataset as an example to demonstrate how to perform the spatial enrichment analysis.

```r
target_geneset <- manual_76k[manual_76k$cluster=='L6 IT CTX',] 
```

According to p_val_adj sort values, the top 200 genes will be selected for the next analysis (see below for details on why the top 200 genes were chosen).

```r
target_geneset <- target_geneset[order(target_geneset$p_val_adj),] 
target_geneset <- target_geneset[target_geneset$p_val_adj<0.05,] 
nrow(target_geneset) #Output the number of significant genes
```

```
## [1] 1907
```

```r
target_geneset <- target_geneset[1:200,] #Choose the top 200 marker genes as the target genes
target_geneset
```

```
## # A tibble: 200 × 4
##    p_val p_val_adj cluster   gene    
##    <dbl>     <dbl> <chr>     <chr>   
##  1     0         0 L6 IT CTX Dmkn    
##  2     0         0 L6 IT CTX Osr1    
##  3     0         0 L6 IT CTX Blnk    
##  4     0         0 L6 IT CTX Bmp3    
##  5     0         0 L6 IT CTX Col6a1  
##  6     0         0 L6 IT CTX Slc26a4 
##  7     0         0 L6 IT CTX Sulf1   
##  8     0         0 L6 IT CTX Adamts13
##  9     0         0 L6 IT CTX Aspg    
## 10     0         0 L6 IT CTX Dnah14  
## # … with 190 more rows
```
<details>
  <summary>Why do we recommend the top 200 marker genes for analysis?</summary>

Although the p values of significant genes are less than 0.05, if we use all significant genes for enrichment analysis, it will not only take a long time to calculate, occupy a large GPU memory, but also cause high false positives. However, few significant genes may cause information loss. After internal testings, we recommend selecting the top 100 genes (a total of fewer than 1,000 significant genes), the top 200 genes (a total of 1,000 to 5,000 significant genes), and the top 500 genes (a total of more than 5,000 significant genes) as the target genes. There are 1,907 significant genes in L6 IT CTX and we retain the top 200 genes. If there are more significant genes in your datasets, please consider selecting 800, 1000, or more genes.<br>
</details>
<br>

<details>
  <summary>Why not consider setting the threshold value of p value to determine the top number？</summary>

The reason for not using a fixed p value to determine the top number is that the p values of different cell types vary greatly, particularly in cell types with sparse p value distributions. We should pay more attention to the relative difference between p-values for a particular cell type. In the future, we will use other algorithms to solve this issue.<br>
</details>
<br>

**How to judge which cluster(s) is the target gene set significantly enriched to?**<br>

Here we define a judgment baseline, randomly sample 10,000 gene sets with the same length as the target gene set, construct their expression matrices, compare the sum of gene weighted expression between the target gene set and 10,000 random gene sets, and calculate the probability that the weighted sum of the target gene set is greater than the random gene sets.<br>

```r
target_delta <- delta[delta$gene %in% target_geneset$gene,] 
target_gama <- getGama(target_delta, type = 'target') 
random_sample_delta <- getRandomSample(delta, target_geneset$gene) 
random_sample_gama <- getGama(random_sample_delta, type = 'random', cores = 30)
enrichment_output <- getEnrichmentOutput(target_gama, random_sample_gama) 
colnames(enrichment_output)<- c("frequency","possibility","pval","cluster","fdr","bonferroni","graphdata")
```
You can type `?function name()` to learn more about functions. The `enrichment_output` dataframe is shown as below.

<table class="table table-responsive-{sm|md|lg|xl" style="font-size: 5px; width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;position: sticky; top:0; background-color: #FFFFFF;text-align: center;">   </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;text-align: center;"> frequency </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;text-align: center;"> possibility </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;text-align: center;"> pval </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;text-align: center;"> cluster </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;text-align: center;"> fdr </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;text-align: center;"> bonferroni </th>
   <th style="text-align:right;position: sticky; top:0; background-color: #FFFFFF;text-align: center;"> graphdata </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;font-weight: bold;color: red !important;text-align: center;"> 9 </td>
   <td style="text-align:right;font-weight: bold;color: red !important;text-align: center;"> 9999 </td>
   <td style="text-align:right;font-weight: bold;color: red !important;text-align: center;"> 0.9999 </td>
   <td style="text-align:right;font-weight: bold;color: red !important;text-align: center;"> 0.0001 </td>
   <td style="text-align:right;font-weight: bold;color: red !important;text-align: center;"> 8 </td>
   <td style="text-align:right;font-weight: bold;color: red !important;text-align: center;"> 0.0019000 </td>
   <td style="text-align:right;font-weight: bold;color: red !important;text-align: center;"> 0.0019 </td>
   <td style="text-align:right;font-weight: bold;color: red !important;text-align: center;"> 2.7212464 </td>
  </tr>
  <tr>
   <td style="text-align:left;text-align: center;"> 17 </td>
   <td style="text-align:right;text-align: center;"> 9908 </td>
   <td style="text-align:right;text-align: center;"> 0.9908 </td>
   <td style="text-align:right;text-align: center;"> 0.0092 </td>
   <td style="text-align:right;text-align: center;"> 16 </td>
   <td style="text-align:right;text-align: center;"> 0.0874000 </td>
   <td style="text-align:right;text-align: center;"> 0.1748 </td>
   <td style="text-align:right;text-align: center;"> 0.7574586 </td>
  </tr>
  <tr>
   <td style="text-align:left;text-align: center;"> 13 </td>
   <td style="text-align:right;text-align: center;"> 9699 </td>
   <td style="text-align:right;text-align: center;"> 0.9699 </td>
   <td style="text-align:right;text-align: center;"> 0.0301 </td>
   <td style="text-align:right;text-align: center;"> 12 </td>
   <td style="text-align:right;text-align: center;"> 0.1562750 </td>
   <td style="text-align:right;text-align: center;"> 0.5719 </td>
   <td style="text-align:right;text-align: center;"> 0.2426799 </td>
  </tr>
  <tr>
   <td style="text-align:left;text-align: center;"> 18 </td>
   <td style="text-align:right;text-align: center;"> 9671 </td>
   <td style="text-align:right;text-align: center;"> 0.9671 </td>
   <td style="text-align:right;text-align: center;"> 0.0329 </td>
   <td style="text-align:right;text-align: center;"> 17 </td>
   <td style="text-align:right;text-align: center;"> 0.1562750 </td>
   <td style="text-align:right;text-align: center;"> 0.6251 </td>
   <td style="text-align:right;text-align: center;"> 0.2040505 </td>
  </tr>
  <tr>
   <td style="text-align:left;text-align: center;"> 8 </td>
   <td style="text-align:right;text-align: center;"> 9512 </td>
   <td style="text-align:right;text-align: center;"> 0.9512 </td>
   <td style="text-align:right;text-align: center;"> 0.0488 </td>
   <td style="text-align:right;text-align: center;"> 7 </td>
   <td style="text-align:right;text-align: center;"> 0.1854400 </td>
   <td style="text-align:right;text-align: center;"> 0.9272 </td>
   <td style="text-align:right;text-align: center;"> 0.0328266 </td>
  </tr>
  <tr>
   <td style="text-align:left;text-align: center;"> 15 </td>
   <td style="text-align:right;text-align: center;"> 9160 </td>
   <td style="text-align:right;text-align: center;"> 0.9160 </td>
   <td style="text-align:right;text-align: center;"> 0.0840 </td>
   <td style="text-align:right;text-align: center;"> 14 </td>
   <td style="text-align:right;text-align: center;"> 0.2460500 </td>
   <td style="text-align:right;text-align: center;"> 1.0000 </td>
   <td style="text-align:right;text-align: center;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;text-align: center;"> 4 </td>
   <td style="text-align:right;text-align: center;"> 8976 </td>
   <td style="text-align:right;text-align: center;"> 0.8976 </td>
   <td style="text-align:right;text-align: center;"> 0.1024 </td>
   <td style="text-align:right;text-align: center;"> 3 </td>
   <td style="text-align:right;text-align: center;"> 0.2460500 </td>
   <td style="text-align:right;text-align: center;"> 1.0000 </td>
   <td style="text-align:right;text-align: center;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;text-align: center;"> 1 </td>
   <td style="text-align:right;text-align: center;"> 8794 </td>
   <td style="text-align:right;text-align: center;"> 0.8794 </td>
   <td style="text-align:right;text-align: center;"> 0.1206 </td>
   <td style="text-align:right;text-align: center;"> 0 </td>
   <td style="text-align:right;text-align: center;"> 0.2460500 </td>
   <td style="text-align:right;text-align: center;"> 1.0000 </td>
   <td style="text-align:right;text-align: center;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;text-align: center;"> 19 </td>
   <td style="text-align:right;text-align: center;"> 8719 </td>
   <td style="text-align:right;text-align: center;"> 0.8719 </td>
   <td style="text-align:right;text-align: center;"> 0.1281 </td>
   <td style="text-align:right;text-align: center;"> 18 </td>
   <td style="text-align:right;text-align: center;"> 0.2460500 </td>
   <td style="text-align:right;text-align: center;"> 1.0000 </td>
   <td style="text-align:right;text-align: center;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;text-align: center;"> 10 </td>
   <td style="text-align:right;text-align: center;"> 8705 </td>
   <td style="text-align:right;text-align: center;"> 0.8705 </td>
   <td style="text-align:right;text-align: center;"> 0.1295 </td>
   <td style="text-align:right;text-align: center;"> 9 </td>
   <td style="text-align:right;text-align: center;"> 0.2460500 </td>
   <td style="text-align:right;text-align: center;"> 1.0000 </td>
   <td style="text-align:right;text-align: center;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;text-align: center;"> 7 </td>
   <td style="text-align:right;text-align: center;"> 8402 </td>
   <td style="text-align:right;text-align: center;"> 0.8402 </td>
   <td style="text-align:right;text-align: center;"> 0.1598 </td>
   <td style="text-align:right;text-align: center;"> 6 </td>
   <td style="text-align:right;text-align: center;"> 0.2760182 </td>
   <td style="text-align:right;text-align: center;"> 1.0000 </td>
   <td style="text-align:right;text-align: center;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;text-align: center;"> 16 </td>
   <td style="text-align:right;text-align: center;"> 7468 </td>
   <td style="text-align:right;text-align: center;"> 0.7468 </td>
   <td style="text-align:right;text-align: center;"> 0.2532 </td>
   <td style="text-align:right;text-align: center;"> 15 </td>
   <td style="text-align:right;text-align: center;"> 0.4009000 </td>
   <td style="text-align:right;text-align: center;"> 1.0000 </td>
   <td style="text-align:right;text-align: center;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;text-align: center;"> 11 </td>
   <td style="text-align:right;text-align: center;"> 7126 </td>
   <td style="text-align:right;text-align: center;"> 0.7126 </td>
   <td style="text-align:right;text-align: center;"> 0.2874 </td>
   <td style="text-align:right;text-align: center;"> 10 </td>
   <td style="text-align:right;text-align: center;"> 0.4200462 </td>
   <td style="text-align:right;text-align: center;"> 1.0000 </td>
   <td style="text-align:right;text-align: center;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;text-align: center;"> 5 </td>
   <td style="text-align:right;text-align: center;"> 6731 </td>
   <td style="text-align:right;text-align: center;"> 0.6731 </td>
   <td style="text-align:right;text-align: center;"> 0.3269 </td>
   <td style="text-align:right;text-align: center;"> 4 </td>
   <td style="text-align:right;text-align: center;"> 0.4436500 </td>
   <td style="text-align:right;text-align: center;"> 1.0000 </td>
   <td style="text-align:right;text-align: center;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;text-align: center;"> 6 </td>
   <td style="text-align:right;text-align: center;"> 6300 </td>
   <td style="text-align:right;text-align: center;"> 0.6300 </td>
   <td style="text-align:right;text-align: center;"> 0.3700 </td>
   <td style="text-align:right;text-align: center;"> 5 </td>
   <td style="text-align:right;text-align: center;"> 0.4686667 </td>
   <td style="text-align:right;text-align: center;"> 1.0000 </td>
   <td style="text-align:right;text-align: center;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;text-align: center;"> 14 </td>
   <td style="text-align:right;text-align: center;"> 5543 </td>
   <td style="text-align:right;text-align: center;"> 0.5543 </td>
   <td style="text-align:right;text-align: center;"> 0.4457 </td>
   <td style="text-align:right;text-align: center;"> 13 </td>
   <td style="text-align:right;text-align: center;"> 0.5292687 </td>
   <td style="text-align:right;text-align: center;"> 1.0000 </td>
   <td style="text-align:right;text-align: center;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;text-align: center;"> 2 </td>
   <td style="text-align:right;text-align: center;"> 4957 </td>
   <td style="text-align:right;text-align: center;"> 0.4957 </td>
   <td style="text-align:right;text-align: center;"> 0.5043 </td>
   <td style="text-align:right;text-align: center;"> 1 </td>
   <td style="text-align:right;text-align: center;"> 0.5606056 </td>
   <td style="text-align:right;text-align: center;"> 1.0000 </td>
   <td style="text-align:right;text-align: center;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;text-align: center;"> 3 </td>
   <td style="text-align:right;text-align: center;"> 4689 </td>
   <td style="text-align:right;text-align: center;"> 0.4689 </td>
   <td style="text-align:right;text-align: center;"> 0.5311 </td>
   <td style="text-align:right;text-align: center;"> 2 </td>
   <td style="text-align:right;text-align: center;"> 0.5606056 </td>
   <td style="text-align:right;text-align: center;"> 1.0000 </td>
   <td style="text-align:right;text-align: center;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;text-align: center;"> 12 </td>
   <td style="text-align:right;text-align: center;"> 1107 </td>
   <td style="text-align:right;text-align: center;"> 0.1107 </td>
   <td style="text-align:right;text-align: center;"> 0.8893 </td>
   <td style="text-align:right;text-align: center;"> 11 </td>
   <td style="text-align:right;text-align: center;"> 0.8893000 </td>
   <td style="text-align:right;text-align: center;"> 1.0000 </td>
   <td style="text-align:right;text-align: center;"> 0.0000000 </td>
  </tr>
</tbody>
</table>

The **pval** (p value) indicates that the target gene set (L6 IT CTX) is significantly enriched to the 8, 16, 12, 17, and 7 clusters. More precisely, refer to the **bonferroni** (adjusted p value), the target gene set is significantly enriched in the cluster 8 (pval less than 0.05 is significant). However, should L6 IT CTX cell types be located within cluster 8? We recommend conducting the following visualization analysis to determine the accuracy of enrichment results.<br>

## Step4: Visualizing spatial enrichment results

Using the `SpatialFeaturePlot()` function in [Seurat](https://satijalab.org/seurat/index.html), we can visualize the enrichment results of SpatialE in order to better comprehend their biological significance.

We will compare the SpatialE with the [SPOTlight](https://marcelosua.github.io/SPOTlight/) and [MIA](https://pubmed.ncbi.nlm.nih.gov/31932730/) (Multimodal intersection analysis).The analysis results of SPOTlight and MIA have been packaged in SpatialE.

```r
spot_info <- mouse_brain@meta.data
spot_info[1,ncol(spot_info)+1] <- NA
colnames(spot_info)[ncol(spot_info)] <- 'SpatialE_L6_IT_CTX' # Enter the name you want to define
for (i in 1:length(levels(spot_info$seurat_clusters))) {
  for (j in 1:length(spot_info$seurat_clusters)) {
    if (spot_info$seurat_clusters[j] == enrichment_output$cluster[i]) {
      spot_info[j,ncol(spot_info)] <- enrichment_output$graphdata[i]
    }
  }
}
mouse_brain@meta.data <- spot_info
SpatialFeaturePlot(mouse_brain, features = c('SpatialE_L6_IT_CTX','SPOTlight_L6_IT_CTX','MIA_L6_IT_CTX'), alpha = c(0.3,1))
```

<img src="https://github.com/TJ-zhanglab/SpatialE/blob/main/PNG/mapping_results.png" width="1000"/>

Referring to the [Allen Brain map](http://atlas.brain-map.org/atlas?atlas=1#atlas=1&structure=12998&resolution=16.75&x=5536&y=4146&zoom=-4&plate=100960236&z=5), we can conclude that SpatialE correctly mapped l6 IT CTX cell type to the L6 cortex in the coronal slice of the mouse brain. SPOTlight demonstrated that it would be enriched to the hippocampus and hypothalamus. MIA showed that it was significantly enriched in L2/3 cortex besides L6 cortex. SpatialE achieves more specific and accurate enrichment results than the other two methods.<br>

