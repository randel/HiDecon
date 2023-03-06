# HiDecon
Hierarchical cellular deconvolution
![\textbf{Hi}erarchical \\ \\ \textbf{Decon}volution \\ \\ (\textbf{HiDecon)}](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctextbf%7BHi%7Derarchical%20%5C%20%5C%20%5Ctextbf%7BDecon%7Dvolution%20%5C%20%5C%20%28%5Ctextbf%7BHiDecon%29%7D "\textbf{Hi}erarchical \ \ \textbf{Decon}volution \ \ (\textbf{HiDecon)}")
is a cellular deconvolution method that uses single-cell RNA sequencing
references and a hierarchical cell type tree, which models the
similarities among cell types and cell differentiation relationships, to
estimate cellular fractions in bulk data. HiDecon has outstanding
performance on estimating related cell types as well as rare cell types.

## Installation

``` r
if (!"devtools" %in% rownames(installed.packages())) {
  install.packages('devtools')
}

# install the HiDecon R package
library(devtools)
if (!"HiDecon" %in% rownames(installed.packages())) {
  install_github('randel/HiDecon')
}


# load
library(HiDecon)
```

## Examples

In this part, we will show how to use HiDecon package using example data
in the package.

### Hierarchical tree and the cell type mapping matrix

Hierarchical tree is specified from well-studied cell lineages or can be
learnt from hierarchical clustering of scRNA-seq data. The hierarchical
tree used to guide the deconvolution of example data is shown as below:

<figure>
<embed src="tree.pdf" style="width:35.0%" />
<figcaption aria-hidden="true">Hierarchical tree constructed from cell
lineage relationship and used to guide HiDecon.</figcaption>
</figure>

### Example data
