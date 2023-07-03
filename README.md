
**Hi**erarchical  **Decon**volution  (**HiDecon)**
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
<img src="tree.png" style="width:35.0%" alt="Hierarchical tree constructed from cell lineage relationship and used to guide HiDecon." /><figcaption aria-hidden="true">Hierarchical tree constructed from cell lineage relationship and used to guide HiDecon.</figcaption>
</figure>

### Example data

-   bulk.dat: Gene by sample matrix. Simulated bulk tissue data of 126
    samples and 1572 genes.

-   ref.dat: Gene by cell matrix. This single cell reference used has
    100 cells and 1572 genes.

-   ref.type: Cell type labels of the single cell data ref.dat.

-   B: List of the same length of the hierarchical tree. It contains
    cell type mapping matrices from layer 0 to 1 and L-1 to L.

-   order.type: the cell type order of bottom layer of the given
    hierarchical tree.

``` r
data("B", "bulk.dat", "order.type", "ref.dat", "ref.type", package = "HiDecon")
```

### HiDecon using prespecified parameter

``` r
example.res <- HiDecon(bulk = bulk.dat, ref = ref.dat, B = B, cell_type = ref.type, type_order = order.type,
                    lambda = 40)
```

### HiDecon using selected parameter

``` r
select.res <- select_HiDecon(bulk = bulk.dat, ref = ref.dat, B =B, cell_type = ref.type, type_order = order.type)
```

Check *λ* chosen by tuning parameter selection method and the mCCC on
different *λ*s.

``` r
select.res$lambda
#> [1] 40
select.res$mCCC
#>        10        20        30        40        50        60        70        80 
#> 0.8693923 0.8882590 0.8959806 0.8968308 0.8936180 0.8890917 0.8840554 0.8790285 
#>        90       100       110       120       130       140       150       160 
#> 0.8742005 0.8695607 0.8651292 0.8609163 0.8569639 0.8532691 0.8497291 0.8462979 
#>       170       180       190       200 
#> 0.8430338 0.8398991 0.8368969 0.8339872
```
