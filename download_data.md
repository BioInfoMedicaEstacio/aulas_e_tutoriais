R Notebook
================

<http://genomicsclass.github.io/book/pages/GEOquery.html>

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("GEOquery")
```

    ## Bioconductor version 3.11 (BiocManager 1.30.10), R 4.0.0 (2020-04-24)

    ## Installing package(s) 'GEOquery'

    ## package 'GEOquery' successfully unpacked and MD5 sums checked
    ## 

``` r
library(GEOquery)
```

    ## Loading required package: Biobase

    ## Loading required package: BiocGenerics

    ## Loading required package: parallel

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:parallel':
    ## 
    ##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
    ##     clusterExport, clusterMap, parApply, parCapply, parLapply,
    ##     parLapplyLB, parRapply, parSapply, parSapplyLB

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, append, as.data.frame, basename, cbind, colnames,
    ##     dirname, do.call, duplicated, eval, evalq, Filter, Find, get, grep,
    ##     grepl, intersect, is.unsorted, lapply, Map, mapply, match, mget,
    ##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
    ##     rbind, Reduce, rownames, sapply, setdiff, sort, table, tapply,
    ##     union, unique, unsplit, which, which.max, which.min

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Setting options('download.file.method.GEOquery'='auto')

    ## Setting options('GEOquery.inmemory.gpl'=FALSE)

``` r
### This will download a ~70 Mb
gse <- getGEO("GSE21942", GSEMatrix = TRUE)
```

    ## Found 1 file(s)

    ## GSE21942_series_matrix.txt.gz

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   ID_REF = col_character()
    ## )

    ## See spec(...) for full column specifications.

    ## File stored at:

    ## C:\Users\Joao\AppData\Local\Temp\Rtmpq4SXEL/GPL570.soft

    ## Warning: 62 parsing failures.
    ##   row     col           expected    actual         file
    ## 54614 SPOT_ID 1/0/T/F/TRUE/FALSE --Control literal data
    ## 54615 SPOT_ID 1/0/T/F/TRUE/FALSE --Control literal data
    ## 54616 SPOT_ID 1/0/T/F/TRUE/FALSE --Control literal data
    ## 54617 SPOT_ID 1/0/T/F/TRUE/FALSE --Control literal data
    ## 54618 SPOT_ID 1/0/T/F/TRUE/FALSE --Control literal data
    ## ..... ....... .................. ......... ............
    ## See problems(...) for more details.

``` r
show(gse)
```

    ## $GSE21942_series_matrix.txt.gz
    ## ExpressionSet (storageMode: lockedEnvironment)
    ## assayData: 54675 features, 29 samples 
    ##   element names: exprs 
    ## protocolData: none
    ## phenoData
    ##   sampleNames: GSM545818 GSM545819 ... GSM545846 (29 total)
    ##   varLabels: title geo_accession ... disease state:ch1 (33 total)
    ##   varMetadata: labelDescription
    ## featureData
    ##   featureNames: 1007_s_at 1053_at ... AFFX-TrpnX-M_at (54675 total)
    ##   fvarLabels: ID GB_ACC ... Gene Ontology Molecular Function (16 total)
    ##   fvarMetadata: Column Description labelDescription
    ## experimentData: use 'experimentData(object)'
    ##   pubMedIds: 22021740 
    ## Annotation: GPL570

O comando abaixo fará o download do dado cru referente ao experimento
GSE21653.

<https://www.affymetrix.com/support/developer/powertools/changelog/gcos-agcc/cel.html>

``` r
filePaths = getGEOSuppFiles("GSE21942")
filePaths
```

<https://baderlab.github.io/Cytoscape_workflows/EnrichmentMapPipeline/supplemental_protocol2_microarray.html>

``` r
library("Biobase")

library(limma)
```

    ## 
    ## Attaching package: 'limma'

    ## The following object is masked from 'package:BiocGenerics':
    ## 
    ##     plotMA
