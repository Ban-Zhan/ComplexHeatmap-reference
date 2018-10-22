

# A List of Heatmaps {#a-list-of-heatmaps}

The main feature of **ComplexHeatmap** package is it supports to concatenate a list of 
heatmaps/annotations horizontally or vertically so that it makes it possible to visualize
the associations from various sources of information. In this chapter, we mainly introduce
the horizontal concatenation because this is the major case we will use in the analysis. In
the end we show some examples of vertical concatenation. The concept behind basically is similar.

For the horizontal concatenation, the number of rows for all heatmaps/annotations should be the same.
In following we first introduce the concatenation of heatmaps and later we will show how to concatenate
heatmaps with annotations.

In following example, there are three matrices where the third heatmap is a vector and it will
be transformed as a one-column matrix. The one-column heatmap is sometimes useful when you
concatenate a list of heatmaps that it can show e.g. annotations for each row or some scores
of each row. e.g. if rows are genes, the the whether genes are protein coding gene can be 
represented ..., or p-values or foldchange from differnetial...

To concatenate heatmaps, simply use `+` operator.


```r
set.seed(123)
mat1 = matrix(rnorm(80, 2), 8, 10)
mat1 = rbind(mat1, matrix(rnorm(40, -2), 4, 10))
rownames(mat1) = paste0("R", 1:12)
colnames(mat1) = paste0("C", 1:10)

mat2 = matrix(runif(60, max = 3, min = 1), 6, 10)
mat2 = rbind(mat2, matrix(runif(60, max = 2, min = 0), 6, 10))
rownames(mat2) = paste0("R", 1:12)
colnames(mat2) = paste0("C", 1:10)

le = sample(letters[1:3], 12, replace = TRUE)
names(le) = paste0("R", 1:12)

ind = sample(12, 12)
mat1 = mat1[ind, ]
mat2 = mat2[ind, ]
le = le[ind]

ht1 = Heatmap(mat1, name = "rnorm")
ht2 = Heatmap(mat2, name = "runif")
ht3 = Heatmap(le, name = "letters")

ht1 + ht2 + ht3
```

<img src="04-a_list_of_heatmaps_files/figure-html/heatmap_list_default-1.png" width="960" style="display: block; margin: auto;" />

Under default mode, dendrograms from the second heatmap will be removed and
row orders will be same as the first one. also row names which are put on the right
side of the heatmap for hte first two heatmaps are removed as well.

The returned value of addition of two heatmaps is a `HeatmapList` object.
Directly calling `ht_list` object will call `draw()` method with default
settings. With explicitly calling `draw()` method, you can have more controls
e.g. on the legend and titles.


```r
ht_list = ht1 + ht2 + ht3
class(ht_list)
```

```
## [1] "HeatmapList"
## attr(,"package")
## [1] "ComplexHeatmap"
```

You can append any number of heatmaps to the heatmap list. Also you can append a heatmap list to a heatmap list.


```r
ht1 + ht_list
ht_list + ht1
ht_list + ht_list
```

`NULL` can be added to the heatmap list. It would be convinient when users want to construct a heatmap list through a `for` loop.


```r
ht_list = NULL  ## Heatmap(...) + NULL gives you a HeatmapList object
for(s in sth) {
    ht_list = ht_list + Heatmap(...)
}
```

## Titles

A heatmap list also has titles which are independent to the heatmap titles.


```r
col_rnorm = colorRamp2(c(-3, 0, 3), c("green", "white", "red"))
col_runif = colorRamp2(c(0, 3), c("white", "orange"))
col_letters = c("a" = "pink", "b" = "purple", "c" = "blue")
ht1 = Heatmap(mat1, name = "rnorm", col = col_rnorm,
    row_title = "Heatmap 1", column_title = "Heatmap 1")
ht2 = Heatmap(mat2, name = "runif", col = col_runif,
    row_title = "Heatmap 2", column_title = "Heatmap 2")
ht3 = Heatmap(le, name = "letters", col = col_letters)
ht_list = ht1 + ht2 + ht3

draw(ht_list, row_title = "Three heatmaps, row title", row_title_gp = gpar(col = "red"),
    column_title = "Three heatmaps, column title", column_title_gp = gpar(fontsize = 16))
```

<img src="04-a_list_of_heatmaps_files/figure-html/heatmap_list_title-1.png" width="960" style="display: block; margin: auto;" />


## Size of heatmaps

The width for some (not all) heatmaps can be set to a fixed width.


```r
ht2 = Heatmap(mat2, name = "runif", col = col_runif, width = unit(4, "cm"))
ht3 = Heatmap(le, name = "letters", col = col_letters, width = unit(5, "mm"))
ht1 + ht2 + ht3
```

<img src="04-a_list_of_heatmaps_files/figure-html/heatmap_list_size-1.png" width="960" style="display: block; margin: auto;" />

or the width can be set as relative values. Please not in this case, `width` for all heatmaps
should be set (relative width and fixed width can be mixed).


```r
ht1 = Heatmap(mat1, name = "rnorm", col = col_rnorm, width = unit(4, "cm"))
ht2 = Heatmap(mat2, name = "runif", col = col_runif, width = unit(6, "cm"))
ht3 = Heatmap(le, name = "letters", col = col_letters, width = unit(1, "cm"))
ht1 + ht2 + ht3
```

<img src="04-a_list_of_heatmaps_files/figure-html/heatmap_list_relative_size-1.png" width="960" style="display: block; margin: auto;" />

```
## Since all heatmaps/annotations have absolute units, the total width of the plot is 151mm
```


```r
ht1 = Heatmap(mat1, name = "rnorm", col = col_rnorm,
    show_row_names = FALSE, width = 6)
ht2 = Heatmap(mat2, name = "runif", col = col_runif,
    show_row_names = FALSE, width = 4)
ht3 = Heatmap(le, name = "letters", col = col_letters, width = 1)
ht1 + ht2 + ht3
```

<img src="04-a_list_of_heatmaps_files/figure-html/unnamed-chunk-5-1.png" width="960" style="display: block; margin: auto;" />

## Gap between heatmaps



```r
ht1 = Heatmap(mat1, name = "rnorm", col = col_rnorm)
ht2 = Heatmap(mat2, name = "runif", col = col_runif)
ht3 = Heatmap(le, name = "letters", col = col_letters)
draw(ht_list, ht_gap = unit(1, "cm"))
```

<img src="04-a_list_of_heatmaps_files/figure-html/heatmap_list_gap-1.png" width="960" style="display: block; margin: auto;" />

```r
draw(ht_list, ht_gap = unit(c(3, 10), "mm"))
```

<img src="04-a_list_of_heatmaps_files/figure-html/heatmap_list_gap-2.png" width="960" style="display: block; margin: auto;" />

## Auto adjustment to the main heatmap

There are some automatic adjustment if more than one heatmaps are plotted. There should be a main
heatmap which by default is the first one. Some settings for the remaining heatmaps will be modified
to the settings in the main heatmap. The adjustment are:

- row clusters are removed.
- row titles are removed.
- if the main heatmap is split by rows, all remaining heatmaps will also be split by same levels as
  the main one.

The main heatmap can be specified by `main_heatmap` argument. The value can be a numeric index or
the name of the heatmap (of course, you need to set the heatmap name when you create the `Heatmap`
object).


```r
ht1 = Heatmap(mat1, name = "rnorm", col = col_rnorm, km = 2)
ht2 = Heatmap(mat2, name = "runif", col = col_runif)
ht3 = Heatmap(le, name = "letters", col = col_letters)
```



```r
ht2 + ht1 + ht3
```

<img src="04-a_list_of_heatmaps_files/figure-html/unnamed-chunk-6-1.png" width="480" style="display: block; margin: auto;" />

here although `ht1` is the second heatmap, we specify `ht1` to be
the main heatmap by explicitely setting `main_heatmap` argument


```r
ht_list = ht2 + ht1 + ht3
draw(ht_list, main_heatmap = "rnorm")
```

<img src="04-a_list_of_heatmaps_files/figure-html/unnamed-chunk-7-1.png" width="480" style="display: block; margin: auto;" />


```r
ht_list = ht2 + ht1 + ht3
draw(ht_list, main_heatmap = "rnorm", row_dend_side = "right", row_sub_title_side = "left")
```

<img src="04-a_list_of_heatmaps_files/figure-html/unnamed-chunk-8-1.png" width="480" style="display: block; margin: auto;" />

If there is no row clustering in the main heatmap, all other heatmaps have no row clustering neither.


```r
ht1 = Heatmap(mat1, name = "rnorm", col = col_rnorm, cluster_rows = FALSE)
ht2 = Heatmap(mat2, name = "runif", col = col_runif)
ht3 = Heatmap(le, name = "letters", col = col_letters)
ht1 + ht2 + ht3
```

<img src="04-a_list_of_heatmaps_files/figure-html/heatmap_list_auto_adjust_no_row_cluster-1.png" width="960" style="display: block; margin: auto;" />

## control row ... in draw() function


Since the main heamtap controls the row order of all heatmaps, the parameters which co...


- cluster_rows
- clustering_distance_rows
- clustering_method_rows 
- row_dend_width
- show_row_dend
- row_dend_reorder
- row_dend_gp
- row_order 

And for splitting rows 

- row_gap
- row_km
- row_split


```r
ht1 = Heatmap(mat1, name = "rnorm", col = col_rnorm, row_km = 2)
ht2 = Heatmap(mat2, name = "runif", col = col_runif)
ht3 = Heatmap(le, name = "letters", col = col_letters)
ht_list = ht1 + ht2 + ht3
draw(ht_list, row_km = 1, row_split = le, cluster_rows = FALSE)
```

<img src="04-a_list_of_heatmaps_files/figure-html/unnamed-chunk-9-1.png" width="480" style="display: block; margin: auto;" />

## Annotation as components are adjusted


```r
ha1 = HeatmapAnnotation(foo1 = 1:10, annotation_name_side = "left")
ht1 = Heatmap(mat1, name = "rnorm", col = col_rnorm, top_annotation = ha1)
ht2 = Heatmap(mat2, name = "runif", col = col_runif)
ht3 = Heatmap(le, name = "letters", col = col_letters)
ht1 + ht2 + ht3
```

<img src="04-a_list_of_heatmaps_files/figure-html/unnamed-chunk-10-1.png" width="480" style="display: block; margin: auto;" />


```r
ha1 = HeatmapAnnotation(foo1 = 1:10, bar1 = anno_points(1:10), annotation_name_side = "left")
ha2 = HeatmapAnnotation(bar2 = anno_barplot(1:10))
ht1 = Heatmap(mat1, name = "rnorm", col = col_rnorm, top_annotation = ha1)
ht2 = Heatmap(mat2, name = "runif", col = col_runif, top_annotation = ha2)
ht3 = Heatmap(le, name = "letters", col = col_letters)
ht_list = ht1 + ht2 + ht3
draw(ht_list, ht_gap = unit(c(6, 2), "mm"))
```

<img src="04-a_list_of_heatmaps_files/figure-html/unnamed-chunk-11-1.png" width="480" style="display: block; margin: auto;" />


```r
ha1 = HeatmapAnnotation(foo1 = 1:10, annotation_name_side = "left")
ha2 = HeatmapAnnotation(bar2 = anno_barplot(1:10))
ht1 = Heatmap(mat1, name = "rnorm", col = col_rnorm, top_annotation = ha1)
ht2 = Heatmap(mat2, name = "runif", col = col_runif, top_annotation = ha2)
ht3 = Heatmap(le, name = "letters", col = col_letters)
ht_list = ht1 + ht2 + ht3
draw(ht_list, ht_gap = unit(c(6, 2), "mm"))
```

<img src="04-a_list_of_heatmaps_files/figure-html/unnamed-chunk-12-1.png" width="480" style="display: block; margin: auto;" />


```r
ha1 = HeatmapAnnotation(foo1 = 1:10, bar1 = anno_points(1:10), annotation_name_side = "left")
ht1 = Heatmap(mat1, name = "rnorm", col = col_rnorm, bottom_annotation = ha1)
ht2 = Heatmap(mat2, name = "runif", col = col_runif)
ht3 = Heatmap(le, name = "letters", col = col_letters)
ht_list = ht1 + ht2 + ht3
draw(ht_list)
```

<img src="04-a_list_of_heatmaps_files/figure-html/unnamed-chunk-13-1.png" width="480" style="display: block; margin: auto;" />


## concatenate with annotations

For horizontal concatenation with the annotaions.



```r
ha1 = rowAnnotation(foo = 1:12, bar = anno_barplot(1:12, width = unit(4, "cm")))
ht1 = Heatmap(mat1, name = "rnorm", col = col_rnorm, row_km = 2)
ht1 + ha1
```

<img src="04-a_list_of_heatmaps_files/figure-html/unnamed-chunk-14-1.png" width="480" style="display: block; margin: auto;" />


```r
Heatmap(mat1, name = "rnorm", col = col_rnorm, row_km = 2) + 
    rowAnnotation(foo = 1:12) +
    rowAnnotation(bar = anno_barplot(1:12, width = unit(4, "cm")))
```

<img src="04-a_list_of_heatmaps_files/figure-html/unnamed-chunk-15-1.png" width="480" style="display: block; margin: auto;" />


```r
Heatmap(mat1, name = "rnorm", col = col_rnorm, row_km = 2) + 
    rowAnnotation(foo = 1:12) +
    rowAnnotation(bar = anno_barplot(1:12, width = unit(4, "cm"))) +
    Heatmap(mat2, name = "runif", col = col_runif)
```

<img src="04-a_list_of_heatmaps_files/figure-html/unnamed-chunk-16-1.png" width="480" style="display: block; margin: auto;" />


## Only annotations


```r
rowAnnotation(foo = 1:12) +
    rowAnnotation(bar = anno_barplot(1:12, width = unit(4, "cm")))
```

<img src="04-a_list_of_heatmaps_files/figure-html/unnamed-chunk-17-1.png" width="480" style="display: block; margin: auto;" />

```
## Since all heatmaps/annotations have absolute units, the total width of the plot is 62mm
```


```r
rowAnnotation(bar = anno_barplot(1:12, width = unit(4, "cm"))) + NULL
```

<img src="04-a_list_of_heatmaps_files/figure-html/unnamed-chunk-18-1.png" width="480" style="display: block; margin: auto;" />

```
## Since all heatmaps/annotations have absolute units, the total width of the plot is 44mm
```

## vertical align


```r
mat1t = t(mat1)
mat2t = t(mat2)
ht1 = Heatmap(mat1t, name = "rnorm", col = col_rnorm)
ht2 = Heatmap(mat2t, name = "runif", col = col_runif)
ht3 = Heatmap(rbind(letters = le), name = "letters", col = col_letters)
ht_list = ht1 %v% ht2 %v% ht3
draw(ht_list)
```

<img src="04-a_list_of_heatmaps_files/figure-html/unnamed-chunk-19-1.png" width="480" style="display: block; margin: auto;" />

```r
draw(ht_list, column_km = 2)
```

<img src="04-a_list_of_heatmaps_files/figure-html/unnamed-chunk-19-2.png" width="480" style="display: block; margin: auto;" />


```r
ht1 = Heatmap(mat1t, name = "rnorm", col = col_rnorm)
ht2 = Heatmap(mat2t, name = "runif", col = col_runif)
ht3 = Heatmap(rbind(letters = le), name = "letters", col = col_letters)
ha = HeatmapAnnotation(foo = anno_barplot(1:12, height = unit(2, "cm")))
ht_list = ht1 %v% ha %v% ht2 %v% ht3
draw(ht_list, column_km = 2)
```

<img src="04-a_list_of_heatmaps_files/figure-html/unnamed-chunk-20-1.png" width="480" style="display: block; margin: auto;" />



```r
ht1 = Heatmap(mat1t, name = "rnorm", col = col_rnorm, row_km = 2)
ht2 = Heatmap(mat2t, name = "runif", col = col_runif, row_km = 2)
ht3 = Heatmap(rbind(letters = le), name = "letters", col = col_letters)
ha = HeatmapAnnotation(foo = anno_barplot(1:12, height = unit(2, "cm")))
ht_list = ht1 %v% ha %v% ht2 %v% ht3
draw(ht_list, column_km = 2)
```

<img src="04-a_list_of_heatmaps_files/figure-html/unnamed-chunk-21-1.png" width="480" style="display: block; margin: auto;" />

## Retrieve orders and dendrograms

`row_order`, `column_order`, `row_dend` and `column_dend` can be used to retrieve corresponding information from
the heatmaps. The usage is straightforward by following example:


```r
ht1 = Heatmap(mat1, name = "rnorm", col = col_rnorm)
ht2 = Heatmap(mat2, name = "runif", col = col_runif)
ht_list = ht1 + ht2
ht_list = draw(ht_list)
```

<img src="04-a_list_of_heatmaps_files/figure-html/unnamed-chunk-22-1.png" width="480" style="display: block; margin: auto;" />

```r
row_order(ht_list)
```

```
##  [1]  8  6 10 11  5  7  9  2 12  4  1  3
```

```r
column_order(ht_list)
```

```
## $rnorm
##  [1]  2  7  5  6 10  1  9  8  4  3
## 
## $runif
##  [1]  1  3  8  7  6  9  4 10  2  5
```


```r
ht1 = Heatmap(mat1, name = "rnorm", col = col_rnorm)
ht2 = Heatmap(mat2, name = "runif", col = col_runif, column_km = 2)
ht_list = ht1 + ht2
ht_list = draw(ht_list, row_km = 2)
```

<img src="04-a_list_of_heatmaps_files/figure-html/unnamed-chunk-23-1.png" width="480" style="display: block; margin: auto;" />

```r
row_order(ht_list)
```

```
## $`1`
## [1] 12  4  1  3
## 
## $`2`
## [1]  8  6 10 11  5  7  9  2
```

```r
column_order(ht_list)
```

```
## $rnorm
## $rnorm[[1]]
##  [1]  2  7  5  6 10  1  9  8  4  3
## 
## 
## $runif
## $runif$`1`
## [1] 1 3 8 7 6 9
## 
## $runif$`2`
## [1]  4 10  2  5
```

Same logic for vertical ... which we will not show here

## Change graphic parameters simultaneously

`ht_opt()` can set graphic parameters for dimension names and titles as global settings.


```r
ht_opt
```

```
##                   Option Value
##     heatmap_row_names_gp  NULL
##  heatmap_column_names_gp  NULL
##     heatmap_row_title_gp  NULL
##  heatmap_column_title_gp  NULL
##          legend_title_gp  NULL
##    legend_title_position  NULL
##         legend_labels_gp  NULL
##       legend_grid_height  NULL
##        legend_grid_width  NULL
##            legend_border  NULL
##           heatmap_border  NULL
##        annotation_border  NULL
##              fast_hclust FALSE
##                  verbose FALSE
##                  show_vp FALSE
##         anno_simple_size   5mm
##       DENDROGRAM_PADDING 0.5mm
##          DIMNAME_PADDING   1mm
##            TITLE_PADDING 2.5mm
##      COLUMN_ANNO_PADDING   1mm
##         ROW_ANNO_PADDING   1mm
```


```r
ht_opt(heatmap_column_names_gp = gpar(fontface = "italic"), 
    heatmap_column_title_gp = gpar(fontsize = 10),
    legend_border = "black",
    heatmap_border = TRUE,
    annotation_border = TRUE
)
ht1 = Heatmap(mat1, name = "ht1", column_title = "Heatmap 1",
    top_annotation = HeatmapAnnotation(foo = 1:10))
ht2 = Heatmap(mat2, name = "ht2", column_title = "Heatmap 2",
    top_annotation = HeatmapAnnotation(bar = 1:10))
ht1 + ht2
```

<img src="04-a_list_of_heatmaps_files/figure-html/unnamed-chunk-25-1.png" width="960" style="display: block; margin: auto;" />

```r
ht_opt(RESET = TRUE)
```

Following are global settings supported by `ht_global_opt()`. By this function, you can also control settings
for the legends.

## Session info


```r
sessionInfo()
```

```
## R version 3.4.4 (2018-03-15)
## Platform: x86_64-apple-darwin15.6.0 (64-bit)
## Running under: macOS High Sierra 10.13.2
## 
## Matrix products: default
## BLAS: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRblas.0.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/3.4/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
## 
## attached base packages:
## [1] grid      methods   stats     graphics  grDevices utils     datasets 
## [8] base     
## 
## other attached packages:
## [1] dendsort_0.3.3        dendextend_1.8.0      circlize_0.4.5       
## [4] ComplexHeatmap_1.99.0
## 
## loaded via a namespace (and not attached):
##  [1] shape_1.4.4         modeltools_0.2-22   GetoptLong_0.1.7   
##  [4] tidyselect_0.2.5    xfun_0.3            kernlab_0.9-27     
##  [7] purrr_0.2.5         lattice_0.20-35     colorspace_1.3-2   
## [10] viridisLite_0.3.0   htmltools_0.3.6     stats4_3.4.4       
## [13] yaml_2.2.0          rlang_0.2.2         pillar_1.3.0       
## [16] glue_1.3.0          prabclus_2.2-6      RColorBrewer_1.1-2 
## [19] fpc_2.1-11.1        bindrcpp_0.2.2      plyr_1.8.4         
## [22] bindr_0.1.1         robustbase_0.93-3   stringr_1.3.1      
## [25] munsell_0.5.0       gtable_0.2.0        mvtnorm_1.0-8      
## [28] GlobalOptions_0.1.1 evaluate_0.12       knitr_1.20         
## [31] flexmix_2.3-14      class_7.3-14        DEoptimR_1.0-8     
## [34] trimcluster_0.1-2.1 Rcpp_0.12.19        scales_1.0.0       
## [37] backports_1.1.2     diptest_0.75-7      gridExtra_2.3      
## [40] rjson_0.2.20        ggplot2_3.0.0       digest_0.6.18      
## [43] stringi_1.2.4       bookdown_0.7        dplyr_0.7.7        
## [46] rprojroot_1.3-2     tools_3.4.4         magrittr_1.5       
## [49] lazyeval_0.2.1      tibble_1.4.2        cluster_2.0.7-1    
## [52] crayon_1.3.4        whisker_0.3-2       pkgconfig_2.0.2    
## [55] MASS_7.3-51         viridis_0.5.1       assertthat_0.2.0   
## [58] rmarkdown_1.10      R6_2.3.0            mclust_5.4.1       
## [61] nnet_7.3-12         compiler_3.4.4
```
