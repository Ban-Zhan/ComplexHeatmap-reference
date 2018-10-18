
# Other High-level Plots {#other-high-level-plots}

## Density heatmap {#density-heatmap}



```r
set.seed(123)
m = cbind(matrix(rnorm(10*100), ncol = 10),
	      matrix(runif(10*100, min = -2, max = 2) + 0.5, ncol = 10))
colnames(m) = paste0("C", 1:ncol(m))
densityHeatmap(m)
```

<img src="08-other-high-level-plots_files/figure-html/unnamed-chunk-2-1.png" width="576" style="display: block; margin: auto;" />


```r
densityHeatmap(m, ylim = c(-2, 2), title = "Distribution as heatmap", ylab = "some values")
```

<img src="08-other-high-level-plots_files/figure-html/unnamed-chunk-3-1.png" width="576" style="display: block; margin: auto;" />


```r
densityHeatmap(m, column_order = sample(20, 20))
```

<img src="08-other-high-level-plots_files/figure-html/unnamed-chunk-4-1.png" width="576" style="display: block; margin: auto;" />



```r
densityHeatmap(m, col = topo.colors(10))
```

<img src="08-other-high-level-plots_files/figure-html/unnamed-chunk-5-1.png" width="576" style="display: block; margin: auto;" />



```r
ha1 = HeatmapAnnotation(dist = c(rep("rnorm", 10), rep("runif", 10)))
ha2 = HeatmapAnnotation(foo = anno_points(rnorm(20)))
densityHeatmap(m, top_annotation = ha1, bottom_annotation = ha2)
```

<img src="08-other-high-level-plots_files/figure-html/unnamed-chunk-6-1.png" width="576" style="display: block; margin: auto;" />















