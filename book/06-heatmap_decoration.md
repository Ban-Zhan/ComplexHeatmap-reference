

# Heatmap Decoration {#heatmap-decoration}

Each component of the heatmap/heatmap list has a name (unique id). You can go to any viewport 
to add graphics in by specifying the heatmap/annotation name.

First generate a figure that almost contains all types of heatmap components.


```r
mat = matrix(rnorm(80, 2), 8, 10)
mat = rbind(mat, matrix(rnorm(40, -2), 4, 10))
rownames(mat) = paste0("R", 1:12)
colnames(mat) = paste0("C", 1:10)

ha_column1 = HeatmapAnnotation(points = anno_points(rnorm(10)), 
    annotation_name_side = "left")
ht1 = Heatmap(mat, name = "ht1", km = 2, column_title = "Heatmap 1", 
    top_annotation = ha_column1, row_names_side = "left")

ha_column2 = HeatmapAnnotation(type = c(rep("a", 5), rep("b", 5)),
    col = list(type = c("a" = "red", "b" = "blue")))
ht2 = Heatmap(mat, name = "ht2", row_title = "Heatmap 2", column_title = "Heatmap 2",
    bottom_annotation = ha_column2, column_km = 2)

ht_list = ht1 + ht2 + rowAnnotation(bar = anno_barplot(rowMeans(mat), width = unit(2, "cm")))
draw(ht_list, row_title = "Heatmap list", column_title = "Heatmap list")
list_components()
```

```
##  [1] "ROOT"                         "global"                      
##  [3] "global_layout"                "main_heatmap_list"           
##  [5] "heatmap_ht1"                  "ht1_heatmap_body_wrap"       
##  [7] "ht1_heatmap_body_1_1"         "ht1_heatmap_body_2_1"        
##  [9] "ht1_column_title_1"           "ht1_row_title_1"             
## [11] "ht1_row_title_2"              "ht1_dend_row_1"              
## [13] "ht1_dend_row_2"               "ht1_dend_column_1"           
## [15] "ht1_row_names_1"              "ht1_row_names_2"             
## [17] "ht1_column_names_1"           "annotation_points_1"         
## [19] "heatmap_ht2"                  "ht2_heatmap_body_wrap"       
## [21] "ht2_heatmap_body_1_1"         "ht2_heatmap_body_1_2"        
## [23] "ht2_heatmap_body_2_1"         "ht2_heatmap_body_2_2"        
## [25] "ht2_column_title_1"           "ht2_dend_column_1"           
## [27] "ht2_dend_column_2"            "ht2_column_names_1"          
## [29] "ht2_column_names_2"           "annotation_type_1"           
## [31] "annotation_type_2"            "heatmap_heatmap_annotation_2"
## [33] "annotation_bar_1"             "annotation_bar_2"            
## [35] "global_column_title"          "global_row_title"            
## [37] "heatmap_legend"               "annotation_legend"
```

<img src="06-heatmap_decoration_files/figure-html/access_components-1.png" width="960" style="display: block; margin: auto;" />

The components (viewports) that have names are:

- `heatmap_@{heatmap_name}_@{i}_@{j}`: the viewport which contains a single heatmap
- `annotation_@{annotation_name}_@{i}`: for row annotations
- `@{heatmap_name}_heatmap_body_@{i}`: the heatmap body.
- `@{heatmap_name}_column_title_@{i}`: column title for a single heatmap.
- `@{heatmap_name}_row_title_@{i}`: since a heatmap body may be splitted into several parts. `@{i}` is the index of the row slice.
- `@{heatmap_name}_dend_row_@{i}`: dendrogram for ith row slice.
- `@{heatmap_name}_dend_column_@{i}`: dendrogram on columns
- `@{heatmap_name}_row_names_@{i}`: the viewport which contains row names.
- `@{heatmap_name}_column_names_@{i}`: the viewport which contains column names.

## decorate_*() functions

Basically, you can go to these components by `seekViewport()`, but to hide the
details that is too low-level, **ComplexHeatmap** package provides
`decorate_*` family functions which makes it easy to add graphics into
different components.

Following code add annotation names, mark one grid in the heatmap and seperate the first column clusters with two rectangles.


```r
ht_list = draw(ht_list, row_title = "Heatmap list", column_title = "Heatmap list", 
    heatmap_legend_side = "right", annotation_legend_side = "left")

decorate_heatmap_body("ht1", {
    grid.text("outlier", 1.5/10, 2.5/4, default.units = "npc")
    grid.lines(c(0.5, 0.5), c(0, 1), gp = gpar(lty = 2, lwd = 2))
}, slice = 2)

decorate_column_dend("ht1", {
    tree = column_dend(ht_list)$ht1[[1]]
    ind = cutree(as.hclust(tree), k = 2)[order.dendrogram(tree)]

    first_index = function(l) which(l)[1]
    last_index = function(l) { x = which(l); x[length(x)] }
    x1 = c(first_index(ind == 1), first_index(ind == 2)) - 1
    x2 = c(last_index(ind == 1), last_index(ind == 2))
    grid.rect(x = x1/length(ind), width = (x2 - x1)/length(ind), just = "left",
        default.units = "npc", gp = gpar(fill = c("#FF000040", "#00FF0040"), col = NA))
})

decorate_row_names("ht1", {
    grid.rect(gp = gpar(fill = "#FF000040"))
}, slice = 2)

decorate_row_title("ht1", {
    grid.rect(gp = gpar(fill = "#00FF0040"))
}, slice = 1)

decorate_annotation("points", {
    grid.lines(c(0, 1), unit(c(0, 0), "native"), gp = gpar(col = "red"))
})
```

<img src="06-heatmap_decoration_files/figure-html/components-1.png" width="960" style="display: block; margin: auto;" />

For annotations which are created by `anno_points()`, `anno_barplot()` and `anno_boxplot()`, "native" unit
can be used in the decoration code.

## Examples

### Barplot for single column heatmap


```r
library(circlize)
bed = generateRandomBed(nr = 1000)
prop = runif(nrow(bed))
col_fun = colorRamp2(c(0, 1), c("white", "orange"))
split = c(rep("group1", 400), rep("group2", nrow(bed) - 400))
ht = Heatmap(prop, name = "prop", col = col_fun, width = unit(2, "cm"),
    top_annotation = HeatmapAnnotation(barplot = anno_empty(height = unit(4, "cm"))))
ht = draw(ht, row_split = split)
```

```
## Since all heatmaps/annotations have absolute units, the total width of the plot is 56mm
```

```r
ro = row_order(ht)
w = bed[, 3] - bed[, 2]
p = sapply(ro, function(index) {
    sum(w[index]*prop[index])/sum(w[index])
})
decorate_annotation("barplot", {
    pushViewport(viewport(xscale = c(0.5, 2.5), yscale = c(0, max(p)*1.1)))
    grid.rect(x = 1, y = 0, width = 0.8, height = p[1], just = "bottom",
        gp = gpar(fill = "orange"), default.units = "native")
    grid.rect(x = 2, y = 0, width = 0.8, height = p[2], just = "bottom",
        gp = gpar(fill = "orange"), default.units = "native")
    grid.yaxis()
    grid.text("mean proprotion", x = unit(-1.5, "cm"),rot = 90, just = "bottom")
    popViewport()
})
```

<img src="06-heatmap_decoration_files/figure-html/unnamed-chunk-2-1.png" width="576" style="display: block; margin: auto;" />

### add decoration for a customized heatmap


```r
customized_heatmap = function(m, ...) {
    post_fun = function(ht) {
        m = ht@matrix
        co = column_order(ht)
        cm = colMeans(m)[co]
        nc = ncol(m)
        ind = which(cm < 0)
        if(length(ind)) {
            decorate_heatmap_body(ht@name, {
                grid.rect((ind - 1)/nc, 0.5, 1/nc, 1, gp = gpar(col = "black", lwd = 2, fill = NA),
                    just = "left")
            })
        }
    }
    x = colMeans(m)
    Heatmap(m, ...,
        top_annotation = HeatmapAnnotation(mean = anno_barplot(x, gp = gpar(fill = ifelse(x > 0, "red", "green")), width = unit(4, "cm"))),
        post_fun = post_fun)
}
m1 = matrix(rnorm(100), nr = 10)
m2 = matrix(rnorm(100), nr = 10)
customized_heatmap(m1) %v% customized_heatmap(m2)
```

<img src="06-heatmap_decoration_files/figure-html/unnamed-chunk-3-1.png" width="576" style="display: block; margin: auto;" />

### Other probably usage of decoration

