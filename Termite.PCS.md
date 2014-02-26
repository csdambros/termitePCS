# This file describes the analysis of the termite phylogenetic community structure in the Amazonian Forest

### Load Packages


```r

library(ape)
library(picante)
```

```
## Loading required package: vegan
```

```
## Loading required package: permute
```

```
## This is vegan 2.0-7
```

```
## Loading required package: nlme
```


### Import data


```r

isoptera <- read.csv("Isoptera.csv", row.names = 1)

# Removes repeated records of the same species in the same subplot (5x2
# quadrant)

isoptera <- isoptera[!duplicated(isoptera[, c("Taxon", "Location", "Grid", "Trail", 
    "Plot", "Distance_beginning", "Side")]), ]

```



### Uses taxonomic information to build the basic phylo tree


```r

sp.ls <- unique(isoptera[, c("Taxon", "Genus", "Subfamily", "Family", "Order")])

sp.ls.int <- as.matrix(data.frame(lapply(sp.ls, as.integer)))

sp.ls.int[is.na(sp.ls.int)] <- 0

# Calculates the taxonomic distance among all species pairs (the code can
# be improved - get rid of the for loop )

sp.ls.dist <- (sp.ls.int[, 1]) %*% t(sp.ls.int[, 1] * 0 + 1) != (sp.ls.int[, 
    1] * 0 + 1) %*% t(sp.ls.int[, 1])

for (i in 2:ncol(sp.ls.int)) {
    sp.ls.dist <- as.matrix(sp.ls.dist) + as.matrix((sp.ls.int[, i]) %*% t(sp.ls.int[, 
        i] * 0 + 1) != (sp.ls.int[, i] * 0 + 1) %*% t(sp.ls.int[, i]))
}

dimnames(sp.ls.dist) <- list(sp.ls[, 1], sp.ls[, 1])

write.table(sp.ls.dist, "taxodist.csv")

sp.ls.dist <- as.dist(sp.ls.dist)

termite.taxophylo <- as.phylo(hclust(sp.ls.dist))

write.tree(termite.taxophylo, "termite.taxophylo.tre")

# cophenetic(termite.taxophylo) returns the original taxonomic distance
```




```r

plot(termite.taxophylo)
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1.png) 




```r

attach(isoptera)


termite.plot <- tapply(Frequency, list(paste(Location, Grid, Trail, Plot), Taxon), 
    sum)
termite.plot[is.na(termite.plot)] <- 0

termite.plot.PA <- (termite.plot > 0) + 0  # Easier and faster than ifelse(termite.plot>0,1,0)

detach(isoptera)
```




```r

mpd(termite.plot.PA, cophenetic(termite.taxophylo))
```

```
##   [1] 2.578 2.911 3.103 3.200 3.029 2.872 2.835 2.923 2.883 2.952 1.810
##  [12] 2.848 2.859 2.300 2.933 2.980 2.893 2.705 3.500 3.667 2.667 3.000
##  [23] 3.067 2.808 2.833 3.022 2.722 3.089 2.889 2.806 3.255 2.911 2.952
##  [34] 3.071 3.273 2.524    NA 3.536 3.194 3.167 3.090 2.933 3.267 3.500
##  [45] 2.933 2.606 2.745 3.205 3.133 2.964 2.879 2.848 2.648 2.978 3.111
##  [56] 2.733 3.091 3.088 3.533 3.056 2.952 3.000 3.071 3.111 2.764 2.762
##  [67] 2.861 2.806 3.000 3.571 2.964 3.133 3.278 2.667 3.321 2.758 2.619
##  [78] 3.067 2.694 3.222 2.758 2.667 3.000 2.691 3.200 2.643 2.867 3.143
##  [89] 2.667 2.697 3.048 2.864 3.067 3.291 3.150 3.077 3.041 3.043 3.200
## [100] 3.048 2.778 2.927 2.935 3.278 3.095 2.933 3.150 2.983 2.869 3.091
## [111] 3.000 3.055 2.810 2.961 3.077 3.067 2.912 3.045 3.176 2.857 3.008
## [122] 2.987 3.000 2.927 2.731 3.067 3.011 3.055 2.850 3.125 3.200 2.619
## [133] 2.758 2.910 2.758 2.737 3.025 2.953 3.013 2.846 2.859 2.843 2.955
## [144] 3.077 2.912 3.178 2.978 3.076
```

```r

ses.mpd <- ses.mpd(termite.plot.PA, cophenetic(termite.taxophylo))
```




```r

# Using functional traits (trophic group)



# Using ecological traits (like in Helmus et al. 2007a,b)

```


