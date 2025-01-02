This repository provides materials for an article, "A renewal-equation approach to estimating $R_t$ and infectious disease case counts in the presence of reporting delays" by Sumali Bajaj, Robin Thompson and Ben Lambert.

We use [renv](https://rstudio.github.io/renv/articles/renv.html) to manage the environment, and after cloning this article and opening the R project file, the R environment can be reproduced by running the following from the RStudio console:

```r
renv::restore()
```

We use the [targets R package](https://books.ropensci.org/targets/) to create reproducible data analysis pipelines. To rerun all analyses, type the following in the Rstudio console:


```r
targets::tar_make()
```

To generate only a specific target, consult the `_targets.R` file to find a particular target name (e.g. `file_plot_projections_waning_over_png` which generates the `projections_waning.png` file) and run:

```r
targets::tar_make(file_plot_projections_waning_over_png)
```
