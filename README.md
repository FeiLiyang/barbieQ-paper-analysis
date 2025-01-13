# public_barcode_count

A [workflowr][] project.

[workflowr]: https://github.com/workflowr/workflowr

Liyang's analysis on the public datasets of barcode count. 

The repository is structured using workflowr. 

Each dataset corresponds to a subdirectory under ./data, ./analysis, ./docs, ./output.

*How to knit the Rmd's?*

The default `wflow::build()` is sabotaged by the unconventional set-up of subdirectories.

Use `rmarkdown::render()` to render the Rmd's into html's.

Example command in console where the root directory is always the project `"."`:

`rmarkdown::render("analysis/index.Rmd", output_dir = "docs")`

`rmarkdown::render("analysis/WuC/filtering_explore_alter.Rmd", output_dir = "docs/WuC")`

*How to load the datasets?*

We load datasets in an Rmd file, and the working directory is the directory of the current Rmd file.

This is diffrent from that in console.

(../../ moves two levels up; ../ moves one level up; ./ refers to current directory)

Example of load a data in `"analysis/WuC/filtering_explore_alter.Rmd"`:

`load("../../output/WuC/tagged_bq.rda")`

*How to dispatch custom html setting?*

The `scroll-plot` style and the `.color-tabs` style are saved in `"header.html"`.

These styles are configured by each `"_site.yml"` under each subdirectory of `"analysis/"`.

Example when dispatching in `"analysis/WuC/testscroll.Rmd"`:

 - flanking the r code trunk by ::: scroll-plot \n :::
 
 - adding {.color-tabs .tabset .tabset-fade .tabset-pills} following the title / subtitle
