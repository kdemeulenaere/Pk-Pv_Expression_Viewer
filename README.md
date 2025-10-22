# '*P. knowlesi* - *P. vivax* Expression Viewer' shiny App

`app.R` contains the R code to produce a Shiny App in which expression patterns of *P. knowlesi* A1-H.1 (De Meulenaere et al., 2025) and *P. vivax* ([Zhu et al., 2016](https://doi.org/10.1038/srep20498)) genes can be visualised over the asexual blood stages: the '*P. knowlesi* - *P. vivax* Expression Viewer' (available [here](https://interactive.itg.be/app/mal-pk-pv-expression-viewer)).

The directory `/datasets` contains all datasets used in the script:
- `Pk_expr.txt`: *P. knowlesi* A1-H.1 (De Meulenaere et al., 2025) expression of the 5 sampled time points, and all interpolated time points.
- `Pv_expr.txt`: *P. vivax* smru1 ([Zhu et al., 2016](https://doi.org/10.1038/srep20498)) expression of the 7 sampled time points, and all interpolated time points.
- `Pk_unfiltered.txt`: all *P. knowlesi* genes, and whether they passed filtering or not (sufficient transcriptional variation over the IDC, sufficient expression levels).
- `Pv_unfiltered.txt`: all *P. vivax* genes, and whether they passed filtering or not (sufficient transcriptional variation over the IDC, sufficient expression levels).
- `ortho.txt`: every row represents a *P. vivax* - *P. knowlesi* orthologue pair.
- `sim.txt`: similarity statistics for the *P. vivax* - *P. knowlesi* orthologue pairs.

The '*P. knowlesi* - *P. vivax* Expression Viewer' App can be accessed [here](https://interactive.itg.be/app/mal-pk-pv-expression-viewer).
- In the first tab of the App, a *P. vivax/P. knowlesi* gene is entered, the ortholog(s) are searched, and an expression plot for the input gene + ortholog(s) is given together with similarity statistics of the expression patterns.
- In the second tab of the App, up to 4 *P. knowlesi/P. vivax* genes are entered, and one expression plot is given for all input genes together.

**If you use the App, code, or *P. knowlesi* datasets, please cite: De Meulenaere et al., 2025, XXX. doi: XXX.** (yet unpublished)

This is a project from [ITMmalaria](https://github.com/ITMmalaria).
