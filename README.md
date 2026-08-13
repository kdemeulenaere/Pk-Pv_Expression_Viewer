# '*P. knowlesi* Orthologue Expression Viewer' shiny App

`app.R` contains the R code to produce a Shiny App in which expression patterns of *P. knowlesi* A1-H.1 ([De Meulenaere et al., 2026](https://doi.org/10.64898/2026.02.09.704036)), *P. vivax* ([Zhu et al., 2016](https://doi.org/10.1038/srep20498)) and *P. knowlesi* ([Subudhi et al., 2020](https://www.nature.com/articles/s41467-020-16593-y)) genes can be visualised over the asexual blood stages: the ['*P. knowlesi* Orthologue Expression Viewer' web tool](https://interactive.itg.be/app/mal-pk-pv-expression-viewer).

The directory `/datasets` contains all datasets used in the script:
- `Pk_expr.txt`: *P. knowlesi* A1-H.1 ([De Meulenaere et al., 2026](https://doi.org/10.64898/2026.02.09.704036)) expression of the 5 sampled time points, and all interpolated time points.
- `Pv_expr.txt`: *P. vivax* natural isolate smru1 ([Zhu et al., 2016](https://doi.org/10.1038/srep20498)) expression of the 7 sampled time points, and all interpolated time points.
- `Pf_expr.txt`: *P. falciparum* strain II3 (derived from 3D7) ([Subudhi et al., 2020](https://www.nature.com/articles/s41467-020-16593-y)) expression of the 25 sampled time points, and all interpolated time points.
- `Pk_unfiltered.txt`: all *P. knowlesi* genes, and whether they passed filtering or not (sufficient transcriptional variation over the IDC, sufficient expression levels).
- `Pv_unfiltered.txt`: all *P. vivax* genes, and whether they passed filtering or not (sufficient transcriptional variation over the IDC, sufficient expression levels).
- `Pf_unfiltered.txt`: all *P. falciparum* genes, and whether they passed filtering or not (sufficient transcriptional variation over the IDC, sufficient expression levels).
- `ortho_PkPv.txt`: every row represents a *P. vivax* - *P. knowlesi* orthologue pair.
- `ortho_PkPf.txt`: every row represents a *P. falciparum* - *P. knowlesi* orthologue pair.
- `sim_PkPv.txt`: similarity statistics for the *P. vivax* - *P. knowlesi* orthologue pairs.
- `sim_PkPf.txt`: similarity statistics for the *P. falciparum* - *P. knowlesi* orthologue pairs.

The '*P. knowlesi* Orthologue Orthologue Expression Viewer' App can be accessed [here](https://interactive.itg.be/app/mal-pk-pv-expression-viewer).
- In the first tab of the App, a *P. vivax/P. knowlesi* gene is entered, the orthologue(s) are searched, and an expression plot for the input gene + orthologue(s) is given together with similarity statistics of the expression patterns.
- In the second tab of the App, a *P. falciparum/P. knowlesi* gene is entered, the orthologue(s) are searched, and an expression plot for the input gene + orthologue(s) is given together with similarity statistics of the expression patterns.
- In the third tab of the App, up to 4 *P. knowlesi/P. vivax/P. falciparum* genes are entered, and one expression plot is given for all input genes together.

**If you use the App, code, or *P. knowlesi* datasets, please cite: [De Meulenaere et al., 2026, doi: 10.64898/2026.02.09.704036](https://doi.org/10.64898/2026.02.09.704036).**

This is a project from [ITMmalaria](https://github.com/ITMmalaria).
