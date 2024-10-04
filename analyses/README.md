## Analyses

This subdirectory contains the code to run the analyses. For Quarto (.qmd) files, the html rendered version is also available. They are numbered according to the order in which they must be rum but some analyses are independent. The corresponding results are stored in the subfolder `outputs/xxx/` where `xxx/` is the analysis subfolder of interest.

-   `01_clean_data/` contains the code to format and clean network data.
-   `02_niche_example_figure/` contains the code to reproduce the figure exemplifying niche measures from reciprocal scaling.
-   `03_simulation/` contains the code to simulate interaction networks and analyze the results to evaluate reciprocal scaling (main text and appendices). 
-   `04_data_analysis/` contains the code to preform reciprocal scaling on the Peru1 network from the ANDEAN frugivory dataset (Dehling et al., 2021).

## References

Dehling, D. M., Bender, I. M. A., Blendinger, P. G., Muñoz, M. C., Quitián, M., Saavedra, F., Santillán, V., Böhning-Gaese, K., Neuschulz, E.-L., & Schleuning, M. (2021). ANDEAN frugivory: Data on plant-bird interactions and functional traits of plant and bird species from montane forests along the Andes (Version 2, p. 187397 bytes) [dataset]. Dryad. <https://doi.org/10.5061/DRYAD.WM37PVMN5>
