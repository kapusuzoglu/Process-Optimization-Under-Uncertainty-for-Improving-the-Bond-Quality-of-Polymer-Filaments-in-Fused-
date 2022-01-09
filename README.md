## Cite Paper
Kapusuzoglu, Berkcan, et al. "Process Optimization Under Uncertainty for Improving the Bond Quality of Polymer Filaments in Fused Filament Fabrication." Journal of Manufacturing Science and Engineering 143.2 (2020).

Please, cite this repository using: 

    @article{berkcan2020Optimization,
        author = {Kapusuzoglu, Berkcan and Sato, Matthew and Mahadevan, Sankaran and Witherell, Paul},
        title = "{Process Optimization under Uncertainty for Improving the Bond Quality of Polymer Filaments in Fused Filament Fabrication}",
        journal = {Journal of Manufacturing Science and Engineering},
        pages = {1-46},
        year = {2020},
        month = {08},
        issn = "1087-1357",
        doi={10.1115/1.4048073},
        url = "[https://doi.org/10.1115/1.4048073](https://asmedigitalcollection.asme.org/manufacturingscience/article-abstract/143/2/021007/1086236/Process-Optimization-Under-Uncertainty-for?redirectedFrom=fulltext)"}

# Process-Design-Optimization-under-Uncertainty
This paper develops a computational framework to optimize the process parameters such that the bond quality between extruded polymer filaments is maximized in fused filament fabrication (FFF). A one-dimensional heat transfer analysis providing an estimate of the temperature profile of the filaments is coupled with a sintering neck growth model to assess the bond quality that occurs at the interfaces between adjacent filaments. Predicting the variability in the FFF process is essential for achieving proactive quality control of the manufactured part; however, the models used to predict the variability are affected by assumptions and approximations. This paper systematically quantifies the uncertainty in the bond quality model prediction due to various sources of uncertainty, both aleatory and epistemic, and includes the uncertainty in the process parameter optimization. Variance-based sensitivity analysis based on Sobol' indices is used to quantify the relative contributions of the different uncertainty sources to the uncertainty in the bond quality. A Gaussian process (GP) surrogate model is constructed to compute and include the model error within the optimization. Physical experiments are conducted to show that the proposed formulation for process parameter optimization under uncertainty results in high bond quality between adjoining filaments of the FFF product.

## Surrogate Model
Build a Gaussian process (GP) surrogate model to represent the model discrepancy using 4 inputs (temperature, speed, x and y coordinates of the bonds) and 1 output (model discrepancy)

## GSA
Global sensitivity analysis (GSA) is done for model parameters throughout the whole printing process.

## ProcessDesign
Perform Bayesian (surrogate) optimization
Process design optimization under uncertainty is performed using the posterior distributions of calibrated model parameters and the GP model representing the model discrepancy.

## GeometricAreaModel
Estimate the total void area of the manufactured parts.
