# extendedSurface
extendedSurface is an approach to fitting complex multi Ornstein-Uhlenbeck models to comparative datasets. It combines methods implemented in [SURFACE](https://www.otago.ac.nz/ecoevotago/code/surface.html) and [OUwie](https://www.jeremybeaulieu.org/r.html).

## Installation
```R
install.packages("devtools")
devtools::install_github("mongiardino/extendedSurface")
```

## Description
Multi Ornstein-Uhlenbeck (OUM) models are used to model the evolution of traits across a macroevolutionary adaptive landscape. Many approaches exist to fit these types of models, but they generally suffer from one of several shortcomings:
 1) They only work on ultrametric trees, even when paleontological data has been shown to improve accuracy of OU models (Ho & Ané 2014).
 2) They require *a priori* hypotheses on the number and location of regimes shifts (Beaulieau et al. 2012).
 3) They favor overly complex models (i.e., too many adaptive peaks; Khabbazian et al. 2016).
 
extendedSurface was developed in an attempt to solve these different issues. The method starts by employing the output of surface (Ingram & Mahler 2013), an approach that can employ trees with fossils and does not require any *a priori* input besides the trait data and tree topology. However, the method has been found to favor models that are too complex, a result that might stem from the fact that it assumes all regimes share a common rate of evolution and strength of attraction (Mongiardino Koch & Thompson 2020). This assumption is not only problematic from a model-fitting perspective, but is also likely biologically unrealistic. However, relaxing this assumption is not always feasible, as model fit with other approaches sush as OUwie (Beaulieu et al. 2012) might fail when attempting to optimize all of these parameters independently (Benson et al. 2017).

extendedSurface provides an easy interface between surface and OUwie, so that parameters of OUM models found with the former can be estimated with the latter. This allows to test the relative fit of models in which the regimes have different rates of evolution (OUMV), different strengths of attraction (OUMA), or both (OUMVA). Furthermore, it also allows optimizing the state at the root of the tree, which is otherwise assumed to be the optimal value of the ancestral regime (models OUMVZ, OUMAZ, OUMVAZ). All of this can be done both with the model preferred by the backwards phase of surface, as well as with simpler models obtained after independent regimes are merged (the extended phase of surface). Even though these simpler OUM models are suboptimal by surface standards, they have less number of parameters and model fitting with OUwie is more likely to be successful. Furthermore, simpler OUMV/OUMA/OUMVA models can fit the data drastically better that more complex OUM models. For example, Mongiardino Koch & Thompson (2020) found a 4-peak OUMVAZ model to strongly outperform the 7-peak OUM model favored by surface. These functions therefore provide a way to explore OUMV/OUMA/OUMVA(Z) models in non-ultrametric trees, without specifying *a prior* regime-shift hypotheses, and reducing problems of over-fitting.

![model_comparison](https://github.com/mongiardino/extendedSurface/blob/master/images/model_comparison.jpg)
**Fig. 1:** Comparison of the relative fit of multi Ornstein-Uhlenbeck models obtained for a comparative dataset of echinoid body size using Surface and extendedSurface (highlighted in yellow). Modified from Mongiardino Koch & Thompson (2020).

## Usage
The package comes with exemplary data that includes a time-calibrated phylogeny of echinoids, body size and measurement errors for all tips, as well as the outputs of running the forward and backwards phases of Surface on this data.
```R
data(echinoid_data)
surfaceExtended(bwd_surface = echinoid_data$bwd_surface,
                data = echinoid_data$size,
                tree = echinoid_data$tree,
                error = echinoid_data$error,
                models = 'OUMVAZ', limit = 4, plot = T,
                fwd_surface = echinoid_data$fwd_surface)
```

## Author
Nicolás Mongiardino Koch. Department of Earth & Planetary Sciences, Yale University.

**_Citation:_** Mongiardino Koch, N & Thompson, J. R. 2020. A Total-Evidence Dated Phylogeny of Echinoids and the Evolution of Body Size across Adaptive Landscape. bioRxiv, https://doi.org/10.1101/2020.02.13.947796.

## References
Beaulieu J.M., Jhuwueng D.‐C., Boettiger C., O'Meara B.C. 2012. Modeling stabilizing selection: expanding the Ornstein–Uhlenbeck model of adaptive evolution. Evolution, 66:2369–2383.

Benson R.B.J., Hunt G., Carrano M.T., Campione N. (2018), Cope's rule and the adaptive landscape of dinosaur body size evolution. Palaeontology, 61:13-48.

Ho L.S.T, Ané C. 2014. Intrinsic inference difficulties for trait evolution with Ornstein‐Uhlenbeck models. Methods in Ecology & Evolution, 5:1133–1146.

Ingram T., Mahler D.L. 2013. SURFACE: detecting convergent evolution from comparative data by fitting Ornstein‐Uhlenbeck models with stepwise Akaike Information Criterion. Methods in Ecology & Evolution, 4:416–425.

Khabbazian M., Kriebel R., Rohe K., Ané, C. 2016. Fast and accurate detection of evolutionary shifts in Ornstein‐Uhlenbeck models. Methods in Ecology & Evolution, 7:811–824.
