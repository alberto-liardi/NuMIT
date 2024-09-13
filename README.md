NuMIT: Null Models for Information Theory
=========================================

This repository provides an implementation for creating null models within the context of information theory.
These models enable a mathematically robust normalisation procedure, allowing for comparisons of information-theoretic 
quantities across different datasets.


Download and Installation
-------------------------

To download the repository, either download the `.zip` file or use `git clone`. 
Before using the package, run the `startup.m` script to initialise the repository. 

* If MATLAB is opened from the repository's root directory, `startup.m` will be automatically executed.
* Otherwise, run `startup` from the MATLAB command window.

Do not add the entire NuMIT directory tree to the path, as it may cause errors. 
The `startup` script will automatically set the necessary paths.

**Note:** This package contains a copy of [MVGC2](https://github.com/SacklerCentre/MVGC2) and 
[partial-info-decomp](https://github.com/robince/partial-info-decomp) libraries. 
Full credit for both packages goes to the original authors, please refer to the original 
repositories and related papers [[1](https://www.sciencedirect.com/science/article/pii/S0165027013003701?via%3Dihub)] 
[[2](https://www.mdpi.com/1099-4300/19/7/318)] for more details.


Usage
-----

The package offers a null model-based normalisation procedure for Partial Information Decomposition (PID), 
Integrated Information Decomposition (PhiID), and Integrated Information measures ($\Phi_{WMS}$ and $\Phi_R$) 
for three different system types: Discrete, Gaussian, and Vector Autoregression (VAR). 

For a given a set of measures (e.g. PID or PhiID atoms) calculated from real data (and using one of the metrics and systems above), 
NuMIT normalisation generates a null distribution of systems with the same total mutual information as 
the original system, but otherwise random. The original measures are then compared to the null distribution 
by taking their quantiles, resulting in the NuMIT-normalised measures. 

For a more detailed description, please refer to the main article.

**Note**: At the moment, only these combinations are available:
* Discrete: PID
* Gaussian: PID
* VAR: PID, PhiID, Phis
  
### Discrete
PID atoms are compared to those randomly generated from systems of the form $T=f(X,Y)+\epsilon$, where $f$ represents a 
2-bit logic gate, and $\epsilon$ is a noise term which flips the value of $f(X,Y)$ with probability $p_{\epsilon}$.

The original PID atoms can be fed to `NuMIT_PID` by specifying the Discrete model (`model('name',"Discrete",...)`), which 
will generate the null distribution through `Null_model_Discrete`. 

Note that these functions can be called separately, if one wants to examine the properties of the null distribution. 
A complete usage example is provided in `Discrete_test.m`. 

**Note**: Only MMI definition is available for the discrete PID.

### Gaussian
The Gaussian model considered is of the form $T=AS+\epsilon$, where $S$ represents jointly Gaussian-distributed sources, $A$ a 
matrix of the coefficients, and $\epsilon$ a white-noise term. 

Similarly to the discrete case, original PID atoms can be given as input to `NuMIT_PID`, specifying `model('name',"Gauss",...)` 
in the input. These are then compared to a null distribution (constructed in `Null_model_Gauss`) and their quantiles are returned.
See `Gaussian_test.m` for more details.

**Note**: MMI, DEP, and CCS definitions are available for the Gaussian PID by setting e.g. `model(...,'red_fun',"MMI",...)` in the model 
specification.

### Vector Autoregression
In addition to providing the NuMIT normalisation, the package also offers a pipeline to fit and calculate PID, PhiID, and Phis 
information quantities for multivariate time series data using a VAR model (see `VAR_analysis`). This procedure randomly samples
`N` channels and `E` epochs from the total number of channels of the data, fits a VAR model on their time series, and computes 
PID/PhiID/Phis measures by separating the channels into two random groups. This procedure is repeated `N_run` times, each time 
sampling different channels and epochs. 

Once the _raw_ information quantities are computed, the NuMIT normalisation can be achieved by either calling `NuMIT_PID`, 
`NuMIT_PhiID`, or `NuMIT_Phis`, depending on the measures calculated, and specifying `model('name',"VAR",...)` as input. 
See `VAR_PID_example.m` and `VAR_PhiID_example.m` for a more detailed working example.  

**Note**: MMI, DEP, and CCS definitions are available for the VAR PID by setting e.g. `model(...,'red_fun',"MMI",...)` in the model 
specification. For VAR PhiID only MMI is currently implemented.



Licence
-------

This software is distributed under the modified 3-clause BSD Licence.


Further references
---------------

* P. Williams and R. Beer (2010). Nonnegative decomposition of multivariate
  information. [arXiv:1004.2515](https://arxiv.org/abs/1004.2515)
  
* P. Mediano\*, F. Rosas\*, _et al._ (2019). Beyond integrated information: A
  taxonomy of information dynamics phenomena.
  [arXiv:1909.02297](https://arxiv.org/abs/1909.02297)

* A: Barrett, (2015). Exploration of synergistic and redundant information sharing
  in static and dynamical Gaussian systems.
  [doi:10.1103/PhysRevE.91.052802](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.91.052802)

* R. Ince, (2017). Measuring multivariate redundant information with pointwise
  common change in surprisal. [doi:10.3390/e19070318](https://doi.org/10.3390/e19070318)

* R. G. James, _et al._ (2018). Unique information via dependency constraints.
  [doi:10.1088/1751-8121/aaed53](https://iopscience.iop.org/article/10.1088/1751-8121/aaed53/pdf)


\(C\) Alberto Liardi, 2024
