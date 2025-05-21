# Cn-chart
Multivariate Control Chart using Copulas 

This repository is related to the proposed methodology $c_n$-Chart, which proposes a new control chart based on the combination of marginal densities combined with a copula density for the tolerance region estimation.

This study proposes a purely data-driven alternative for a traditional control tolerance region that also generalizes the multidimensional setting ($d \geq 2$). The developed methodology adopted the joint distribution estimated by combining the process marginal distributions linked with a copula function. This approach allows visual change detection, e.g., in the process mean vector in multivariate contexts, although the application of the $c_n$ multivariate control chart (or simply, $c_n$-Chart) is not restricted to only two-dimensional problems, its result will always be a projection on a two-dimensional (2D) graph.


The R script associated with the SUPPLEMENTARY MATERIAL (Regression Extensions) includes two case studies: (i) water quality control, measuring the relationship between the potential of hydrogen (pH) levels and the phosphate (PO4) concentration (Subsection \ref{sec:ap_biv}); and (ii) brass-steel thermostats monitoring, considering three response variables: deflection, curvature, and resistivity. The `script_copula.R` is related to the Water Quality Control Case, and the `script_copula2.R` to the Brass-Steel Thermostat Monitoring study.
