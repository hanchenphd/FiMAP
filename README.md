# FiMAP (Fast Identity-by-Descent Mapping Test)
## Description
FiMAP is an R package for performing efficient identity-by-descent (IBD) mapping test using IBD segments to identify genomic regions associated with complex traits in large biobank-scale cohorts.
## Installing
FiMAP imports R packages <a href="https://CRAN.R-project.org/package=CompQuadForm">CompQuadForm</a>, <a href="https://CRAN.R-project.org/package=foreach">foreach</a>, <a href="https://CRAN.R-project.org/view=HighPerformanceComputing">parallel</a>, <a href="https://CRAN.R-project.org/package=Matrix">Matrix</a>, <a href="https://CRAN.R-project.org/package=GMMAT">GMMAT</a>. In addition, FiMAP requires <a href="https://CRAN.R-project.org/package=doMC">doMC</a> to run parallel computing. These dependencies should be installed before installing FiMAP.

For optimal computational performance, it is recommended to use an R version configured with the Intel Math Kernel Library (or other fast BLAS/LAPACK libraries). See the <a href="https://software.intel.com/en-us/articles/using-intel-mkl-with-r">instructions</a> on building R with Intel MKL.

To install FiMAP from GitHub, please use
```
devtools::install_github("hanchenphd/FiMAP", ref = "main")
```
## Version
The current version is 0.1.0 (July 1, 2021).
## License
This software is licensed under GPL-3.
## Contact
Please refer to the R help document of FiMAP for specific questions about each function. For comments, suggestions, bug reports and questions, please contact Han Chen (Han.Chen.2 AT uth.tmc.edu). For bug reports, please include an example to reproduce the problem without having to access your confidential data.
## References
<p>If you use FiMAP, please cite
<li>Chen H, Naseri A, Zhi D. (2021) FiMAP: A Fast Identity-by-Descent Mapping
Test for Biobank-scale Cohorts. <em>medRxiv</em> 2021.06.30.21259773
[<b>Preprint</b>]. July 5, 2021. DOI: <a href="https://doi.org/10.1101/2021.06.30.21259773">10.1101/2021.06.30.21259773</a>.</li></p>



