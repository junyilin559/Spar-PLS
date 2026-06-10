# Sparsification Subsampling for Partial Least Squares Regression
This repository includes the implementation of our work **"Sparsification Subsampling for Partial Least Squares Regression"**.

## Introduction
A brief introduction about the folders and files:
* `src/`: the source file for different PLS methods;
    * `computing.cpp`: Calculates probabilities and sampling functions.
    * `Spar-PLS.R`: Implementation of Spar-PLS.
    * `RandomPLS.R`: Implementation of row-sampling PLS.
    * `RandomPLS.R`: Implementation of row-sampling PLS.
    * `IFPLS.R`: Implementation of IFPLS.
    * `MCDPLS.R`: Implementation of MCDPLS.
    * `Utils.R`: Implementation of SGDPLS and CIPLS, by referencing `SGDPLS.py` and `CIPLS.py`.
* `simu/`: simulation scripts;
    * `BMSE.R`: comparison of FULL, UNIF, LEV, IF, MCD and SPAR w.r.t. BMSE;
    * `QMSE.R`: comparison of FULL, UNIF, LEV, IF, MCD and SPAR w.r.t. QMSE;
    * `PMSE.R`: comparison of FULL, UNIF, LEV, IF, MCD and SPAR w.r.t. PMSE;
    * `time.R`: Time comparison of FULL, UNIF, LEV, IF, MCD and SPAR;
    * `SmallSample.R`: comparison of FULL, UNIF, LEV, and SPAR w.r.t. BMSE, QMSE, PMSE on small sample settings;
* `real_data/`: real data analysis script;
    * `WEC.R`: comparison of FULL, UNIF, LEV, IF, MCD and SPAR w.r.t. PMSE;
    * `WEC_time.R`: Time comparison of FULL, UNIF, LEV, IF, MCD and SPAR;
    * `WEC_full.R`: comparison of FULL, UNIF, LEV, IF, MCD, CIPLS, SGDPLS, SPAR w.r.t. PMSE and time;
    * `LEV_compare.R`: comparison of FULL, UNIF, LEV, LEV*, SPAR(COV), SPAR(LEV), SPAR(LEV*) w.r.t. PMSE and time;

## Reproducibility
For simulation studies in Section 4 and the Supplementary Material,
* you can run `BMSE.R` to reproduce the results in Figure 3;
* you can run `QMSE.R` to reproduce the results in Figure 4;
* you can run `PMSE.R` to reproduce the results in Figure 5;
* you can run `time.R` to reproduce the results in Table 2;
* you can run `SmallSample.R` to reproduce the results in Section C.2;

For the real data example in Section 6,
* Download the [Wave Energy Converters dataset](https://archive.ics.uci.edu/dataset/494/wave+energy+converters) and store them in the `real_data/data/` path;

* you can run `WEC.R` to reproduce the results in Figure 7;
* you can run `WEC_time.R` to reproduce the results in Table 3;
* you can run `WEC_full.R` to reproduce the results in Section C.1;
* you can run `LEV_compare.R` to reproduce the results in Section C.3;

## Citation

If you found this repository useful, please cite the following.

```bibtex
@article{lin2026sparpls,
      author = {Junyi Lin and Mengyu Li and Cheng Meng and Yongdao Zhou},
      title = {Sparsification Subsampling for Partial Least Squares Regression},
      journal = {Journal of Computational and Graphical Statistics},
      volume = {0},
      number = {ja},
      pages = {1--30},
      year = {2026},
      publisher = {Taylor \& Francis},
      doi = {10.1080/10618600.2026.2686429},
      URL = {https://www.tandfonline.com/doi/abs/10.1080/10618600.2026.2686429},
}
```
