Unless otherwise stated the source code and scripts in this project are licensed under the GNU General Public License Version 3, see [LICENSE.md](https://github.com/mattiaghilardi/SIMMmultitroph_simulations/tree/main/LICENSE.md).

This project includes a modified version of the function `run_model()` of the `MixSIAR` R package:

- Original Authors: Brian Stock, Brice Semmens, Eric Ward, Andrew Parnell, Andrew Jackson, Donald Phillips
- Project Repository: https://github.com/brianstock/MixSIAR
- Original License: GPL-3 License
- Modified by: Mattia Ghilardi
- Modified on: 2026
- Modifications: 
  - function `run_model()` has been modified to run chains in parallel and renamed `run_model_parallel()`
  - `run_model_parallel()` has an additional argument 'jags.seed' for model reproducibility

All modified code resides in the [R/](https://github.com/mattiaghilardi/SIMMmultitroph_simulations/tree/main/R) directory.
