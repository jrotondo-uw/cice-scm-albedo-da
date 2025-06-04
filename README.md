# Sea Ice Data Assimilation with DART

This repository contains scripts and workflows for assimilating sea ice data using the Data Assimilation Research Testbed (DART). The focus of this project is on comparing the assimilation of narrow-band sea ice albedo (SIAL) with traditional sea ice concentration (SIC) and sea ice thickness (SIT) within a bounded data assimilation (DA) framework. The implementation includes adaptive inflation techniques to improve state estimates.

## Project Overview

The assimilation framework leverages relationships among key Icepack variables (SIC, SIT, SIAL) to investigate their impact on sea ice prediction. The current methodology is based on ensemble data assimilation techniques using the Data Assimilation Research Testbed (DART) provided by the National Center for Atmospheric Research (NCAR).

## Key Components

- **Synthetic SIC and Albedo Data Generation**: Scripts for creating synthetic observations of SIC and albedo.
- **Observation Processing**: Workflow for formatting observation sequences compatible with DART.
- **DART Configuration**: Setup for bounded DA and adaptive inflation techniques.
- **Visualization and Analysis**: Scripts and notebooks for plotting SIC-based and albedo-based predictions, comparing ensemble means, and evaluating assimilation performance.

## Repository Structure

- **`scripts/`**: Python and shell scripts for generating synthetic observations, preparing input files, running DA cycles, and plotting results.
- **`obs_seq/`**: Pre-generated observation sequence files used in assimilation experiments.
- **`Icepack_mods/`**: Modifications to Icepack code to allow assimilation of narrow-band sea ice albedo within the single-column model (SCM) framework.
- **`Dart_mods/`**: Changes to DART code (especially `model_mod.f90`) to support albedo assimilation in the cice-scm configuration. The pointer variable `e` in `model_mod.f90` must be manually set to `1`, `2`, etc. to specify which model observation type is assimilated before rebuilding DART.
- **`data/templates/`**: Contains the DART namelists (`input.nml`, `model_nml`, `obs_kind_nml`, etc.) used for configuration of the filter, inflation, and model integration.
- **`figures/`**: Output plots generated during development and evaluation.
- **`Jupyter_Notebooks/paper_figures.ipynb`**: Contains all code to reproduce the final figures used in the manuscript.
- **`data/`**: Contains only the final post-processed datasets used to generate manuscript figures. The full raw synthetic observations and model output from free runs are not included due to storage limitations.

## Data Availability

The final processed data used to produce the figures in the manuscript is archived on Zenodo:  
➡️ [https://doi.org/10.5281/zenodo.15571204](https://doi.org/10.5281/zenodo.15571204)

Raw synthetic observations and full free-run model outputs are not included due to size constraints, but the workflow to regenerate them is available in `scripts/`.

## Dependencies

- **DART cice-scm**:  
  [https://github.com/NCAR/DART/tree/cice-scm/models/cice-scm](https://github.com/NCAR/DART/tree/cice-scm/models/cice-scm)

- **Icepack v1.4.0**:  
  [https://github.com/CICE-Consortium/Icepack/tree/Icepack1.4.0](https://github.com/CICE-Consortium/Icepack/tree/Icepack1.4.0)

Use the modified files in `Icepack_mods/` and `Dart_mods/` to enable assimilation of albedo observations. Be sure to update the hardcoded pointer `e = ...` in `model_mod.f90` before building DART to ensure you’re assimilating the intended model variable.

## Getting Started

1. Clone this repository:
   ```bash
   git clone https://github.com/jrotondo-uw/cice-scm-albedo-da.git
   cd cice-scm-albedo-da
   ```

2. Follow the workflows in `scripts/` to regenerate observations and run assimilation cycles.

3. Use `Jupyter_Notebooks/paper_figures.ipynb` to recreate the final figures and compare the performance of SIC-, SIT-, and SIAL-based assimilation schemes.

---

For questions or collaboration, feel free to reach out via the GitHub issues page.
