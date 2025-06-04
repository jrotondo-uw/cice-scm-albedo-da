# Sea Ice Data Assimilation with DART

This repository contains scripts and workflows for assimilating sea ice data using the Data Assimilation Research Testbed (DART). The focus of this project is on comparing the assimilation of narrow-band sea ice albedo (SIAL) with traditional sea ice concentration (SIC) and sea ice thickness (SIT) within a bounded data assimilation (DA) framework. The implementation includes adaptive inflation techniques to improve state estimates.

## Project Overview

The assimilation framework leverages relationships among key Icepack variables (SIC, SIT, SIAL) to investigate their impact on sea ice prediction. The current methodology is based on ensemble data assimilation techniques using the Data Assimilation Research Testbed (DART) provided by the National Center for Atmospheric Research (NCAR).

## Key Components

- **Synthetic SIC and Albedo Data Generation**: Scripts for creating synthetic observations of SIC and albedo.
- **Observation Processing**: Workflow for handling and formatting observation sequences for DART.
- **DART Configuration**: Setup files for running bounded DA schemes and adaptive inflation.
- **Visualization and Analysis**: Scripts for plotting SIC-based and albedo-based predictions, comparing ensemble means, and evaluating assimilation performance.

## Repository Structure

- **`scripts/`**: Shell and Python scripts for preprocessing, running the model, and analyzing results.
- **`namelists/`**: Configuration files for model and assimilation parameters.
- **`src/`**: Source code modifications and custom modules for DART and CICE-SCM.
- **`results/`**: Output from assimilation experiments and diagnostic figures.
- **`data/`**: Includes only the final processed datasets used to create the figures in the manuscript.

## Data Availability

The final processed data used to produce the figures in the manuscript is archived on Zenodo:  
➡️ [https://doi.org/10.5281/zenodo.15571204](https://doi.org/10.5281/zenodo.15571204)  

Raw synthetic observations and model free runs are not included due to size constraints, but the entire workflow for regenerating these datasets is available in the `scripts/` directory.

## Dependencies

- **DART cice-scm**: [https://github.com/NCAR/DART/tree/cice-scm/models/cice-scm](https://github.com/NCAR/DART/tree/cice-scm/models/cice-scm)
- **Icepack v1.4.0**: [https://github.com/CICE-Consortium/Icepack/tree/Icepack1.4.0](https://github.com/CICE-Consortium/Icepack/tree/Icepack1.4.0)

Please ensure these are cloned and built as described in the DART and Icepack documentation.

## Getting Started

1. Clone this repository:
   ```bash
   git clone https://github.com/jrotondo-uw/cice-scm-albedo-da.git
   cd cice-scm-albedo-da
   ```

2. Review and run the scripts in `scripts/` to generate or assimilate sea ice observations.

3. Use the visualization scripts to replicate figure outputs and diagnostics.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Acknowledgments

This work uses the DART data assimilation framework developed by NCAR and Icepack developed by the CICE Consortium. We thank the developers and community for their support and resources.
