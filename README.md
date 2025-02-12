# Sea Ice Data Assimilation with DART

This repository contains scripts and workflows for assimilating sea ice data using the Data Assimilation Research Testbed (DART). The focus of this project is on comparing the assimilation of broadband albedo (SIAL) with traditional sea ice concentration (SIC) and sea ice thickness (SIT) within a bounded data assimilation (DA) framework. The implementation includes adaptive inflation techniques to improve state estimates.

## Project Overview

The assimilation framework leverages relationships among key Icepack variables (SIC, SIT, SIAL) to investigate their impact on sea ice prediction. The current methodology is based on ensemble data assimilation techniques using the Data Assimilation Research Testbed (DART) provided by the National Center for Atmospheric Research (NCAR).

## Key Components

- **Synthetic SIC and Albedo Data Generation**: Scripts for creating synthetic observations of SIC and albedo.
- **Observation Processing**: Workflow for handling and formatting observation sequences for DART.
- **DART Configuration**: Setup files for running bounded DA schemes and adaptive inflation.
- **Visualization and Analysis**: Scripts for plotting SIC-based and albedo-based predictions, comparing ensemble means, and evaluating assimilation performance.

## Dependencies

- Python (NumPy, Matplotlib, xarray, SciPy, etc.)
- DART (Data Assimilation Research Testbed)
- Computing access (for running large-scale ensemble experiments)

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.com/jrotondo-uw/cice-scm-albedo-da.git
   cd cice-scm-albedo-da.git
   ```
2. Install necessary Python dependencies:
   ```bash
   pip install -r requirements.txt
   ```
3. Set up DART and configure paths accordingly.

## Running Assimilation Experiments

1. **Prepare observations**:
   ```bash
   python preprocess_observations.py
   ```
2. **Run DART assimilation**:
   Follow the standard DART workflow for setting up and executing assimilation cycles.
3. **Visualize results**:
   ```bash
   python plot_results.py
   ```

## Acknowledgments

This work is conducted as part of research at the University of Washington, leveraging NCAR’s Derecho supercomputer.

## License

This project is licensed under the MIT License. See the `LICENSE` file for details.

## Contact

For questions or contributions, please open an issue or reach out to Joey Rotondo, jrotondo@uw.edu.
