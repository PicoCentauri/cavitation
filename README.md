# Water cavitation  results from the kinetic competition of bulk, surface and surface-defect nucleation events

[![arXiv](https://img.shields.io/badge/arXiv-2410.17626-B31B1B.svg)](https://arxiv.org/abs/2410.17626)

The repository includes codes and simulations parameter files used for the work "Water
cavitation  results from the kinetic competition of bulk, surface and surface-defect
nucleation events" available at https://arxiv.org/abs/2410.17626.

## Installation

In a fresh Python environment (virtualenv or conda), run the
following command to install the code and all dependencies:

```bash
pip install -r requirements.txt
```

## Usage

The [simulations](simulations) folder contains *GROMACS* structure, topology and parameter
files to run constant pressure rate simulations.

Additionally, the [analysis](analysis) folder contains Python functions to analyze
the simulations to extract the prefactor and further analysis based on the prefactor.
