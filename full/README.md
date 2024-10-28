# Quantum Key Distribution System Modeling

**Author**: Ibrahim Almosallam  
**Contact**: [ibrahim@almosallam.org](mailto:ibrahim@almosallam.org)

## Overview

This repository contains the code and resources accompanying the research paper, *"Comprehensive Modeling of Quantum Key Distribution Systems"*. The work proposes a Bayesian inference-based approach to enhance the estimation accuracy of the tagged photon fraction, \(\Delta\), within Quantum Key Distribution (QKD) systems. This refined estimation aims to improve secure key rates, lower system costs, and extend operational distances.

## Repository Structure

- **`main.pdf`**: The manuscript detailing the theory, methodology, and results.
- **`session.m`**: A main script to run and experiment with the QKD model, allowing for custom parameter adjustments and insights.
- **`experiments/`**: Folder containing experiment scripts used to generate results and figures for the manuscript.

### Experiment Scripts

Each experiment script in the `experiments/` folder is dedicated to a specific part of the analysis:

- **`experiment1_validate_iid.m`**: Validates the i.i.d. assumption in the detection model.
- **`experiment2_validate_hmm.m`**: Tests the hidden Markov model for dependency in detection events.
- **`experiment3_deltas.m`**: Studies the estimation accuracy for different tagged photon fractions, \(\Delta\).
- **`experiment4_rates.m`**: Compares secure key rates across various parameter configurations.
- **`experiment5_eve_iid.m`**: Models and analyzes the behavior of an eavesdropper under i.i.d. assumptions.
- **`experiment6_eve_hmm.m`**: Simulates eavesdropping effects using the HMM-based dependency model.
- **`experiment7_full_iid.m`**: Full system simulation with i.i.d. assumptions for validation against idealized conditions.
- **`experiment8_full_hmm.m`**: Full system simulation incorporating HMM-based modeling of dependencies.

## Requirements

This code is developed in MATLAB (version 2023 or later). Ensure the following packages or toolboxes are available:

- **Statistics and Machine Learning Toolbox** (for Bayesian and statistical inference)
- **Signal Processing Toolbox** (for advanced modeling and signal handling)

## Getting Started

1. **Setup**: Clone this repository and ensure all MATLAB dependencies are available.
2. **Running `session.m`**: Use `session.m` as the main entry point to explore the core QKD model, allowing for real-time adjustments of parameters.
3. **Experiments**: Run each experiment script in the `experiments/` folder individually. Results generated from these scripts will replicate the figures and findings in the manuscript.

### Running Experiments

To run a specific experiment, navigate to the `experiments/` folder in MATLAB and execute the desired script. Each experiment is self-contained and can be run independently. For example, to validate the i.i.d. assumption:

```matlab
cd experiments
run('experiment1_validate_iid.m')

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

## Contact

For questions or further details, please contact [ibrahim@almosallam.org](mailto:ibrahim@almosallam.org).

