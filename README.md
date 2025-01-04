# README

This repository contains the MATLAB implementation described in the paper [Overcoming Intensity Limits for Long-Distance Quantum Key Distribution](https://arxiv.org/abs/2412.20265). The code is organized into two main sections:

1. **`eve/`** – Demonstrates inference when only Eve’s parameters $\theta_E$ (i.e., $d_{AE}$, $p_{EB}$, $k$, $\Delta$) are random.  
2. **`full/`** – Extends the approach to a fully flexible Bayesian model, allowing *any* QKD parameter (e.g., fiber attenuation $\alpha$, detector efficiencies $p_c$, misalignment $p_e$, etc.) to be random.

Each folder contains:
- A demo script (`demo.m`) showing how to run a sample inference procedure (see each folder’s `README.md` for details).  
- An `experiments/` folder with scripts to reproduce the results presented in the paper.

## Requirements

- **MATLAB** 2024b (Version: 24.2) or later  
- **Optimization Toolbox** (Version: 24.2)  
- **Statistics and Machine Learning Toolbox** (Version: 24.2)

---

## License

This code is released under the MIT License. See the [LICENSE](LICENSE) file for details.

---

## Contact

For inquiries, feedback, or collaboration requests, please email  [Ibrahim Almosallam](mailto:ibrahim@almosallam.org)  or open an issue in this repository.