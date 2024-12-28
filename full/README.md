# Fully Bayesian

The script (`demo.m`) demonstrates a **fully flexible** Bayesian inference approach where **any** QKD-related parameter (including the parameters in $\theta_A$ and $\theta_B$) can be treated as random variables. You can specify which variables are fixed vs. random, then jointly perform inference on all random variables. Below is an outline of the major steps the script performs, followed by a detailed listing of the variables used in each step.

---

## Script Outline

1. **Environment Setup**  
   Initializes the MATLAB environment (path, clearing workspace, setting a seed for reproducibility) and configures file-handling switches (load or save data, save plots, etc.).

2. **Variable Sets**  
   Declares a comprehensive list of all possible QKD variables (e.g., $\lambda$, $\alpha$, $d_{AB}$, $p_{a_0}, p_{a_1}, \ldots$, $d_{AE}$, $\Delta$, etc.), then partitions them into *fixed* vs. *random* sets.

3. **System Parameters**  
   Defines baseline values for Alice ($\theta_A$), Bob ($\theta_B$), and Eve ($\theta_E$)—just like the partial approach, but you may let any subset become random by including them in the “random” set.

4. **Priors and Noise**  
   Configures how random variables are initialized, including potential noise added to them. Also sets up the shape parameters or bounds used later for MAP and MCMC.

5. **Simulation**  
   Either simulates a fresh QKD session (generating click counts) using the chosen baseline parameters, or loads previously saved data if `loadData = true`.

6. **MAP Estimation**  
   Performs a maximum a posteriori optimization for whichever parameters are designated as random, comparing the results to the known ground-truth values.

7. **MCMC Sampling**  
   From the MAP point, samples each random variable’s posterior distribution via the selected sampler (e.g., `'srss'`), storing histograms and progress information if requested.

8. **Key Rate Estimation**  
   Evaluates error/gain probabilities and obtains secure key rates for each sample. Compares these results to the “true” asymptotic rates from the ground-truth values.

9. **Saving Data and Displaying Results**  
   Stores simulation outputs, MAP estimates, and posterior samples in a `.mat` file (if specified).

10. **Displaying Results**  
   Plots histograms of the inferred parameters and the resulting key rates, optionally saving figures to PDF.

---

## Detailed Variables by Section

### 1. Environment Setup

- **`seed`**: An integer controlling `rng(seed)` for reproduciblity.

- **`file_path`**: If non-empty, the script will save new results or load previously stored data (`.mat`).

- **`loadData`** (boolean)  
  - `true`: Skips simulation and inference, loading data directly from `file_path`.  
  - `false`: Performs the full simulation/inference pipeline and saves the resulting data to `file_path` (if specified).

- **`savePlots`** (boolean)  
  - `true`: Outputs figures (e.g., parameter histograms, key rates) to PDF.  
  - `false`: Displays them without saving.

---

### 2. Variable Sets

- **`variables`**: A comprehensive cell array listing all potential QKD parameters (intensities, attenuation, misalignment, after-pulsing, etc.).
- **`varR`**: The subset of `variables` designated as *random* (e.g., $d_{AE}$, $p_{EB}$, $k$, $\Delta$ etc).
- **`varF`**: The complementary subset (fixed). Computed by `setdiff(variables, varR)` to ensure disjoint sets.
- **`noise`**: Control optional variance injection in random variables (e.g., to reflect measurement noisee).
- **`sigma`**: The standard deviation of the prior distributions for the random variables (a `sigma` multiple of the expected value).

---

### 3. System Parameters

- **Bob’s ($\theta_B$)**  
  - **`pa0`, `pa1`**: After-pulsing probabilities for Bob’s two detectors (with slight ±10% variations).  
  - **`pc0`, `pc1`**: Detector efficiencies (again ±10% variations).  
  - **`pd0`, `pd1`**: Dark count probabilities (±10%).  
  - **`pe`**: Misalignment probability.  
  - These remain fixed and are grouped as $\theta_B = \\{pa0, pa1, pc0, pc1, pd0, pd1, pe\\}$.

- **Alice’s ($\theta_A$)**  
  - **`lambdas`**: A vector of intensity levels, derived from a heuristic function (e.g., `get_lambdas`) that spans $\lambda_{\min}$ to a `limit`.  
  - **`alpha`**: The fiber attenuation coefficient (dB/km).  
  - **`dAB`**: Distance between Alice and Bob (km).  
  - These are grouped as $\theta_A = \\{\Lambda, \alpha, d_{AB}\\}$.

- **Eve’s ($\theta_E$)**  
  - **`dAE`**: Eve’s distance from Alice (km).  
  - **`pEB`**: Eve’s channel efficiency, computed via a heuristic (e.g., `get_pEB`) so Eve remains inconspicuous for a chosen $k$.  
  - **`k`**: Number of photons intercepted per pulse.  
  - **`Delta`**: Fraction of pulses intercepted by Eve.  
  - These four are the **random** parameters we infer: $\theta_E = \\{d_{AE}, p_{EB}, k, \Delta\\}$.

- **`N`, `Nl`, `limit`**  
  Pulse count, number of intensities, and maximum intensity cap.

- **`algorithm`, `maxIters`**  
  MAP optimization method (e.g., `'quasi-newton'`), plus iteration limit.

- **MCMC**  
  - `Ns`, `Nb`, `method` (e.g. `'slice'`), `display`, `chunkSize`, etc.

---

### 4. Priors and Noise

- **`processVars(...)`**: Allocates each variable to fixed vs. random sets (taking into consideration `sigma`).

- **`sampleParameters(...)`**: Allows adding `noise` to the baseline, shaping initial random variables used in simulation.

- **`processParams(...)`**: Prepares the prior hyper-parameters $thetaP$
  - **$\alpha$, $\beta$**: Shape parameters for Beta/Gamma priors on each of Eve’s parameters.  
  - **`ub`, `lb`**: Vectors of upper/lower bounds for $d_{AE}$, $p_{EB}$, $k$, $\Delta$.  
  - Grouped as $\theta_P = \\{\alpha, \beta, \text{ub}, \text{lb}\\}$.

---

### 5. Simulation

- **`simulate(N, thetaA, thetaB, thetaE)`**: Produces synthetic click counts $\mathbf{C}$, given the total pulses `N`, plus the system parameters for Alice, Bob, and Eve (the “true” $\theta_E$ values).

---

### 6. MAP Estimation

- **`MAP(...)`**: Uses an optimization routine (e.g., `'quasi-newton'`) up to `maxIters` iterations to find the $\theta_R$ maximizing the posterior given $\mathbf{C}$.  
- **`thetaE_MAP`**: The resulting MAP estimate, which is then compared to the ground-truth $\theta_R$.

---

### 7. MCMC Sampling

- **`sample(...)`**  
  - `method`: e.g., `'slice'`, `'cmss'`, or `'srss'`.  
  - `Ns`: total samples.  
  - `Nb`: burn-in samples.  
  Initializes from `thetaR_MAP` and iterates to draw from the posterior distribution of parameters considered as random variables

- **`display`**, **`chunkSize`**  
  If `display = true`, updates histograms every `chunkSize` samples.
  
- **`varR, varF`**  
  Remain your key index to which variables are random vs. fixed.

---

### 8. Key Rate Estimation

- **`EQ`, `Qs`**: Derived error and gain probabilities for each intensity level.  
- **`delta`**: The resulting error rate (QBER).  
- **`keyRate(Delta, delta, Qs)`**: Computes the secure key rate from these probabilities.  
- **`K_theory`, `Ks`**: The theoretical key rate (using true $\Delta$) vs. that from posterior draws (using the final MCMC samples).

---

### 9. Saving Data

- If `loadData = false` and `file_path` is non-empty:
  - Saves relevant items in `file_path.mat`.
- If `loadData = true`:
  - Skips the simulation and MCMC steps, loading them from the same `file_path.mat`.

---

### 10. Displaying Results

- **`displayHistograms(...)`**  
  Generates posterior histograms, credible intervals, or key-rate plots.  
- Figures may be printed to PDF if `savePlots = true`.

---

## Conclusion

In **`demo.m`**, you have the freedom to treat **any** QKD parameter (Alice’s, Bob’s, or Eve’s) as random or fixed. By adjusting the sets `varR` and `varF`, plus the prior/heuristic definitions, you can explore a wide range of scenarios—from partial to fully Bayesian inferences.