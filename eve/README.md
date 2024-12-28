# Inference of Eve

This script (`demo.m`) demonstrates parameter inference when **only Eve’s parameters** $\theta_E$ are treated as random, while all other system parameters remain fixed. Below is an outline of the major steps the script performs, followed by a detailed listing of the variables used in each step.

---

## Script Outline

1. **Environment Setup**  
   Initializes the MATLAB environment (path, random seed, etc.) so that results are reproducible and other required files can be accessed.

2. **System Parameters**  
   Defines and configures $\theta_A$ (Alice’s parameters), $\theta_B$ (Bob’s parameters), and $\theta_E$ (Eve’s parameters). Also applies heuristics (e.g., to set $\lambda$ values, or find the maximum $k$).

3. **Prior Configuration**  
   Sets bounds and shape parameters for the Bayesian priors on Eve’s random parameters $\theta_E$.

4. **Simulation**  
   Generates synthetic click counts based on the fixed system parameters (Alice/Bob) and Eve’s “true” values.

5. **MAP Estimation**  
   Performs a maximum a posteriori search for $\theta_E$, then compares the results with Eve’s ground-truth values used in the simulation.

6. **MCMC Sampling**  
   Initializes from the MAP point and explores the posterior distribution of $\theta_E$ using a chosen sampler (e.g., slice sampling).

7. **Key Rate Estimation**  
   Computes error/gain probabilities and uses them to derive secure key rates under both the true and inferred parameters.

8. **Saving Data**  
   Stores simulation outputs, MAP estimates, and posterior samples in a `.mat` file (if specified).

9. **Displaying Results**  
   Plots histograms of Eve’s inferred parameters and the resulting key rates, optionally saving figures to PDF.

---

## Detailed Variables by Section

### 1. Environment Setup

- **`seed`**: An integer controlling `rng(seed)`. Setting this ensures reproducible randomness.

- **`file_path`**: A string specifying where to save/load `.mat` data (e.g. `'demo_data'`). If empty (`[]`), results are not saved.

- **`loadData`** (boolean)  
  - `true`: Skip simulation/MCMC and load existing data from `file_path`.  
  - `false`: Run the full process and save new results.

- **`savePlots`** (boolean)  
  - `true`: Export final plots (PDF).  
  - `false`: No file output for figures.

---

### 2. System Parameters

#### 2.1 Alice’s Parameters ($\theta_A$)

- **`lambdas`**: A vector of intensity levels, derived from a heuristic function (e.g., `get_lambdas`) that spans $\lambda_{\min}$ to a `limit`.  
- **`alpha`**: The fiber attenuation coefficient (dB/km).  
- **`dAB`**: Distance between Alice and Bob (km).  
- These are grouped as $\theta_A = \\{\Lambda, \alpha, d_{AB}\\}$.

#### 2.2 Bob’s Parameters ($\theta_B$)

- **`pa0`, `pa1`**: After-pulsing probabilities for Bob’s two detectors (with slight ±10% variations).  
- **`pc0`, `pc1`**: Detector efficiencies (again ±10% variations).  
- **`pd0`, `pd1`**: Dark count probabilities (±10%).  
- **`pe`**: Misalignment probability.  
- These remain fixed and are grouped as $\theta_B = \\{pa0, pa1, pc0, pc1, pd0, pd1, pe\\}$.

#### 2.3 Eve’s Parameters ($\theta_E$)

- **`dAE`**: Eve’s distance from Alice (km).  
- **`pEB`**: Eve’s channel efficiency, computed via a heuristic (e.g., `get_pEB`) so Eve remains inconspicuous for a chosen $k$.  
- **`k`**: Number of photons intercepted per pulse.  
- **`Delta`**: Fraction of pulses intercepted by Eve.  
- These four are the **random** parameters we infer: $\theta_E = \\{d_{AE}, p_{EB}, k, \Delta\\}$.

---

### 3. Prior Configuration

- **$\alpha$, $\beta$**: Shape parameters for Beta/Gamma priors on each of Eve’s parameters.  
- **`ub`, `lb`**: Vectors of upper/lower bounds for $d_{AE}$, $p_{EB}$, $k$, $\Delta$.  
- Grouped as $\theta_P = \\{\alpha, \beta, \text{ub}, \text{lb}\\}$.

---

### 4. Simulation

- **`simulate(N, thetaA, thetaB, thetaE)`**: Produces synthetic click counts $\mathbf{C}$, given the total pulses `N`, plus the system parameters for Alice, Bob, and Eve (the “true” $\theta_E$ values).

---

### 5. MAP Estimation

- **`MAP(...)`**: Uses an optimization routine (e.g., `'quasi-newton'`) up to `maxIters` iterations to find the $\theta_E$ maximizing the posterior given $\mathbf{C}$.  
- **`thetaE_MAP`**: The resulting MAP estimate, which is then compared to the ground-truth $\theta_E$.

---

### 6. MCMC Sampling

- **`sample(method, thetaE_MAP, Ns, Nb, ...)`**  
  - `method`: e.g., `'slice'`, `'cmss'`, or `'srss'`.  
  - `Ns`: total samples.  
  - `Nb`: burn-in samples.  
  Initializes from `thetaE_MAP` and iterates to draw from the posterior distribution of Eve’s parameters.  

- **`display`**, **`chunkSize`**  
  If `display = true`, updates histograms every `chunkSize` samples.

---

### 7. Key Rate Estimation

- **`EQ`, `Qs`**: Derived error and gain probabilities for each intensity level.  
- **`delta`**: The resulting error rate (QBER).  
- **`keyRate(Delta, delta, Qs)`**: Computes the secure key rate from these probabilities.  
- **`K_theory`, `Ks`**: The theoretical key rate (using true $\Delta$) vs. that from posterior draws (using the final MCMC samples).

---

### 8. Saving Data

- If `loadData = false` and `file_path` is non-empty:
  - Saves **`samples`, `thetaE_MAP`, `K_theory`, `Ks`**, and other relevant items in `.mat` format.
- If `loadData = true`:
  - Skips the simulation and MCMC steps, loading them from the same file_path.

---

### 9. Displaying Results

- **`displayHistograms(...)`**  
  Generates posterior histograms, credible intervals, or key-rate plots.  
- Figures may be printed to PDF if `savePlots = true`.

---

## Conclusion

By following these steps, **`demo.m`** provides a workflow for simulating, inferring, and analyzing Eve’s QKD attack parameters under the assumption that Bob’s and Alice’s parameters are fixed. You can adjust any of the above variables (e.g., number of pulses `N`, intensities `lambdas`, or priors) to explore various intrusion scenarios. For a more general approach that also randomizes Alice/Bob parameters, refer to the **`full/`** folder.
