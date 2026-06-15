# Brain Dynamic Metric Toolbox (MATLAB)

This repository contains a suite of MATLAB functions designed to extract complex, nonlinear biomarkers from time-series data (such as neurophysiological signals, EEG, or other physiological responses). 

The toolkit focuses on quantifying system dynamics, latency, and chaotic behavior using three core metrics: the **Perturbative Integration Latency Index (PILI)**, the **Lyapunov Exponent (LE)**, and the **Autocorrelation Function Width (ACFW)**.

---

## 🛠 Features Included

### 1. Perturbative Integration Latency Index (`getPILI.m`)
Calculates the PILI by analyzing the area under the curve (AUC) of a normalized, baseline-corrected signal. It fits a trapezoidal integration over the time span of the response to quantify how a system recovers or responds to perturbations over time.
* **Outputs:** `biomarkers.pili` (normalized AUC) and `biomarkers.pili2` (absolute AUC).

### 2. Lyapunov Exponent (`getLE.m`)
Estimates the largest Lyapunov Exponent to measure the chaotic nature and predictability of the system. It normalizes the data via z-scoring and relies on Wolf's algorithm for phase-space reconstruction to calculate the exponential divergence of nearby trajectories.
* **Outputs:** `biomarkers.LE_est`
* *Note:* This function requires external helper functions (`basgen_new.m` and `fet.m` from https://au.mathworks.com/matlabcentral/fileexchange/48084-wolf-lyapunov-exponent-estimation-from-a-time-series) to generate the database and evolve the trajectories.

### 3. Autocorrelation Function Width (`getACFW.m`)
Measures the "memory" of a signal by calculating its autocorrelation function. It detrends and normalizes the data, then counts the number of lags where the autocorrelation coefficient remains strictly above a specific threshold (0.6).
* **Outputs:** `biomarkers.arE`

---

## ⚙️ Prerequisites and Dependencies

* **MATLAB:** Requires a standard MATLAB installation. 
* **Toolboxes:** The `autocorr` function used in `getACFW.m` may require the Econometrics Toolbox or Signal Processing Toolbox depending on your MATLAB version.
* **External Functions:** `getLE.m` requires an implementation of Wolf's algorithm. You must ensure that `basgen_new.m` and `fet.m` are present in your MATLAB path.

---

## 🚀 Usage

All three functions update and return a central `biomarkers` structure, making it easy to loop through multiple trials, channels, or datasets (indicated by the index `d`).

### Common Input Parameters
* `fs`: Sampling frequency (Hz) — *Required for `getPILI` only*.
* `response`: A 2D matrix of your time-series data `[trials/channels x time]`.
* `d`: The current iteration or dimension index to store the calculated biomarker in the structure.
* `biomarkers`: A MATLAB struct used to accumulate the results.
* `nonartefact_length`: An array of indices representing the clean, artifact-free portion of the signal to be analyzed.

---

## 📂 Data Formatting Notes

* **Baseline Correction:** `getPILI` automatically attempts to remove the DC offset by calculating the baseline from the final 10% of the provided data block. 
* **Detrending:** `getACFW` performs linear detrending natively before computing the autocorrelation to ensure drifting baselines do not artificially inflate the memory of the signal.
