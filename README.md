# Brain Dynamic Metric Toolbox (MATLAB)

This repository contains a suite of MATLAB functions designed to extract complex, nonlinear biomarkers from time-series data. 

**Note on Applicability:** While this suite (BDMT - Brain Dynamics Metric Tool) was originally presented within the context of a Transcranial Magnetic Stimulation (TMS) case study, it is highly adaptable and can be readily applied to other modalities of brain stimulation or physiological responses wherever an estimation of brain dynamics is required.

---

## 🛠 Features Included

### 1. Perturbative Integration Latency Index (`getPILI.m`)
Calculates the PILI to quantify the residual energy and system recovery following a perturbation. 
* **Method:** The data is normalized (between -1 and 1) and baseline-corrected. The script computes the absolute value and calculates the area under the curve (AUC) using trapezoidal integration, with the signal at 0 acting as the baseline.
* **Outputs:** `biomarkers.pili` (normalized AUC) and `biomarkers.pili2` (absolute AUC).

### 2. Lyapunov Exponent (`getLE.m`)
Estimates the largest Lyapunov Exponent (LE) to measure the chaotic nature and oscillatory stability of the neural system without introducing normalization bias into the phase space.
* **Method:** Requires dimension ($m$) and lag ($\tau$) inputs to embed the data into a new phase-space dimension. By default, these are set to $m = 3$ and $\tau = 1$. For custom datasets, these parameters can be optimized using the Cao embedding algorithm (e.g., an `optim_m_tau` function).
* **Outputs:** `biomarkers.LE_est`
* **Attribution:** This function relies on helper scripts (`basgen_new.m` and `fet.m`) to generate the database and evolve the trajectories based on Wolf's algorithm (https://au.mathworks.com/matlabcentral/fileexchange/48084-wolf-lyapunov-exponent-estimation-from-a-time-series). These scripts are modified versions of the originals provided by Taehyeun Park (The Cooper Union, EE’15).

### 3. Autocorrelation Function Width (`getACFW.m`)
Measures the "memory" and interdependency of a signal.
* **Method:** The script detrends and normalizes the data, then uses MATLAB's native `autocorr` function. The width of the autocorrelation is determined by measuring the length (number of lags) where the correlation coefficient remains above a specific threshold (e.g., falling to half of its maximum value, or strictly >0.6).
* **Outputs:** `biomarkers.arE`

---

## ⚙️ Prerequisites and Dependencies

* **MATLAB:** Requires a standard MATLAB installation. 
* **Toolboxes:** The `autocorr` function used in `getACFW.m` requires the Econometrics Toolbox or Signal Processing Toolbox.
* **External Functions:** `getLE.m` requires `basgen_new.m` and `fet.m` to be present in your MATLAB path.

---

## 🚀 Usage & Downstream Workflow

All three functions update and return a central `biomarkers` structure, making it easy to loop through multiple subjects, trials, or channels.

### Standard Pipeline
1. **Extraction:** Loop your dataset through `getPILI.m`, `getLE.m`, and `getACFW.m` to populate the `biomarkers` struct.
2. **Export:** Export the final struct array to a `.csv` or `.xlsx` (Excel) format using MATLAB's `writetable()` or `writestruct()`.
3. **Statistical Analysis:** The exported data is formatted to be easily imported into **R**. From there, you can perform subsequent statistical analyses, such as applying Linear Mixed Models (LMM), to investigate significant interactions between stimulation parameters, time, and frequency.

### Common Input Parameters
* `fs`: Sampling frequency (Hz) — *Required for `getPILI` only*.
* `response`: A 2D matrix of your time-series data `[trials/channels x time]`.
* `d`: The current iteration or dimension index to store the calculated biomarker.
* `biomarkers`: A MATLAB struct used to accumulate the results.
* `nonartefact_length`: An array of indices representing the clean, artifact-free portion of the signal.
