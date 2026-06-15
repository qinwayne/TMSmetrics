%% Main Test Script for Biomarker Functions
clear; clc;

% 1. Setup Dummy Parameters
fs = 1000;                     % Sampling frequency: 1000 Hz
num_trials = 5;                % Number of simulated channels/trials
num_samples = 2000;            % Total length of the time-series response
d = 1;                         % Dimension/Iteration index

% Create an empty struct to store results
biomarkers = struct();

% Define the artifact-free index range 
% (e.g., ignoring the first 200 samples due to stimulation artifact)
nonartefact_length = 201:num_samples;

% 2. Generate Dummy Data
% Creating a simulated signal: a decaying sine wave with added noise
t = (0:num_samples-1) / fs;
clean_signal = exp(-2*t) .* sin(2*pi*10*t); 
noise = 0.1 * randn(num_trials, num_samples);

% Our dummy response matrix [trials x time]
response = repmat(clean_signal, num_trials, 1) + noise;

%% 3. Test Function 1: getPILI
fprintf('Running getPILI...\n');
try
    biomarkers = getPILI(fs, response, d, biomarkers, nonartefact_length);
    fprintf('  Success! PILI: %.4f | PILI2: %.4f\n', biomarkers.pili(d), biomarkers.pili2(d));
catch ME
    fprintf('  Error in getPILI: %s\n', ME.message);
end

%% 4. Test Function 2: getLE
fprintf('Running getLE...\n');
try
    biomarkers = getLE(response, d, biomarkers, nonartefact_length);
    fprintf('  Success! LE_est: %.4f\n', biomarkers.LE_est(1, d));
catch ME
    fprintf('  Error in getLE: %s\n', ME.message);
end

%% 5. Test Function 3: getACFW
fprintf('Running getACFW...\n');
try
    biomarkers = getACFW(response, d, biomarkers, nonartefact_length);
    fprintf('  Success! arE: %.4f\n', biomarkers.arE(d));
catch ME
    fprintf('  Error in getACFW: %s\n', ME.message);
end

%% 6. Display Final Struct
fprintf('\n--- Final Biomarkers Struct ---\n');
disp(biomarkers);
