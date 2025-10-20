%% Script 04: FBA vs pFBA vs Sampling Comparison

clc;
clear;
close all;

%% Setup directories
current_path = pwd;
results_dir = fullfile(current_path, '04_sampling_comp'); % 
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

% sampling_dir = fullfile('C:\Users\Sun Jin\Documents\revision\tissueGEM_revision\v2\code_to_upload\3_Flux_analysis_orig_v2_f', 'new_different_20251001','C_sampling');
sampling_dir = fullfile(current_path,'models');

%% Load analysis results
disp('Loading FBA and pFBA results...');
load(fullfile(current_path,'01_results','model_constrained_out.mat'));
load(fullfile(current_path,'02_comp_fba_pFBA','FBA_pFBA_comparison.mat'));

%% Load sampling results
disp('Loading sampling results...');
try
    nsd_sampling = load(fullfile(sampling_dir, 'nsd.mat'));
    hsd_sampling = load(fullfile(sampling_dir, 'hsd.mat'));
    disp('Sampling data loaded successfully.');
catch ME
    disp(['Error loading sampling data: ' ME.message]);
    disp(['Please check if files exist at: ' sampling_dir]);
    return;
end

%% Extract data structures
models = {'NSD', 'HSD'};
sampling_data = {nsd_sampling, hsd_sampling};

%% Basic data inspection
disp('=== Data Structure Inspection ===');

% Check FBA/pFBA data
for i = 1:2
    disp([newline models{i} ' Model:']);
    disp(['  FBA solution size: ' num2str(length(comparison_results{i}.FBA.solution.x)) ' reactions']);
    disp(['  pFBA solution size: ' num2str(length(comparison_results{i}.pFBA.solution.x)) ' reactions']);
    disp(['  Model reactions: ' num2str(length(model_constrained_out{i}.rxns))]);
end

% Check sampling data structure
disp([newline 'Sampling Data Structure:']);
for i = 1:2
    disp([newline models{i} ' Sampling:']);
    sampling_fields = fieldnames(sampling_data{i});
    for j = 1:length(sampling_fields)
        field = sampling_fields{j};
        data = sampling_data{i}.(field);
        if isnumeric(data)
            disp(['  ' field ': ' num2str(size(data,1)) 'x' num2str(size(data,2)) ' numeric array']);
        else
            disp(['  ' field ': ' class(data)]);
        end
    end
end

%% Prepare comparison data structure
comparison_data = struct();

for i = 1:2
    model_id = models{i};
    
    % FBA and pFBA fluxes
    fba_fluxes = comparison_results{i}.FBA.solution.x;
    pfba_fluxes = comparison_results{i}.pFBA.solution.x;
    
    % Store in comparison structure
    comparison_data.(model_id).fba_fluxes = fba_fluxes;
    comparison_data.(model_id).pfba_fluxes = pfba_fluxes;
    comparison_data.(model_id).rxn_names = model_constrained_out{i}.rxns;
    
    % Extract sampling data from the 'samples.points' field
    if isfield(sampling_data{i}, 'samples') && isfield(sampling_data{i}.samples, 'points')
        sampling_matrix = sampling_data{i}.samples.points;
        disp(['Found sampling data: ' num2str(size(sampling_matrix,1)) ' reactions x ' num2str(size(sampling_matrix,2)) ' samples for ' model_id]);
    else
        disp(['Warning: Could not find sampling matrix for ' model_id]);
        sampling_matrix = zeros(length(fba_fluxes), 100);
    end
    
    comparison_data.(model_id).sampling_matrix = sampling_matrix;
    comparison_data.(model_id).n_samples = size(sampling_matrix, 2);
    
    % Calculate sampling statistics
    comparison_data.(model_id).sampling_mean = mean(sampling_matrix, 2);
    comparison_data.(model_id).sampling_std = std(sampling_matrix, 0, 2);
    comparison_data.(model_id).sampling_median = median(sampling_matrix, 2);
    comparison_data.(model_id).sampling_min = min(sampling_matrix, [], 2);
    comparison_data.(model_id).sampling_max = max(sampling_matrix, [], 2);
end

%% Summary statistics
disp([newline '=== Summary Statistics ===']);
for i = 1:2
    model_id = models{i};
    data = comparison_data.(model_id);
    
    disp([newline model_id ' Model:']);
    disp(['  Number of reactions: ' num2str(length(data.fba_fluxes))]);
    disp(['  Number of samples: ' num2str(data.n_samples)]);
    disp(['  FBA total flux: ' num2str(sum(abs(data.fba_fluxes)))]);
    disp(['  pFBA total flux: ' num2str(sum(abs(data.pfba_fluxes)))]);
    disp(['  Sampling mean total flux: ' num2str(sum(abs(data.sampling_mean)))]);
    disp(['  Sampling std total flux: ' num2str(sum(data.sampling_std))]);
end

%% Method Correlation Analysis
disp([newline '=== Method Correlation Analysis ===']);

correlation_results = struct();

for i = 1:2
    model_id = models{i};
    data = comparison_data.(model_id);
    
    disp([newline model_id ' Model Correlations:']);
    
    % Get absolute flux values for correlation analysis
    fba_abs = abs(data.fba_fluxes);
    pfba_abs = abs(data.pfba_fluxes);
    sampling_mean_abs = abs(data.sampling_mean);
    
    % Calculate correlations (using only non-zero fluxes)
    % Remove reactions with zero flux in all methods
    nonzero_idx = fba_abs > 1e-6 | pfba_abs > 1e-6 | sampling_mean_abs > 1e-6;
    
    fba_nz = fba_abs(nonzero_idx);
    pfba_nz = pfba_abs(nonzero_idx);
    sampling_nz = sampling_mean_abs(nonzero_idx);
    
    % Correlation coefficients
    corr_fba_pfba = corr(fba_nz, pfba_nz);
    corr_fba_sampling = corr(fba_nz, sampling_nz);
    corr_pfba_sampling = corr(pfba_nz, sampling_nz);
    
    disp(['  FBA vs pFBA correlation: ' num2str(corr_fba_pfba, '%.3f')]);
    disp(['  FBA vs Sampling correlation: ' num2str(corr_fba_sampling, '%.3f')]);
    disp(['  pFBA vs Sampling correlation: ' num2str(corr_pfba_sampling, '%.3f')]);
    
    % Store results
    correlation_results.(model_id).correlations = [corr_fba_pfba, corr_fba_sampling, corr_pfba_sampling];
    correlation_results.(model_id).nonzero_reactions = sum(nonzero_idx);
    correlation_results.(model_id).total_reactions = length(nonzero_idx);
    
    % Agreement analysis - find reactions where methods disagree significantly
    % Define disagreement as >2-fold difference
    fold_threshold = 2;
    
    % FBA vs pFBA disagreement
    fba_pfba_ratio = (fba_nz + 1e-6) ./ (pfba_nz + 1e-6);
    high_disagreement_fba_pfba = sum(fba_pfba_ratio > fold_threshold | fba_pfba_ratio < 1/fold_threshold);
    
    % FBA vs Sampling disagreement
    fba_sampling_ratio = (fba_nz + 1e-6) ./ (sampling_nz + 1e-6);
    high_disagreement_fba_sampling = sum(fba_sampling_ratio > fold_threshold | fba_sampling_ratio < 1/fold_threshold);
    
    % pFBA vs Sampling disagreement
    pfba_sampling_ratio = (pfba_nz + 1e-6) ./ (sampling_nz + 1e-6);
    high_disagreement_pfba_sampling = sum(pfba_sampling_ratio > fold_threshold | pfba_sampling_ratio < 1/fold_threshold);
    
    disp(['  Reactions with >2-fold FBA/pFBA difference: ' num2str(high_disagreement_fba_pfba) ' (' num2str(100*high_disagreement_fba_pfba/length(fba_nz), '%.1f') '%)']);
    disp(['  Reactions with >2-fold FBA/Sampling difference: ' num2str(high_disagreement_fba_sampling) ' (' num2str(100*high_disagreement_fba_sampling/length(fba_nz), '%.1f') '%)']);
    disp(['  Reactions with >2-fold pFBA/Sampling difference: ' num2str(high_disagreement_pfba_sampling) ' (' num2str(100*high_disagreement_pfba_sampling/length(fba_nz), '%.1f') '%)']);
    
    % Store disagreement results
    correlation_results.(model_id).disagreements = [high_disagreement_fba_pfba, high_disagreement_fba_sampling, high_disagreement_pfba_sampling];
    correlation_results.(model_id).disagreement_percent = [high_disagreement_fba_pfba, high_disagreement_fba_sampling, high_disagreement_pfba_sampling] ./ length(fba_nz) * 100;
end

%% Flux Range Comparison
disp([newline '=== Flux Range Comparison ===']);

for i = 1:2
    model_id = models{i};
    data = comparison_data.(model_id);
    
    disp([newline model_id ' Model Flux Ranges:']);
    
    % Calculate flux ranges (exclude near-zero values)
    threshold = 1e-6;
    
    fba_nonzero = abs(data.fba_fluxes(abs(data.fba_fluxes) > threshold));
    pfba_nonzero = abs(data.pfba_fluxes(abs(data.pfba_fluxes) > threshold));
    sampling_nonzero = abs(data.sampling_mean(abs(data.sampling_mean) > threshold));
    
    disp(['  FBA: ' num2str(length(fba_nonzero)) ' active reactions, range: ' num2str(min(fba_nonzero), '%.2e') ' to ' num2str(max(fba_nonzero), '%.2e')]);
    disp(['  pFBA: ' num2str(length(pfba_nonzero)) ' active reactions, range: ' num2str(min(pfba_nonzero), '%.2e') ' to ' num2str(max(pfba_nonzero), '%.2e')]);
    disp(['  Sampling: ' num2str(length(sampling_nonzero)) ' active reactions, range: ' num2str(min(sampling_nonzero), '%.2e') ' to ' num2str(max(sampling_nonzero), '%.2e')]);
    
    % Check for potential cycling reactions (high sampling variance)
    high_variance_idx = data.sampling_std > 100; % Reactions with high variance
    low_pfba_idx = abs(data.pfba_fluxes) < 10;    % But low pFBA flux
    potential_cycling = high_variance_idx & low_pfba_idx;
    
    disp(['  Potential cycling reactions (high variance, low pFBA): ' num2str(sum(potential_cycling))]);
    
    % Store range results
    correlation_results.(model_id).flux_ranges.fba = [min(fba_nonzero), max(fba_nonzero)];
    correlation_results.(model_id).flux_ranges.pfba = [min(pfba_nonzero), max(pfba_nonzero)];
    correlation_results.(model_id).flux_ranges.sampling = [min(sampling_nonzero), max(sampling_nonzero)];
    correlation_results.(model_id).potential_cycling_count = sum(potential_cycling);
end

%% Disease Difference Comparison (HSD - NSD)

% Calculate disease differences for each method
fba_disease_diff = abs(comparison_data.HSD.fba_fluxes) - abs(comparison_data.NSD.fba_fluxes);
pfba_disease_diff = abs(comparison_data.HSD.pfba_fluxes) - abs(comparison_data.NSD.pfba_fluxes);
sampling_disease_diff = abs(comparison_data.HSD.sampling_mean) - abs(comparison_data.NSD.sampling_mean);

% Calculate correlations between disease difference patterns
corr_disease_fba_pfba = corr(fba_disease_diff, pfba_disease_diff);
corr_disease_fba_sampling = corr(fba_disease_diff, sampling_disease_diff);
corr_disease_pfba_sampling = corr(pfba_disease_diff, sampling_disease_diff);

disp('Disease difference correlations:');
disp(['  FBA vs pFBA: ' num2str(corr_disease_fba_pfba, '%.3f')]);
disp(['  FBA vs Sampling: ' num2str(corr_disease_fba_sampling, '%.3f')]);
disp(['  pFBA vs Sampling: ' num2str(corr_disease_pfba_sampling, '%.3f')]);

% Store disease difference results
correlation_results.disease_differences.correlations = [corr_disease_fba_pfba, corr_disease_fba_sampling, corr_disease_pfba_sampling];
correlation_results.disease_differences.fba_diff = fba_disease_diff;
correlation_results.disease_differences.pfba_diff = pfba_disease_diff;
correlation_results.disease_differences.sampling_diff = sampling_disease_diff;

%% Save all results
save(fullfile(results_dir, 'fba_pfba_sampling_comparison.mat'), 'comparison_data', 'correlation_results', 'models');

disp([newline 'Method correlation analysis complete.']);
disp(['Results saved to: ' fullfile(results_dir, 'fba_pfba_sampling_comparison.mat')]);
