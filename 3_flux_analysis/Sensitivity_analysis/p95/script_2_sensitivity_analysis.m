%% 2. Sensitivity analysis
clc;
clear;
close all;

%% Settings
output_dir = '2_results';
pathway = pwd;
input_path = pathway;
output_path = fullfile(pathway, output_dir);

if ~exist(output_path, 'dir')
    mkdir(output_path);
end

%% Load the data
fprintf('Loading exported pyruvate flux data (V3)...\n');
export_data_path = fullfile(input_path, '1_results_pyr_check', 'pyruvate_c_reactions_full_data_v3.mat');

if ~exist(export_data_path, 'file')
    error('Cannot find: %s\nPlease run the V3 export script first!', export_data_path);
end

load(export_data_path);
fprintf('Loaded data successfully!\n');

%% Prepare data for R analysis
fprintf('\nPreparing data for R analysis...\n');

% Extract all data keys
data_keys = fieldnames(export_data);
n_rows = length(data_keys);

% Initialize arrays for Individual_Reactions sheet
Enzyme = cell(n_rows, 1);
Pyruvate_Reaction = cell(n_rows, 1);
Baseline_Median = zeros(n_rows, 1);
Baseline_Mean = zeros(n_rows, 1);
Baseline_StdDev = zeros(n_rows, 1);
Perturbed_Median = zeros(n_rows, 1);
Perturbed_Mean = zeros(n_rows, 1);
Perturbed_StdDev = zeros(n_rows, 1);
Flux_Change = zeros(n_rows, 1);

for i = 1:n_rows
    key = data_keys{i};
    data = export_data.(key);
    
    Enzyme{i} = data.target_enzyme;
    Pyruvate_Reaction{i} = data.reaction_id;
    Baseline_Median(i) = data.baseline_median;
    Baseline_Mean(i) = data.baseline_mean;
    Baseline_StdDev(i) = data.baseline_std;
    Perturbed_Median(i) = data.perturbed_median;
    Perturbed_Mean(i) = data.perturbed_mean;
    Perturbed_StdDev(i) = data.perturbed_std;
    Flux_Change(i) = data.perturbed_median - data.baseline_median; 
end

% Create Individual_Reactions table
individual_reactions = table(Enzyme, Pyruvate_Reaction, ...
                             Baseline_Median, Baseline_Mean, Baseline_StdDev, ...
                             Perturbed_Median, Perturbed_Mean, Perturbed_StdDev, ...
                             Flux_Change);

fprintf('Individual Reactions table: %d rows\n', height(individual_reactions));

%% Calculate summary statistics by enzyme
fprintf('\nCalculating summary statistics...\n');

unique_enzymes = unique(Enzyme);
n_enzymes = length(unique_enzymes);

Summary_Enzyme = cell(n_enzymes, 1);
Total_Baseline_Median = zeros(n_enzymes, 1);
Total_Perturbed_Median = zeros(n_enzymes, 1);
Total_Flux_Change = zeros(n_enzymes, 1);
System_Response = zeros(n_enzymes, 1);
N_Reactions = zeros(n_enzymes, 1);

for e = 1:n_enzymes
    enzyme = unique_enzymes{e};
    idx = strcmp(Enzyme, enzyme);
    
    Summary_Enzyme{e} = enzyme;
    
    % Sum across all reactions for this enzyme (NO abs())
    Total_Baseline_Median(e) = sum(Baseline_Median(idx));
    Total_Perturbed_Median(e) = sum(Perturbed_Median(idx));
    Total_Flux_Change(e) = Total_Perturbed_Median(e) - Total_Baseline_Median(e);
    
    % Calculate system response
    if abs(Total_Baseline_Median(e)) > 1e-10
        System_Response(e) = Total_Flux_Change(e) / Total_Baseline_Median(e);
    else
        System_Response(e) = NaN;
    end
   
    N_Reactions(e) = sum(idx);
end

% Create Summary table
summary_table = table(Summary_Enzyme, Total_Baseline_Median, Total_Perturbed_Median, ...
                     Total_Flux_Change, System_Response, N_Reactions, ...
                     'VariableNames', {'Enzyme', 'Total_Baseline', 'Total_Perturbed', ...
                                      'Total_Flux_Change', 'System_Response', 'N_Reactions'});

fprintf('Summary table: %d enzymes\n', height(summary_table));

%% Perform sensitivity analysis
fprintf('\Perform sensitivity...\n');

Detailed_Enzyme = Summary_Enzyme;
Avg_Baseline_Pyruvate = Total_Baseline_Median;
Avg_Perturbed_Pyruvate = Total_Perturbed_Median;
Flux_Difference = Total_Flux_Change;
System_Response_Detailed = System_Response;

% Calculate sensitivity 
perturbation_factor = 0.95;
Sensitivity = zeros(n_enzymes, 1);
for e = 1:n_enzymes
    if abs(System_Response_Detailed(e)) > 1e-10
        Sensitivity(e) = System_Response_Detailed(e) / -(1 - perturbation_factor);
    else
        Sensitivity(e) = 0;
    end
end

Sensitivity_SE = zeros(n_enzymes, 1);
System_Response_StdDev = zeros(n_enzymes, 1);

for e = 1:n_enzymes
    enzyme = unique_enzymes{e};
    idx = strcmp(Enzyme, enzyme);
    
    flux_changes = Flux_Change(idx);
    baseline_fluxes = Baseline_Median(idx);
    
    % System response for each reaction (NO abs() on baseline)
    system_responses = flux_changes ./ (baseline_fluxes + eps);
    sensitivity_across_rxns = system_responses ./ -(1 - perturbation_factor);
    
    % Standard deviation and standard error
    System_Response_StdDev(e) = std(sensitivity_across_rxns);
    Sensitivity_SE(e) = System_Response_StdDev(e) / sqrt(sum(idx));
end

detailed_summary = table(Detailed_Enzyme, Avg_Baseline_Pyruvate, Avg_Perturbed_Pyruvate, ...
                        Flux_Difference, System_Response_Detailed, System_Response_StdDev, ...
                        Sensitivity, Sensitivity_SE, ...
                        'VariableNames', {'Enzyme', 'Avg_Baseline_Pyruvate', 'Avg_Perturbed_Pyruvate', ...
                                         'Flux_Difference', 'System_Response', 'System_Response_StdDev', ...
                                         'Sensitivity', 'Sensitivity_SE'});

fprintf('Detailed summary table: %d enzymes\n', height(detailed_summary));

%% Create Raw_Sample_Data (sample-level data)
fprintf('\nCreating raw sample data...\n');

Raw_Sample_Enzyme = cell(0, 1);
Raw_Sample_Index = [];
Raw_Baseline_Total = [];
Raw_Perturbed_Total = [];
Raw_Sensitivity_Sample = [];

% Create empty table with proper structure
raw_sample_data = table(Raw_Sample_Enzyme, Raw_Sample_Index, Raw_Baseline_Total, ...
                       Raw_Perturbed_Total, Raw_Sensitivity_Sample, ...
                       'VariableNames', {'Enzyme', 'Sample_Index', 'Baseline_Total_Pyruvate', ...
                                        'Perturbed_Total_Pyruvate', 'Sensitivity_Sample'});

fprintf('Raw sample data: %d samples (placeholder)\n', height(raw_sample_data));

%% Export to Excel with multiple sheets
fprintf('\nExporting to Excel...\n');

excel_file = fullfile(output_path, 'pyruvate_consumption_sensitivity_results_v3.xlsx');

% Write each sheet
writetable(detailed_summary, excel_file, 'Sheet', 'Detailed_Summary');
writetable(raw_sample_data, excel_file, 'Sheet', 'Raw_Sample_Data');
writetable(individual_reactions, excel_file, 'Sheet', 'Individual_Reactions');
writetable(summary_table, excel_file, 'Sheet', 'Summary');

fprintf('Exported to: %s\n', excel_file);

%% Display top results
fprintf('\n=== TOP 5 MOST SENSITIVE ENZYMES ===\n');
[~, sort_idx] = sort(Sensitivity, 'descend');
for i = 1:min(5, length(sort_idx))
    idx = sort_idx(i);
    fprintf('%d. %s: Sensitivity = %.4f ± %.4f\n', i, Detailed_Enzyme{idx}, ...
            Sensitivity(idx), Sensitivity_SE(idx));
end

fprintf('\nDONE!\n');