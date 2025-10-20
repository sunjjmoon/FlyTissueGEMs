%% 1. Check and export Pyruvate [c] Consuming Reactions

clc;
clear;
close all;

%% Settings
output_dir = '1_results_pyr_check';

pathway = pwd;
input_path = pathway;
output_path = fullfile(pathway, output_dir);

if ~exist(output_path, 'dir')
    mkdir(output_path);
end

%% Load analysis summary
fprintf('Loading analysis summary...\n');
load(fullfile(input_path, 'analysis_summary_sampling.mat'));
analysis_summary.valid_enzymes = analysis_summary.valid_enzymes(~strcmp(analysis_summary.valid_enzymes, 'MAR04388'));

%% Load individual condition files
fprintf('Loading condition files...\n');
individual_folder = fullfile(input_path, 'individual_conditions');
files = dir(fullfile(individual_folder, '*.mat'));

enzyme_data = struct();
for f = 1:length(files)
    try
        data = load(fullfile(individual_folder, files(f).name));
        condition = data.condition_data;
        enzyme = condition.target_enzyme;
        pert_factor = condition.perturbation_factor;
        
        if ~isfield(enzyme_data, enzyme)
            enzyme_data.(enzyme) = struct();
        end
        
        if pert_factor == 1.0
            enzyme_data.(enzyme).baseline = condition;
        else
            field_name = sprintf('pert_%.0f', pert_factor*100);
            enzyme_data.(enzyme).(field_name) = condition;
        end
    catch
    end
end
fprintf('Loaded %d condition files\n', length(files));

%% Load model to get stoichiometry
fprintf('Loading model for stoichiometry...\n');
model_path = 'E:\Projects\revision\pFBA_FBA_FVA\models';
load(fullfile(model_path, 'model_out_cbra_u.mat'));
hsd_model = model_out_cbra_u{2, 1};
fprintf('Model loaded: %d reactions, %d metabolites\n', length(hsd_model.rxns), length(hsd_model.mets));

%% Find pyruvate [c] and extract stoichiometry
fprintf('\nFinding pyruvate [c] in model...\n');
pyr_c_idx = find(strcmp(hsd_model.mets, 'MAM02819c[c]'));

if isempty(pyr_c_idx)
    error('Pyruvate [c] (MAM02819c) not found in model!');
end

fprintf('Found pyruvate [c] at index %d: %s\n', pyr_c_idx, hsd_model.metNames{pyr_c_idx});

% Get stoichiometry for all reactions
pyr_stoich = full(hsd_model.S(pyr_c_idx, :));

% Find reactions involving pyruvate [c]
pyr_rxn_indices = find(pyr_stoich ~= 0);
fprintf('Found %d reactions involving pyruvate [c]\n', length(pyr_rxn_indices));

%% Identify TRUE pyruvate [c] CONSUMING reactions
fprintf('\nIdentifying TRUE pyruvate [c] consuming reactions...\n');

first_enzyme = fieldnames(enzyme_data);
baseline_model = enzyme_data.(first_enzyme{1}).baseline.subSysModel;

pyruvate_consuming_reactions = {};
reaction_info = struct();

for i = 1:length(pyr_rxn_indices)
    rxn_idx = pyr_rxn_indices(i);
    rxn_id = hsd_model.rxns{rxn_idx};
    stoich = pyr_stoich(rxn_idx);
    is_reversible = hsd_model.rev(rxn_idx) == 1;
    
    % Check if reaction exists in subsystem model
    model_idx = find(strcmp(baseline_model.rxns, rxn_id));
    if isempty(model_idx)
        continue;
    end
    
    % Get reaction equation for interpretation
    rxn_name = hsd_model.rxnNames{rxn_idx};
    
    % Store reaction info
    pyruvate_consuming_reactions{end+1} = rxn_id;
    reaction_info.(rxn_id).stoich = stoich;
    reaction_info.(rxn_id).reversible = is_reversible;
    reaction_info.(rxn_id).rxn_name = rxn_name;
    
    % Determine which flux direction consumes pyruvate [c]
    % Negative stoich: pyruvate is consumed when flux is positive
    % Positive stoich: pyruvate is consumed when flux is negative (reverse)
    if stoich < 0
        reaction_info.(rxn_id).consuming_direction = 'positive';
        reaction_info.(rxn_id).interpretation = sprintf('Positive flux consumes pyruvate (stoich=%.2f)', stoich);
    else
        reaction_info.(rxn_id).consuming_direction = 'negative';
        reaction_info.(rxn_id).interpretation = sprintf('Negative flux consumes pyruvate (stoich=%.2f)', stoich);
    end
    
    % Create pseudo-equation for display
    if stoich < 0
        reaction_info.(rxn_id).equation = sprintf('... + %.1f pyruvate [c] -> ... (%s)', abs(stoich), rxn_name);
        reaction_info.(rxn_id).left_side = 'pyruvate [c] (substrate)';
        reaction_info.(rxn_id).right_side = 'products';
        reaction_info.(rxn_id).pyr_c_on_left = true;
        reaction_info.(rxn_id).pyr_c_on_right = false;
    else
        reaction_info.(rxn_id).equation = sprintf('... -> ... + %.1f pyruvate [c] (%s)', stoich, rxn_name);
        reaction_info.(rxn_id).left_side = 'substrates';
        reaction_info.(rxn_id).right_side = 'pyruvate [c] (product)';
        reaction_info.(rxn_id).pyr_c_on_left = false;
        reaction_info.(rxn_id).pyr_c_on_right = true;
    end
end

fprintf('Found %d potential consuming reactions\n', length(pyruvate_consuming_reactions));

%% Extract flux values and FILTER for actual consumption
fprintf('\nExtracting flux values and filtering for actual consumers...\n');

valid_enzymes = analysis_summary.valid_enzymes;
perturbation_factor = analysis_summary.perturbation_factors(1);
field_name = sprintf('pert_%.0f', perturbation_factor*100);

export_data = struct();
n_exported = 0;

for e = 1:length(valid_enzymes)
    target_enzyme = valid_enzymes{e};
    
    if ~isfield(enzyme_data, target_enzyme)
        continue;
    end
    
    baseline = enzyme_data.(target_enzyme).baseline;
    if ~isfield(enzyme_data.(target_enzyme), field_name)
        continue;
    end
    perturbed = enzyme_data.(target_enzyme).(field_name);
    
    for r = 1:length(pyruvate_consuming_reactions)
        rxn_id = pyruvate_consuming_reactions{r};
        
        idx_base = find(strcmp(baseline.subSysModel.rxns, rxn_id));
        idx_pert = find(strcmp(perturbed.subSysModel.rxns, rxn_id));
        
        if isempty(idx_base) || isempty(idx_pert)
            continue;
        end
        
        baseline_samples = baseline.samples(idx_base, :);
        perturbed_samples = perturbed.samples(idx_pert, :);
        
        baseline_median = median(baseline_samples);
        perturbed_median = median(perturbed_samples);
        
        % Check if flux is actually consuming
        consuming_direction = reaction_info.(rxn_id).consuming_direction;
        stoich = reaction_info.(rxn_id).stoich;
        
        is_consuming = false;
        
        if strcmp(consuming_direction, 'positive')
            % Negative stoich: need positive flux to consume
            if baseline_median > 1e-10
                is_consuming = true;
                % Keep samples as-is (positive = consumption)
            end
        else
            % Positive stoich: need negative flux to consume  
            if baseline_median < -1e-10
                is_consuming = true;
                % Flip sign for analysis (make consumption positive)
                baseline_samples = -baseline_samples;
                perturbed_samples = -perturbed_samples;
                baseline_median = -baseline_median;
                perturbed_median = -perturbed_median;
            end
        end
        
        % ONLY export if actually consuming
        if ~is_consuming
            continue;
        end
        
        % Store data
        data_key = sprintf('%s_%s', target_enzyme, rxn_id);
        n_exported = n_exported + 1;
        
        export_data.(data_key).target_enzyme = target_enzyme;
        export_data.(data_key).reaction_id = rxn_id;
        export_data.(data_key).equation = reaction_info.(rxn_id).equation;
        export_data.(data_key).interpretation = reaction_info.(rxn_id).interpretation;
        export_data.(data_key).baseline_median = baseline_median;
        export_data.(data_key).baseline_mean = mean(baseline_samples);
        export_data.(data_key).baseline_std = std(baseline_samples);
        export_data.(data_key).perturbed_median = perturbed_median;
        export_data.(data_key).perturbed_mean = mean(perturbed_samples);
        export_data.(data_key).perturbed_std = std(perturbed_samples);
        export_data.(data_key).baseline_samples = baseline_samples;
        export_data.(data_key).perturbed_samples = perturbed_samples;
    end
end

fprintf('Exported %d enzyme × reaction combinations (actual consumers only)\n', n_exported);

%% Export to Excel
fprintf('\nExporting to Excel...\n');

data_keys = fieldnames(export_data);
n_rows = length(data_keys);

Target_Enzyme = cell(n_rows, 1);
Reaction_ID = cell(n_rows, 1);
Equation = cell(n_rows, 1);
Interpretation = cell(n_rows, 1);
Baseline_Median = zeros(n_rows, 1);
Baseline_Mean = zeros(n_rows, 1);
Baseline_Std = zeros(n_rows, 1);
Perturbed_Median = zeros(n_rows, 1);
Perturbed_Mean = zeros(n_rows, 1);
Perturbed_Std = zeros(n_rows, 1);
Change_Percent = zeros(n_rows, 1);

for i = 1:n_rows
    key = data_keys{i};
    data = export_data.(key);
    
    Target_Enzyme{i} = data.target_enzyme;
    Reaction_ID{i} = data.reaction_id;
    Equation{i} = data.equation;
    Interpretation{i} = data.interpretation;
    Baseline_Median(i) = data.baseline_median;
    Baseline_Mean(i) = data.baseline_mean;
    Baseline_Std(i) = data.baseline_std;
    Perturbed_Median(i) = data.perturbed_median;
    Perturbed_Mean(i) = data.perturbed_mean;
    Perturbed_Std(i) = data.perturbed_std;
    
    if abs(data.baseline_median) > 1e-10
        Change_Percent(i) = ((data.perturbed_median - data.baseline_median) / data.baseline_median) * 100;
    else
        Change_Percent(i) = NaN;
    end
end

export_table = table(Target_Enzyme, Reaction_ID, Equation, Interpretation, ...
                     Baseline_Median, Baseline_Mean, Baseline_Std, ...
                     Perturbed_Median, Perturbed_Mean, Perturbed_Std, Change_Percent);

excel_file = fullfile(output_path, 'pyruvate_c_consumer_flux.xlsx');
writetable(export_table, excel_file);
fprintf('Exported to: %s\n', excel_file);


%% Create Reaction Summary Sheet with Consumption Rates
fprintf('\nCreating reaction consumption rate summary...\n');

unique_rxns = unique(Reaction_ID);
n_unique_rxns = length(unique_rxns);

Rxn_Summary_ID = cell(n_unique_rxns, 1);
Rxn_Summary_Stoich = zeros(n_unique_rxns, 1);
Rxn_Summary_Direction = cell(n_unique_rxns, 1);
Rxn_Summary_Avg_Baseline = zeros(n_unique_rxns, 1);
Rxn_Summary_Avg_Perturbed = zeros(n_unique_rxns, 1);
Rxn_Summary_N_Enzymes = zeros(n_unique_rxns, 1);

for r = 1:n_unique_rxns
    rxn_id = unique_rxns{r};
    idx = strcmp(Reaction_ID, rxn_id);
    
    Rxn_Summary_ID{r} = rxn_id;
    Rxn_Summary_Stoich(r) = reaction_info.(rxn_id).stoich;
    Rxn_Summary_Direction{r} = reaction_info.(rxn_id).consuming_direction;
    Rxn_Summary_Avg_Baseline(r) = mean(Baseline_Median(idx));
    Rxn_Summary_Avg_Perturbed(r) = mean(Perturbed_Median(idx));
    Rxn_Summary_N_Enzymes(r) = sum(idx);
end

reaction_summary_table = table(Rxn_Summary_ID, Rxn_Summary_Stoich, ...
                              Rxn_Summary_Direction, Rxn_Summary_Avg_Baseline, ...
                              Rxn_Summary_Avg_Perturbed, Rxn_Summary_N_Enzymes, ...
                              'VariableNames', {'Reaction_ID', 'Stoichiometry', ...
                                               'Consuming_Direction', 'Avg_Baseline_Consumption', ...
                                               'Avg_Perturbed_Consumption', 'N_Enzymes_Tested'});

% Sort by consumption rate
[~, sort_idx] = sort(Rxn_Summary_Avg_Baseline, 'descend');
reaction_summary_table = reaction_summary_table(sort_idx, :);

% Add to Excel file
writetable(reaction_summary_table, excel_file, 'Sheet', 'Reaction_Consumption_Summary');
fprintf('Added Reaction_Consumption_Summary sheet to Excel\n');

%% Save MAT file
mat_file = fullfile(output_path, 'pyruvate_c_reactions_full_data_v3.mat');
save(mat_file, 'export_data', 'reaction_info', 'pyruvate_consuming_reactions');
fprintf('Saved MAT file: %s\n', mat_file);

%% Print summary
fprintf('\n=== SUMMARY ===\n');
fprintf('Total consuming reactions identified: %d\n', length(pyruvate_consuming_reactions));
fprintf('Total enzyme × reaction combinations: %d\n', n_rows);
fprintf('Output directory: %s\n', output_path);

fprintf('\n=== TRUE CONSUMING REACTIONS ===\n');
unique_rxns = unique(Reaction_ID);
for r = 1:length(unique_rxns)
    rxn_id = unique_rxns{r};
    info = reaction_info.(rxn_id);
    idx = strcmp(Reaction_ID, rxn_id);
    avg_baseline = mean(Baseline_Median(idx));
    fprintf('%s (avg baseline=%.2f, stoich=%.2f)\n', rxn_id, avg_baseline, info.stoich);
    fprintf('  %s\n', info.equation);
    fprintf('  %s\n\n', info.interpretation);
end

fprintf('DONE!\n');