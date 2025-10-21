%% script_02_z_test_sampling_nad
% Z-Test on NAD(H)-dependent reactions for Sampling 
% Purpose: Perform reaction-level Z-test on NAD(H)-dependent reactions to identify 
% statistically significant changes between NSD and HSD conditions.

clear all;
close all;
clc;

%% ------------------ Set Folder Paths ------------------
pathway = pwd;

% Load Sampling Data for NSD & HSD
model_sampling_in_name = {'NSD', 'HSD'};
model_sampling_in = cell(1, length(model_sampling_in_name));

for i = 1:length(model_sampling_in_name)
    fileName = fullfile(pwd,'C_sampling', model_sampling_in_name{i});
    fileName = strcat(fileName,'.mat');
    model_sampling_in{i} = load(fileName);
end

%% ------------------ Load NAD(H) Reaction List from Excel ------------------
nadh_file = 'E:\Projects\revision\Code_2_upload\NAD_reaction_analysis\01_nad_balance\NADH_Production_Consumption_Analysis.xlsx';

% Read the reaction IDs
nadh_reactions_table = readtable(nadh_file, 'Sheet', 'All_NADH_Reactions');

% Filter out NAD_demand if present
nadh_reaction_ids = nadh_reactions_table.Reaction_ID;
nadh_reaction_ids = nadh_reaction_ids(~strcmp(nadh_reaction_ids, 'NAD_demand'));

fprintf('Loaded %d NAD(H)-dependent reactions from Excel\n', length(nadh_reaction_ids));

%% ------------------ Create Output Directory ------------------
subfolder = fullfile(pathway, '02_z_test_NAD');

if ~exist(subfolder, 'dir')
    mkdir(subfolder);
end

%% ------------------ Get indices of NAD(H) reactions in sampling data ------------------
rxns_nsd = model_sampling_in{1}.samples.rxns;
rxns_hsd = model_sampling_in{2}.samples.rxns;

% Find indices of NAD(H) reactions
[~, nadh_idx_nsd] = ismember(nadh_reaction_ids, rxns_nsd);
[~, nadh_idx_hsd] = ismember(nadh_reaction_ids, rxns_hsd);

% Remove reactions not found in model (idx = 0)
valid_idx = (nadh_idx_nsd > 0) & (nadh_idx_hsd > 0);
nadh_reaction_ids = nadh_reaction_ids(valid_idx);
nadh_idx_nsd = nadh_idx_nsd(valid_idx);
nadh_idx_hsd = nadh_idx_hsd(valid_idx);

fprintf('Found %d NAD(H) reactions in both NSD and HSD models\n', length(nadh_reaction_ids));

%% ------------------ Add Missing NAD+ Biosynthesis Reactions ------------------
% MAR04261 and MAR04269 are NAD+ biosynthesis reactions that should be included
missing_nad_rxns = {'MAR04261'; 'MAR04269'};

% Check if they exist in the model
for i = 1:length(missing_nad_rxns)
    if ismember(missing_nad_rxns{i}, rxns_nsd)
        if ~ismember(missing_nad_rxns{i}, nadh_reaction_ids)
            fprintf('Adding missing NAD+ reaction: %s\n', missing_nad_rxns{i});
            nadh_reaction_ids = [nadh_reaction_ids; missing_nad_rxns{i}];
        end
    else
        fprintf('Warning: %s not found in model\n', missing_nad_rxns{i});
    end
end

% Recalculate indices after adding missing reactions
[~, nadh_idx_nsd] = ismember(nadh_reaction_ids, rxns_nsd);
[~, nadh_idx_hsd] = ismember(nadh_reaction_ids, rxns_hsd);

% Remove reactions not found in model (idx = 0)
valid_idx = (nadh_idx_nsd > 0) & (nadh_idx_hsd > 0);
nadh_reaction_ids = nadh_reaction_ids(valid_idx);
nadh_idx_nsd = nadh_idx_nsd(valid_idx);
nadh_idx_hsd = nadh_idx_hsd(valid_idx);

fprintf('Total NAD(H) reactions after adding missing ones: %d\n', length(nadh_reaction_ids));

%% ------------------ Z-test Calculation for Each NAD(H) Reaction ------------------
n_reactions = length(nadh_reaction_ids);

% Initialize result arrays
Z_scores = zeros(n_reactions, 1);
P_values = zeros(n_reactions, 1);
log2FC = zeros(n_reactions, 1);
flux_nsd_mean = zeros(n_reactions, 1);
flux_hsd_mean = zeros(n_reactions, 1);
flux_nsd_sd = zeros(n_reactions, 1);
flux_hsd_sd = zeros(n_reactions, 1);
flux_nsd_sem = zeros(n_reactions, 1);
flux_hsd_sem = zeros(n_reactions, 1);
log2FC_sd = zeros(n_reactions, 1);
log2FC_sem = zeros(n_reactions, 1);

% Get flux data
flux_all_nsd = model_sampling_in{1}.samples.points;  % reactions × samples
flux_all_hsd = model_sampling_in{2}.samples.points;

fprintf('Computing Z-scores for %d NAD(H) reactions...\n', n_reactions);

for i = 1:n_reactions
    if mod(i, 50) == 0
        fprintf('Processing reaction %d/%d\n', i, n_reactions);
    end
    
    % Get flux samples for this reaction (take absolute values)
    flux_nsd = abs(flux_all_nsd(nadh_idx_nsd(i), :));
    flux_hsd = abs(flux_all_hsd(nadh_idx_hsd(i), :));
    
    % Calculate statistics
    flux_nsd_mean(i) = mean(flux_nsd, 'omitnan');
    flux_hsd_mean(i) = mean(flux_hsd, 'omitnan');
    
    flux_nsd_sd(i) = std(flux_nsd, 'omitnan');
    flux_hsd_sd(i) = std(flux_hsd, 'omitnan');
    
    n_nsd = sum(~isnan(flux_nsd));
    n_hsd = sum(~isnan(flux_hsd));
    
    flux_nsd_sem(i) = flux_nsd_sd(i) / sqrt(n_nsd);
    flux_hsd_sem(i) = flux_hsd_sd(i) / sqrt(n_hsd);
    
    % Calculate log2FC for each sample, then get SD and SEM
    log2FC_samples = log2((flux_hsd + eps) ./ (flux_nsd + eps));
    log2FC(i) = log2((flux_hsd_mean(i) + eps) / (flux_nsd_mean(i) + eps));
    log2FC_sd(i) = std(log2FC_samples, 'omitnan');
    log2FC_sem(i) = log2FC_sd(i) / sqrt(sum(~isnan(log2FC_samples)));
    
    % Calculate Z-score
    var_nsd = var(flux_nsd, 0, 'omitnan');
    var_hsd = var(flux_hsd, 0, 'omitnan');
    
    den = sqrt(var_hsd/max(n_hsd,1) + var_nsd/max(n_nsd,1));
    z = (flux_hsd_mean(i) - flux_nsd_mean(i)) / max(den, eps);
    
    % Calculate two-tailed p-value
    p = 2 * (1 - normcdf(abs(z)));
    
    Z_scores(i) = z;
    P_values(i) = p;
end

%% ------------------ Adjust P-values (Bonferroni) ------------------
valid_p_idx = ~isnan(P_values);
P_adj = nan(size(P_values));

if sum(valid_p_idx) > 0
    P_adj(valid_p_idx) = P_values(valid_p_idx) * sum(valid_p_idx);  % Simple Bonferroni
    P_adj(P_adj > 1) = 1;  % Cap at 1
    
    % Count significant reactions
    sig_count = sum(P_adj < 0.05, 'omitnan');
    total_tested = sum(~isnan(P_adj));
    fraction_sig = sig_count / total_tested;
    
    fprintf('✅ %d out of %d reactions are significant (Bonferroni-adjusted p < 0.05)\n', sig_count, total_tested);
    fprintf('➡️  %.2f%% of reactions are significant\n', fraction_sig * 100);
else
    fprintf('⚠️  No valid p-values calculated\n');
end

%% ------------------ Get additional annotations ------------------
% Get subsystems and gene rules for NAD(H) reactions
subsystems_nsd = model_sampling_in{1}.samples.subSystems(nadh_idx_nsd);
genes_nsd = model_sampling_in{1}.samples.grRules(nadh_idx_nsd);

%% ------------------ Extract compartment information ------------------
fprintf('Extracting compartment information...\n');

% Get metabolite information
mets = model_sampling_in{1}.samples.mets;
S_matrix = model_sampling_in{1}.samples.S;

% Extract compartments for each NAD(H) reaction
compartments = cell(length(nadh_reaction_ids), 1);

for i = 1:length(nadh_reaction_ids)
    % Find metabolites involved in this reaction
    met_idx = find(abs(S_matrix(:, nadh_idx_nsd(i))) > 0);
    
    % Extract compartment codes from metabolite IDs
    met_ids = mets(met_idx);
    comps = regexp(met_ids, '\[([a-z])\]', 'tokens');
    comps = [comps{:}];
    comps = unique([comps{:}]);
    
    if ~isempty(comps)
        compartments{i} = strjoin(comps, '; ');
    else
        compartments{i} = '';
    end
end

% Clean up subsystems (handle nested cells)
subsystem_names = cell(length(subsystems_nsd), 1);
for i = 1:length(subsystems_nsd)
    if iscell(subsystems_nsd{i})
        subsystem_names{i} = subsystems_nsd{i}{1};
    else
        subsystem_names{i} = subsystems_nsd{i};
    end
end

%% ------------------ Save Full Results ------------------
z_results = table(nadh_reaction_ids, subsystem_names, compartments, genes_nsd, ...
    Z_scores, P_values, P_adj, log2FC, ...
    flux_nsd_mean, flux_hsd_mean, ...
    flux_nsd_sd, flux_hsd_sd, ...
    flux_nsd_sem, flux_hsd_sem, ...
    log2FC_sd, log2FC_sem, ...
    'VariableNames', {'Reaction_ID', 'Subsystem', 'Compartments', 'Gene_Rules', ...
    'Z_score', 'P_value', 'P_adj', 'Log2FoldChange', ...
    'NSD_flux_mean', 'HSD_flux_mean', ...
    'NSD_flux_SD', 'HSD_flux_SD', ...
    'NSD_flux_SEM', 'HSD_flux_SEM', ...
    'Log2FC_SD', 'Log2FC_SEM'});

writetable(z_results, fullfile(subfolder, 'Z_test_NAD_reactions_Bonferroni.xlsx'));

% Sort by absolute Z-score (most significant first)
z_results_sorted = sortrows(z_results, 'Z_score', 'descend');
writetable(z_results_sorted, fullfile(subfolder, 'Z_test_NAD_reactions_sorted_by_Zscore.xlsx'));

% Sort by log2FC
z_results_sorted_fc = sortrows(z_results, 'Log2FoldChange', 'descend');
writetable(z_results_sorted_fc, fullfile(subfolder, 'Z_test_NAD_reactions_sorted_by_log2FC.xlsx'));

% Save only significant reactions
z_results_sig = z_results(P_adj < 0.05, :);
writetable(z_results_sig, fullfile(subfolder, 'Z_test_NAD_reactions_significant_only.xlsx'));

fprintf('✅ NAD(H) reaction Z-test results saved (%d reactions, %d significant).\n', ...
    height(z_results), height(z_results_sig));

%% ------------------ Summary Statistics ------------------
fprintf('\n=== NAD(H) Reaction Z-test Analysis Summary ===\n');
fprintf('Total NAD(H) reactions analyzed: %d\n', n_reactions);
fprintf('Reactions with valid p-values: %d\n', sum(~isnan(P_values)));
if sum(valid_p_idx) > 0
    fprintf('Significant reactions (p_adj < 0.05): %d (%.1f%%)\n', sig_count, fraction_sig * 100);
    fprintf('Reactions with increased flux (log2FC > 0): %d\n', sum(log2FC > 0));
    fprintf('Reactions with decreased flux (log2FC < 0): %d\n', sum(log2FC < 0));
end
fprintf('=============================================\n');

%% ------------------ Create visualization-ready data for R ------------------
% Export in format ready for your R dot plot script
r_ready_data = z_results;
r_ready_data.abs_log2FC = abs(z_results.Log2FoldChange);
r_ready_data.is_significant = P_adj < 0.05;
r_ready_data.significance_label = cell(height(z_results), 1);
r_ready_data.significance_label(P_adj < 0.001) = {'***'};
r_ready_data.significance_label(P_adj >= 0.001 & P_adj < 0.01) = {'**'};
r_ready_data.significance_label(P_adj >= 0.01 & P_adj < 0.05) = {'*'};
r_ready_data.significance_label(P_adj >= 0.05) = {'ns'};

writetable(r_ready_data, fullfile(subfolder, 'Z_test_NAD_reactions_for_R_visualization.xlsx'));

