%% script_05_PFI_z_test
% Purpose: Perform Z-test on sampling-derived fluxes within the pathway to identify statistically significant
% changes between NSD and HSD models across metabolic subsystems.
% Use pathway flux as the average of non zero flux magnitudes within each
% pathway, with equal weights on reactions.

clear all;
close all;
clc;

%% ------------------ Set Folder Paths ------------------
% Define base path and load folder
pathway = pwd;

% Load FVA Data
if ~exist('out_all_fvaBounded.mat','var')
    [parentFolder, ~, ~] = fileparts(pathway);
%     load(fullfile(parentFolder,'3_Flux_analysis_orig','2_FBA_FVA_rxn_bound_updated', 'out_all_fvaBounded.mat'));
end

% Load Sampling Data for NSD & HSD
model_sampling_in_name = {'NSD', 'HSD'};
model_sampling_in = cell(1, length(model_sampling_in_name));

for i = 1:length(model_sampling_in_name)
    fileName = fullfile('models', model_sampling_in_name{i});
    fileName = strcat(fileName,'.mat');
    model_sampling_in{i} = load(fileName);
end

% Load Subsystem Definitions
if ~exist('subsysAll.mat','var')
    load(fullfile(pathway, 'subsysAll.mat'));
end

%% ------------------ Create Output Directory ------------------
subfolder = fullfile(pathway, '5_PFI_zTest_sampling');
if ~exist(subfolder, 'dir')
    mkdir(subfolder);
end

%% ------------------ Z-test Calculation Across Pathways ------------------
Z_scores = zeros(length(subsysAll_fruitflyGEM), 1);
P_values = zeros(length(subsysAll_fruitflyGEM), 1);
log2FC = zeros(length(subsysAll_fruitflyGEM), 1);
pathway_flux_nsd = zeros(length(subsysAll_fruitflyGEM), 1);
pathway_flux_hsd = zeros(length(subsysAll_fruitflyGEM), 1);
rxn_counts = zeros(length(subsysAll_fruitflyGEM), 1);

pathway_flux_nsd_sd = zeros(length(subsysAll_fruitflyGEM), 1);
pathway_flux_hsd_sd = zeros(length(subsysAll_fruitflyGEM), 1);
pathway_flux_nsd_sem = zeros(length(subsysAll_fruitflyGEM), 1);
pathway_flux_hsd_sem = zeros(length(subsysAll_fruitflyGEM), 1);
log2FC_sem = zeros(length(subsysAll_fruitflyGEM), 1);
log2FC_sd = zeros(length(subsysAll_fruitflyGEM), 1);

fprintf('Computing z-scores for %d pathways...\n', length(subsysAll_fruitflyGEM));

for j = 1:length(subsysAll_fruitflyGEM)
    if mod(j, 50) == 0
        fprintf('Processing pathway %d/%d\n', j, length(subsysAll_fruitflyGEM));
    end
    
    flux_nsd = [];
    flux_hsd = [];

    for cond = 1:2
        rxns_all = model_sampling_in{cond}.samples.rxns;
        flux_all = model_sampling_in{cond}.samples.points;
        subsystems = model_sampling_in{cond}.samples.subSystems;

        for k = 1:length(subsystems)
            if iscell(subsystems{k})  % Use {} to access the cell content
                subsystems{k} = subsystems{k}{1};  % Use {} again to dig into the inner cell
            end
        end

        idx_rxns = find(ismember(subsystems, subsysAll_fruitflyGEM{j}));
        
        rxn_fluxes = abs(flux_all(idx_rxns, :));
        
        % Remove reactions that are completely inactive (all zero across samples)
        nonzero_idx = any(rxn_fluxes ~= 0, 2);     % keep rows that have any nonzero flux
        rxn_fluxes = rxn_fluxes(nonzero_idx, :);

        weighted_flux = mean(rxn_fluxes, 1, 'omitnan');

        if cond == 1
            flux_nsd = weighted_flux;
            rxn_counts(j) = size(rxn_fluxes,1);
        else
            flux_hsd = weighted_flux;
        end
    end

    if ~isempty(flux_nsd) && ~isempty(flux_hsd)
        % Store pathway fluxes
        pathway_flux_nsd(j) = mean(flux_nsd, 'omitnan');
        pathway_flux_hsd(j) = mean(flux_hsd, 'omitnan');

        % CALCULATE SD and SEM:
        pathway_flux_nsd_sd(j) = std(flux_nsd, 'omitnan');
        pathway_flux_hsd_sd(j) = std(flux_hsd, 'omitnan');
        pathway_flux_nsd_sem(j) = std(flux_nsd, 'omitnan') / sqrt(sum(~isnan(flux_nsd)));
        pathway_flux_hsd_sem(j) = std(flux_hsd, 'omitnan') / sqrt(sum(~isnan(flux_hsd)));

        % Calculate log2FC for each sample point, then get SEM
        log2FC_samples = log2((flux_hsd + eps) ./ (flux_nsd + eps));
        log2FC_sd(j) = std(log2FC_samples, 'omitnan');
        log2FC_sem(j) = std(log2FC_samples, 'omitnan') / sqrt(sum(~isnan(log2FC_samples)));
      
        mean_diff = mean(flux_hsd) - mean(flux_nsd);
        
        % Correct Z-score calculation
         nN = numel(flux_nsd);
         nH = numel(flux_hsd);

%         nN = length(idx_rxns);
%         nH = length(idx_rxns);
        
        vN = var(flux_nsd, 0, 'omitnan');
        vH = var(flux_hsd, 0, 'omitnan');
        den = sqrt(vH/max(nH,1) + vN/max(nN,1));
        z = (mean(flux_hsd,'omitnan') - mean(flux_nsd,'omitnan')) / max(den, eps);

        p = 2 * (1 - normcdf(abs(z)));
        
        Z_scores(j) = z;
        P_values(j) = p;
        log2FC(j) = log2((mean(flux_hsd) + eps) / (mean(flux_nsd) + eps));
    else
        % No valid flux data
        Z_scores(j) = NaN;
        P_values(j) = NaN;
        log2FC(j) = NaN;
        pathway_flux_nsd(j) = NaN;
        pathway_flux_hsd(j) = NaN;
        rxn_counts(j) = 0;
        pathway_flux_nsd_sd(j) = NaN;
        pathway_flux_hsd_sd(j) = NaN;
        pathway_flux_nsd_sem(j) = NaN;
        pathway_flux_hsd_sem(j) = NaN;
        log2FC_sd(j) = NaN;
        log2FC_sem(j) = NaN;
    end
end

%% ------------------ Adjust P-values (Bonferroni) ------------------
valid_p_idx = ~isnan(P_values);
P_adj = nan(size(P_values));

if sum(valid_p_idx) > 0
    P_adj(valid_p_idx) = P_values(valid_p_idx) * sum(valid_p_idx);  % Simple Bonferroni
    P_adj(P_adj > 1) = 1;  % Cap at 1
    
    % Count significant pathways
    sig_count = sum(P_adj < 0.05, 'omitnan');
    total_tested = sum(~isnan(P_adj));
    fraction_sig = sig_count / total_tested;
    
    fprintf('✅ %d out of %d pathways are significant (Bonferroni-adjusted p < 0.05)\n', sig_count, total_tested);
else
    fprintf('⚠️  No valid p-values calculated\n');
end

%% ------------------ Report Non-Significant Pathways ------------------
nonsig_idx = find(P_adj >= 0.05 & ~isnan(P_adj));
nonsig_pathways = subsysAll_fruitflyGEM(nonsig_idx);

% Save non-significant pathways to TXT
nonsig_file = fullfile(subfolder, 'NonSignificant_Pathways_sampling_Bonferroni_0.05.txt');
fid = fopen(nonsig_file, 'w');
fprintf(fid, 'Non-significant pathways from sampling analysis (Bonferroni-adjusted p ≥ 0.05):\n\n');
for i = 1:length(nonsig_pathways)
    fprintf(fid, '%s\n', nonsig_pathways{i});
end
fclose(fid);

%% ------------------ Save Full Results ------------------
z_results = table(subsysAll_fruitflyGEM, Z_scores, P_values, P_adj, log2FC, ...
    pathway_flux_nsd, pathway_flux_hsd, rxn_counts, ...
    pathway_flux_nsd_sd, pathway_flux_hsd_sd, ...
    pathway_flux_nsd_sem, pathway_flux_hsd_sem, log2FC_sd, log2FC_sem, ...
    'VariableNames', {'Pathway', 'Z_score', 'P_value', 'P_adj', 'Log2FoldChange_flux', ...
    'NSD_flux', 'HSD_flux', 'Reaction_count', ...
    'NSD_flux_SD', 'HSD_flux_SD', 'NSD_flux_SEM', 'HSD_flux_SEM', 'Log2FC_SD', 'Log2FC_SEM'});

writetable(z_results, fullfile(subfolder, 'Z_test_sampling_results_Bonferroni.xlsx'));

% Sort by log2FC (descending: most up-regulated pathways first)
z_results_sorted = sortrows(z_results, 'Log2FoldChange_flux', 'descend');
writetable(z_results_sorted, fullfile(subfolder, 'Z_test_sampling_results_Bonferroni_sorted.xlsx'));

% Remove rows with any NaN values and sort
z_results_clean = rmmissing(z_results);
z_results_sorted_cleaned = sortrows(z_results_clean, 'Log2FoldChange_flux', 'descend');

% Save the cleaned tables
writetable(z_results_clean, fullfile(subfolder, 'Z_test_sampling_results_Bonferroni_cleaned.xlsx'));
writetable(z_results_sorted_cleaned, fullfile(subfolder, 'Z_test_sampling_results_Bonferroni_sorted_cleaned.xlsx'));

fprintf('✅ Cleaned sampling Z-test results saved without NaN rows (%d pathways).\n', height(z_results_clean));
