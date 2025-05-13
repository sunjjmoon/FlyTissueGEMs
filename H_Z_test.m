%% =================== Z-Test on Pathway Flux Index (PFI) ===================
% Purpose: Perform Z-test on sampling-derived PFIs to identify statistically significant
% changes between NSD and HSD models across metabolic subsystems.
%
% This version replaces PFVI calculations with Z-test-based significance testing
% between two sampled conditions (NSD and HSD) using sampling-derived flux values.
% ===========================================================================

clear all;
close all;
clc;

%% ------------------ Set Folder Paths ------------------
% Define base path and load folder
pathway = pwd;
loadfolder = fullfile(pathway, 'B_FBA_FVA_rxn_bound_updated');

% Load FVA Data
if ~exist('out_all_fvaBounded.mat','var')
    load(fullfile(loadfolder, 'out_all_fvaBounded.mat'));
end

% Load Sampling Data for NSD & HSD
model_sampling_in_name = {'NSD', 'HSD'};
model_sampling_in = cell(1, length(model_sampling_in_name));

for i = 1:length(model_sampling_in_name)
    fileName = fullfile(pathway, 'C_sampling', model_sampling_in_name{i});
    model_sampling_in{i} = load(fileName);
end

% Load Subsystem Definitions
if ~exist('subsysAll_fruitflyGEM.mat','var')
    load(fullfile(pathway, 'subsysAll_fruitflyGEM.mat'));
end

%% ------------------ Create Output Directory ------------------
subfolder = fullfile(pathway, 'H_z_test');
if ~exist(subfolder, 'dir')
    mkdir(subfolder);
end

%% ------------------ Z-test Calculation Across Pathways ------------------
Z_scores = zeros(length(subsysAll_fruitflyGEM), 1);
P_values = zeros(length(subsysAll_fruitflyGEM), 1);
log2FC = zeros(length(subsysAll_fruitflyGEM), 1);

for j = 1:length(subsysAll_fruitflyGEM)
    flux_nsd = [];
    flux_hsd = [];

    for cond = 1:2
        rxns_all = model_sampling_in{cond}.samples.rxns;
        flux_all = model_sampling_in{cond}.samples.points;
        subsystems = out_all_fvaBounded{cond}.subSystems;

        for k = 1:length(subsystems)
            if iscell(subsystems{k})  % Use {} to access the cell content
                subsystems{k} = subsystems{k}{1};  % Use {} again to dig into the inner cell
            end
        end



        idx_rxns = find(ismember(subsystems, subsysAll_fruitflyGEM{j}));
        rxn_fluxes = abs(flux_all(idx_rxns, :));
        
        weighted_flux = mean(rxn_fluxes, 1, 'omitnan');

        if cond == 1
            flux_nsd = weighted_flux;
        else
            flux_hsd = weighted_flux;
        end
    end

    if ~isempty(flux_nsd) && ~isempty(flux_hsd)
        mean_diff = mean(flux_hsd) - mean(flux_nsd);
        pooled_sd = sqrt(var(flux_nsd)/length(1) + var(flux_hsd)/length(1));
        z = mean_diff / pooled_sd;
        p = 2 * (1 - normcdf(abs(z)));
        
        Z_scores(j) = z;
        P_values(j) = p;
        log2FC(j) = log2(mean(flux_hsd) / mean(flux_nsd));
    end
end

% Save results
z_results = table(subsysAll_fruitflyGEM, Z_scores, P_values, log2FC, ...
    'VariableNames', {'Pathway', 'Z_score', 'P_value', 'Log2FoldChange'});

writetable(z_results, fullfile(subfolder, 'Z_test_results.xlsx'));
fprintf('✅ Z-test results saved to "Z_test_results.xlsx".\n');
