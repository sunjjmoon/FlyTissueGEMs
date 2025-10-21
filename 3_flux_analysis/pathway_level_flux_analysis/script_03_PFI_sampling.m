%% script_03_PFI_sampling

clc;
close all;

%% 1. Set up
pathway = pwd;
results_dir = fullfile(pathway, '3_sampling_PFI');

if ~exist(results_dir, 'dir')
    mkdir(results_dir); 
end

%% 2. Load sampling data directly (NOT from Script 06)
fprintf('=== LOADING SAMPLING DATA DIRECTLY ===\n');

model_sampling_in_name = {'NSD', 'HSD'};
model_sampling_in = cell(1, length(model_sampling_in_name));

for i = 1:length(model_sampling_in_name)
    sampling_file = fullfile(pwd, [model_sampling_in_name{i} '.mat']);
    model_sampling_in{i} = load(sampling_file);
    fprintf('Loaded sampling data for %s\n', model_sampling_in_name{i});
end

%% 3. Extract all unique pathways
all_pathways = {};
for i = 1:length(model_sampling_in)
    subsystems = model_sampling_in{i}.samples.subSystems;
    % Clean subsystems
    for k = 1:length(subsystems)
        if iscell(subsystems{k})
            subsystems{k} = subsystems{k}{1};
        end
    end
    unique_pathways = unique(subsystems);
    unique_pathways = unique_pathways(~cellfun(@isempty, unique_pathways));
    all_pathways = [all_pathways; unique_pathways];
end
all_pathways = unique(all_pathways);
fprintf('Found %d unique pathways\n', length(all_pathways));

%% 4. Calculate pathway-level statistics using sampling distributions
fprintf('\n=== CALCULATING PATHWAY STATISTICS FROM SAMPLING ===\n');

% Initialize storage
nsd_pathway_flux_samples = cell(length(all_pathways), 1);
hsd_pathway_flux_samples = cell(length(all_pathways), 1);

for p = 1:length(all_pathways)
    pathway_name = all_pathways{p};
    
    % Process NSD
    subsystems_nsd = model_sampling_in{1}.samples.subSystems;
    for k = 1:length(subsystems_nsd)
        if iscell(subsystems_nsd{k})
            subsystems_nsd{k} = subsystems_nsd{k}{1};
        end
    end
    idx_rxns_nsd = find(ismember(subsystems_nsd, pathway_name));
    
    if ~isempty(idx_rxns_nsd)
        flux_all_nsd = model_sampling_in{1}.samples.points;
        rxn_fluxes_nsd = abs(flux_all_nsd(idx_rxns_nsd, :));
        % Calculate pathway flux for EACH sample (mean across reactions)
        nsd_pathway_flux_samples{p} = mean(rxn_fluxes_nsd, 1, 'omitnan');
    else
        nsd_pathway_flux_samples{p} = [];
    end
    
    % Process HSD
    subsystems_hsd = model_sampling_in{2}.samples.subSystems;
    for k = 1:length(subsystems_hsd)
        if iscell(subsystems_hsd{k})
            subsystems_hsd{k} = subsystems_hsd{k}{1};
        end
    end
    idx_rxns_hsd = find(ismember(subsystems_hsd, pathway_name));
    
    if ~isempty(idx_rxns_hsd)
        flux_all_hsd = model_sampling_in{2}.samples.points;
        rxn_fluxes_hsd = abs(flux_all_hsd(idx_rxns_hsd, :));
        % Calculate pathway flux for EACH sample (mean across reactions)
        hsd_pathway_flux_samples{p} = mean(rxn_fluxes_hsd, 1, 'omitnan');
    else
        hsd_pathway_flux_samples{p} = [];
    end
end

%% 5. Calculate summary statistics from distributions
mean_nsd = zeros(length(all_pathways), 1);
mean_hsd = zeros(length(all_pathways), 1);
median_nsd = zeros(length(all_pathways), 1);
median_hsd = zeros(length(all_pathways), 1);
std_nsd = zeros(length(all_pathways), 1);
std_hsd = zeros(length(all_pathways), 1);

for p = 1:length(all_pathways)
    if ~isempty(nsd_pathway_flux_samples{p})
        mean_nsd(p) = mean(nsd_pathway_flux_samples{p}, 'omitnan');
        median_nsd(p) = median(nsd_pathway_flux_samples{p}, 'omitnan');
        std_nsd(p) = std(nsd_pathway_flux_samples{p}, 'omitnan');
    else
        mean_nsd(p) = NaN;
        median_nsd(p) = NaN;
        std_nsd(p) = NaN;
    end
    
    if ~isempty(hsd_pathway_flux_samples{p})
        mean_hsd(p) = mean(hsd_pathway_flux_samples{p}, 'omitnan');
        median_hsd(p) = median(hsd_pathway_flux_samples{p}, 'omitnan');
        std_hsd(p) = std(hsd_pathway_flux_samples{p}, 'omitnan');
    else
        mean_hsd(p) = NaN;
        median_hsd(p) = NaN;
        std_hsd(p) = NaN;
    end
end

% Calculate ratios and log2FC
mean_ratio = mean_hsd ./ mean_nsd;
mean_log2fc = log2(mean_ratio);

median_ratio = median_hsd ./ median_nsd;
median_log2fc = log2(median_ratio);

%% 6. Save results
excel_filename = fullfile(results_dir, 'pathway_flux_analysis.xlsx');

% Mean comparison
mean_comparison = table(all_pathways, mean_nsd, mean_hsd, mean_ratio, mean_log2fc, std_nsd, std_hsd, ...
    'VariableNames', {'Pathway', 'NSD', 'HSD', 'HSD_NSD_Ratio', 'Log2FC', 'NSD_SD', 'HSD_SD'});

valid_idx = ~isnan(mean_log2fc) & ~isinf(mean_log2fc);
mean_comparison_sorted = sortrows(mean_comparison(valid_idx, :), 'Log2FC', 'descend');

writetable(mean_comparison_sorted, excel_filename, 'Sheet', 'Mean_HSD_vs_NSD', 'WriteRowNames', false);

fprintf('✅ Script 07 completed with sampling-based pathway statistics\n');
fprintf('   This now matches the Z-test methodology!\n');