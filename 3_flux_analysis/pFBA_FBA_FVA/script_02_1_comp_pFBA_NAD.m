%% script_02_1_comp_pFBA_NAD

clc;
% clear all;

%% 1. Set up
current_path = pwd;
results_dir = fullfile(current_path, '02_1_comp_nad');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

%% 2. Load
% load comp results
load(fullfile(current_path,'02_comp_fba_pFBA','FBA_pFBA_comparison.mat'));

% load muscle gems
file_loc = current_path;
load(strcat(file_loc,'\models\model_out_cbra_u.mat'));
models = model_out_cbra_u;

nad_analysis = cell(length(comparison_results), 1);


%% 3. Define NAD/NADH/NADPH metabolite IDs (adjust these based on your model)
nad_metabolites = {
    'MAM02552c[c]',  % NAD cytosol
    'MAM02552m[m]',  % NAD mitochondria
    'MAM02553c[c]',  % NADH cytosol  
    'MAM02553m[m]',  % NADH mitochondria
    'MAM02554c[c]',  % NADP cytosol
    'MAM02554m[m]',  % NADP mitochondria
    'MAM02555c[c]',  % NADPH cytosol
    'MAM02555m[m]'   % NADPH mitochondria
};


%% Process each model
for i = 1:length(comparison_results)
    model_id = comparison_results{i}.modelID;
    fprintf('Analyzing NAD reactions for model %s\n', model_id);
    
    % Get the original model
    model = models{i,1};
    
    % Find all reactions involving NAD/NADH/NADPH
    nad_rxn_indices = [];
    nad_rxn_names = {};
    
    for j = 1:length(nad_metabolites)
        met_idx = find(strcmp(model.mets, nad_metabolites{j}));
        if ~isempty(met_idx)
            % Find reactions involving this metabolite
            rxn_idx = find(model.S(met_idx, :) ~= 0);
            nad_rxn_indices = [nad_rxn_indices, rxn_idx];
            for k = 1:length(rxn_idx)
                nad_rxn_names{end+1} = model.rxns{rxn_idx(k)};
            end
        end
    end
    
    % Remove duplicates
    nad_rxn_indices = unique(nad_rxn_indices);
    nad_rxn_names = unique(nad_rxn_names);
    
    fprintf('  Found %d NAD-related reactions\n', length(nad_rxn_indices));
    
    % Get fluxes from FBA and pFBA solutions
    fba_solution = comparison_results{i}.FBA.solution;
    pfba_solution = comparison_results{i}.pFBA.solution;
    
    if fba_solution.stat == 1 && pfba_solution.stat == 1
        
        % Extract NAD-related fluxes
        fba_nad_fluxes = fba_solution.x(nad_rxn_indices);
        pfba_nad_fluxes = pfba_solution.x(nad_rxn_indices);
        
        % Calculate total NAD-related flux
        fba_nad_total = sum(abs(fba_nad_fluxes));
        pfba_nad_total = sum(abs(pfba_nad_fluxes));
        
        fba_nad_total_avg = mean(abs(fba_nad_fluxes));
        pfba_nad_total_avg = mean(abs(pfba_nad_fluxes));
        
        % Calculate reduction
        nad_flux_reduction = ((fba_nad_total - pfba_nad_total) / fba_nad_total) * 100;
        nad_flux_reduction_avg = ((fba_nad_total_avg - pfba_nad_total_avg) / fba_nad_total_avg) * 100;
        
        % Find high-flux NAD reactions
        high_flux_threshold = 100; % 10 % to Vmax
        fba_high_nad = find(abs(fba_nad_fluxes) > high_flux_threshold);
        pfba_high_nad = find(abs(pfba_nad_fluxes) > high_flux_threshold);
        
        % Identify potential cycling reactions
        cycling_candidates = [];
        cycling_details = {};
        for j = 1:length(fba_nad_fluxes)
            if abs(fba_nad_fluxes(j)) > 100 && ... % High flux in FBA. 10% to  Vmax
               abs(fba_nad_fluxes(j)) > 2 * abs(pfba_nad_fluxes(j)) && ... % 100 % greater than the pFBA result
               abs(pfba_nad_fluxes(j)) < abs(fba_nad_fluxes(j)) % Actually reduced
                cycling_candidates(end+1) = j;
                cycling_details{end+1} = sprintf('%s: FBA=%.1f → pFBA=%.1f (%.1f%% reduction)', ...
                    nad_rxn_names{j}, fba_nad_fluxes(j), pfba_nad_fluxes(j), ...
                    ((abs(fba_nad_fluxes(j)) - abs(pfba_nad_fluxes(j))) / abs(fba_nad_fluxes(j))) * 100);
            end
        end
        
        % Store results
        nad_analysis{i}.modelID = model_id;
        nad_analysis{i}.nad_rxn_names = nad_rxn_names;
        nad_analysis{i}.nad_rxn_indices = nad_rxn_indices;
        nad_analysis{i}.fba_nad_fluxes = fba_nad_fluxes;
        nad_analysis{i}.pfba_nad_fluxes = pfba_nad_fluxes;
        nad_analysis{i}.fba_nad_total = fba_nad_total;
        nad_analysis{i}.pfba_nad_total = pfba_nad_total;
        
        nad_analysis{i}.fba_nad_total_avg = fba_nad_total_avg;
        nad_analysis{i}.pfba_nad_total_avg = pfba_nad_total_avg;
        
        nad_analysis{i}.nad_flux_reduction = nad_flux_reduction;
        nad_analysis{i}.nad_flux_reduction_avg = nad_flux_reduction_avg;

        nad_analysis{i}.fba_high_flux_rxns = nad_rxn_names(fba_high_nad);
        nad_analysis{i}.pfba_high_flux_rxns = nad_rxn_names(pfba_high_nad);
        nad_analysis{i}.fba_high_flux_values = fba_nad_fluxes(fba_high_nad);
        nad_analysis{i}.pfba_high_flux_values = pfba_nad_fluxes(pfba_high_nad);
        nad_analysis{i}.cycling_candidates = cycling_candidates;
        nad_analysis{i}.cycling_details = cycling_details;
        
        % Display results
        fprintf('  NAD flux analysis:\n');
        fprintf('    FBA total NAD flux: %.1f\n', fba_nad_total);
        fprintf('    pFBA total NAD flux: %.1f\n', pfba_nad_total);
        fprintf('    NAD flux reduction: %.1f%%\n', nad_flux_reduction);
        fprintf('    FBA high-flux NAD reactions: %d\n', length(fba_high_nad));
        fprintf('    pFBA high-flux NAD reactions: %d\n', length(pfba_high_nad));
        fprintf('    Potential cycling reactions: %d\n', length(cycling_candidates));
    else
        fprintf('  Warning: Solution failed for %s\n', model_id);
        nad_analysis{i}.modelID = model_id;
        nad_analysis{i}.error = 'Solution failed';
    end
    fprintf('\n');
end

%% Comparison report
for i = 1:length(nad_analysis)
    if isfield(nad_analysis{i}, 'error')
        continue;
    end
    
    res = nad_analysis{i};
    fprintf('Model: %s\n', res.modelID);
    fprintf('----------------------------------------\n');
    
    % Overall flux comparison
    total_fba = comparison_results{i}.FBA.total_flux;
    total_pfba = comparison_results{i}.pFBA.total_flux;
    nad_fba_percent = (res.fba_nad_total / total_fba) * 100;
    nad_pfba_percent = (res.pfba_nad_total / total_pfba) * 100;
    
    fprintf('Total flux - FBA: %.0f, pFBA: %.0f\n', total_fba, total_pfba);
    fprintf('NAD flux - FBA: %.0f (%.1f%%), pFBA: %.0f (%.1f%%)\n', ...
        res.fba_nad_total, nad_fba_percent, res.pfba_nad_total, nad_pfba_percent);
    fprintf('NAD flux reduction: %.1f%%\n', res.nad_flux_reduction);
    
    % Show top NAD reactions with largest flux changes
%     fprintf('\nTop NAD reactions with largest flux changes:\n');
    flux_changes = abs(res.fba_nad_fluxes) - abs(res.pfba_nad_fluxes);
    [~, sorted_idx] = sort(flux_changes, 'descend');
    
    for j = 1:min(10, length(sorted_idx)) % Show top 10
        rxn_idx = sorted_idx(j);
        if flux_changes(rxn_idx) > 1 % Only show significant changes
%             fprintf('  %s: FBA=%.1f, pFBA=%.1f, Change=%.1f\n', ...
%                 res.nad_rxn_names{rxn_idx}, res.fba_nad_fluxes(rxn_idx), ...
%                 res.pfba_nad_fluxes(rxn_idx), flux_changes(rxn_idx));
        end
    end
    fprintf('\n');
end

%% Summary table
fprintf('=== NAD FLUX SUMMARY ===\n');
for i = 1:length(nad_analysis)
    if isfield(nad_analysis{i}, 'error')
        continue;
    end
    
    res = nad_analysis{i};
    total_fba = comparison_results{i}.FBA.total_flux;
    total_pfba = comparison_results{i}.pFBA.total_flux;
    nad_fba_percent = (res.fba_nad_total / total_fba) * 100;
    nad_pfba_percent = (res.pfba_nad_total / total_pfba) * 100;
    
    fprintf('%s\t%.0f\t%.0f\t%.0f\t%.0f\t%.1f%%\t%.1f%%\t%.1f%%\n', ...
        res.modelID, total_fba, total_pfba, res.fba_nad_total, res.pfba_nad_total, ...
        res.nad_flux_reduction, nad_fba_percent, nad_pfba_percent);
end

%% NADH cycling analysis
fprintf('\n=== ENHANCED NADH CYCLING ANALYSIS ===\n');
for i = 1:length(nad_analysis)
    if isfield(nad_analysis{i}, 'error')
        continue;
    end
    
    res = nad_analysis{i};
    fprintf('Model: %s\n', res.modelID);
    
    if ~isempty(res.cycling_candidates)
        fprintf('  Potential cycling reactions (high FBA flux, significantly reduced in pFBA):\n');
        for j = 1:length(res.cycling_details)
            fprintf('    %s\n', res.cycling_details{j});
        end
        
        % Calculate what fraction of NAD flux reduction comes from these candidates
        cycling_flux_reduction = 0;
        for j = 1:length(res.cycling_candidates)
            idx = res.cycling_candidates(j);
            cycling_flux_reduction = cycling_flux_reduction + ...
                (abs(res.fba_nad_fluxes(idx)) - abs(res.pfba_nad_fluxes(idx)));
        end
        
        total_nad_reduction = res.fba_nad_total - res.pfba_nad_total;
        cycling_contribution = (cycling_flux_reduction / total_nad_reduction) * 100;
        
        fprintf('  Cycling reactions account for %.1f%% of NAD flux reduction\n', cycling_contribution);
    else
        fprintf('  No obvious cycling reactions identified with current criteria\n');
    end
    fprintf('\n');
end

%% Metabolite-specific analysis
fprintf('=== METABOLITE-SPECIFIC ANALYSIS ===\n');
for i = 1:length(nad_analysis)
    if isfield(nad_analysis{i}, 'error')
        continue;
    end
    
    model = models{i,1};
    res = nad_analysis{i};
    fprintf('Model: %s\n', res.modelID);
    
    % Analyze flux through each NAD metabolite
    for met_idx = 1:length(nad_metabolites)
        met_id = nad_metabolites{met_idx};
        met_model_idx = find(strcmp(model.mets, met_id));
        
        if ~isempty(met_model_idx)
            % Find reactions involving this specific metabolite
            rxn_indices = find(model.S(met_model_idx, :) ~= 0);
            
            if ~isempty(rxn_indices)
                % Calculate total flux through this metabolite
                fba_met_flux = sum(abs(comparison_results{i}.FBA.solution.x(rxn_indices)));
                pfba_met_flux = sum(abs(comparison_results{i}.pFBA.solution.x(rxn_indices)));
                
                if fba_met_flux > 10 % Only report if significant flux
                    reduction = ((fba_met_flux - pfba_met_flux) / fba_met_flux) * 100;
                    fprintf('  %s: FBA=%.1f, pFBA=%.1f (%.1f%% reduction)\n', ...
                        met_id, fba_met_flux, pfba_met_flux, reduction);
                end
            end
        end
    end
    fprintf('\n');
end

%% Statistical summary
fprintf('=== STATISTICAL SUMMARY ===\n');
all_nad_reductions = [];
all_cycling_counts = [];

for i = 1:length(nad_analysis)
    if ~isfield(nad_analysis{i}, 'error')
        all_nad_reductions(end+1) = nad_analysis{i}.nad_flux_reduction;
        all_cycling_counts(end+1) = length(nad_analysis{i}.cycling_candidates);
    end
end

if ~isempty(all_nad_reductions)
    fprintf('NAD flux reduction statistics:\n');
    fprintf('  Mean: %.1f%%, Std: %.1f%%\n', mean(all_nad_reductions), std(all_nad_reductions));
    fprintf('  Range: %.1f%% - %.1f%%\n', min(all_nad_reductions), max(all_nad_reductions));
    fprintf('Average potential cycling reactions per model: %.1f\n', mean(all_cycling_counts));
end

%% Save results to results folder
fprintf('Saving results to %s folder...\n', results_dir);
save(fullfile(results_dir, 'NAD_flux_analysis.mat'), 'nad_analysis');

% Also save a summary table as CSV for easy viewing
summary_filename = fullfile(results_dir, 'NAD_summary.csv');
fid = fopen(summary_filename, 'w');
fprintf(fid, 'Model,Total_FBA,Total_pFBA,NAD_FBA,NAD_pFBA,NAD_Reduction_Percent,NAD_Percent_of_Total_FBA,NAD_Percent_of_Total_pFBA,Cycling_Reactions\n');

for i = 1:length(nad_analysis)
    if ~isfield(nad_analysis{i}, 'error')
        res = nad_analysis{i};
        total_fba = comparison_results{i}.FBA.total_flux;
        total_pfba = comparison_results{i}.pFBA.total_flux;
        nad_fba_percent = (res.fba_nad_total / total_fba) * 100;
        nad_pfba_percent = (res.pfba_nad_total / total_pfba) * 100;
        
        fprintf(fid, '%s,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%.1f,%d\n', ...
            res.modelID, total_fba, total_pfba, res.fba_nad_total, res.pfba_nad_total, ...
            res.nad_flux_reduction, nad_fba_percent, nad_pfba_percent, ...
            length(res.cycling_candidates));
    end
end
fclose(fid);
