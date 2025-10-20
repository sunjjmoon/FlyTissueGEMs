%% Script 02: FBA vs  pFBA Comparison
% Compare regular FBA and proper pFBA using existing results

clc;
% clear all;

% initCobraToolbox
% changeCobraSolver('gurobi', 'LP');

%% 1. Set up
pathway = pwd;
results_dir = fullfile(pathway, '02_comp_fba_pFBA');

if ~exist(results_dir, 'dir')
    mkdir(results_dir); 
end

%% 2. Load models
% Load  pFBA results
pathway = pwd;
load_dir = fullfile(pathway, '01_results');
load(fullfile(load_dir,'pFBA_results.mat'));

file_loc = pathway;
load(strcat(file_loc,'\models\model_out_cbra_u.mat'));
models = model_out_cbra_u;

% Initialize results for comparison
comparison_results = cell(length(models), 1);

%% 3. Run FBA and compare with pFBA
for i = 1:length(models)
    model = models{i,1};
    model_id = model.modelID;
    
    fprintf('Processing model %s (%d/%d)\n', model_id, i, length(models));
    
    model_constrained = apply_constraints(model, model_id);  
    FBA_sol = optimizeCbModel(model_constrained, 'max');
    
    pFBA_sol = pFBA_results{i}.solution;
    
    % Calculate 
    if FBA_sol.stat == 1
        FBA_total_flux = sum(abs(FBA_sol.x));
        FBA_total_flux_avg = mean(abs(FBA_sol.x));
        FBA_objective = FBA_sol.f;
    else
        FBA_total_flux = NaN;
        FBA_objective = NaN;
    end
    
    if pFBA_sol.stat == 1
        pFBA_total_flux = pFBA_results{i}.total_flux;
        pFBA_total_flux_avg = pFBA_results{i}.avg_flux;
        pFBA_objective = pFBA_results{i}.objective;
    else
        pFBA_total_flux = NaN;
        pFBA_objective = NaN;
    end
    
    % Calculate flux reduction and objective change
    if ~isnan(FBA_total_flux) && ~isnan(pFBA_total_flux)
        flux_reduction_percent = ((FBA_total_flux - pFBA_total_flux) / FBA_total_flux) * 100;
        flux_reduction_percent_avg = ((FBA_total_flux_avg - pFBA_total_flux_avg) / FBA_total_flux_avg) * 100;

        objective_change_percent = ((pFBA_objective - FBA_objective) / abs(FBA_objective)) * 100;
    else
        flux_reduction_percent = NaN;
        objective_change_percent = NaN;
    end
    
    % Store results
    comparison_results{i,1}.modelID = model_id;
    comparison_results{i,1}.FBA.solution = FBA_sol;
    comparison_results{i,1}.FBA.objective = FBA_objective;
    comparison_results{i,1}.FBA.total_flux = FBA_total_flux;
    comparison_results{i,1}.FBA.total_flux_avg = FBA_total_flux_avg;
    
    comparison_results{i,1}.pFBA.solution = pFBA_sol;
    comparison_results{i,1}.pFBA.objective = pFBA_objective;
    comparison_results{i,1}.pFBA.total_flux = pFBA_total_flux;
    comparison_results{i,1}.pFBA.total_flux_avg = pFBA_total_flux_avg;
    
    comparison_results{i,1}.flux_reduction_percent = flux_reduction_percent;
    comparison_results{i,1}.flux_reduction_percent_avg = flux_reduction_percent_avg;

    comparison_results{i,1}.objective_change_percent = objective_change_percent;
    
    % Display comparison
    fprintf('  Results for %s:\n', model_id);
    fprintf('    FBA  - Objective: %.4f, Total flux: %.1f\n', FBA_objective, FBA_total_flux);
    fprintf('    pFBA - Objective: %.4f, Total flux: %.1f\n', pFBA_objective, pFBA_total_flux);
    fprintf('    Flux reduction: %.1f%%, Objective change: %.3f%%\n', ...
        flux_reduction_percent, objective_change_percent);
    fprintf('\n');
end

%% Additional Analysis: Flux Distribution Comparison
fprintf('\n=== FLUX DISTRIBUTION ANALYSIS ===\n');
for i = 1:length(comparison_results)
    res = comparison_results{i};
    fprintf('Model: %s\n', res.modelID);
    
    if res.FBA.solution.stat == 1 && res.pFBA.solution.stat == 1
        % Calculate flux statistics
        fba_fluxes = abs(res.FBA.solution.x);
        pfba_fluxes = abs(res.pFBA.solution.x);
        
        % Count active reactions (flux > 1e-6)
        fba_active = sum(fba_fluxes > 1e-6);
        pfba_active = sum(pfba_fluxes > 1e-6);
        
        % High flux reactions (flux > 10)
        fba_high = sum(fba_fluxes > 10);
        pfba_high = sum(pfba_fluxes > 10);
        
        % Very high flux reactions (flux > 100)
        fba_very_high = sum(fba_fluxes > 100);
        pfba_very_high = sum(pfba_fluxes > 100);
        
        fprintf('  Active reactions (>1e-6): FBA=%d, pFBA=%d\n', fba_active, pfba_active);
        fprintf('  High flux reactions (>10): FBA=%d, pFBA=%d\n', fba_high, pfba_high);
        fprintf('  Very high flux (>100): FBA=%d, pFBA=%d\n', fba_very_high, pfba_very_high);
        
        % Calculate flux correlation
%         correlation = corr(fba_fluxes, pfba_fluxes);
%         fprintf('  Flux correlation: %.3f\n', correlation);
    end
    fprintf('\n');
end

%% Identify Major Flux Changes
fprintf('=== REACTIONS WITH LARGEST FLUX CHANGES ===\n');
% load('model_out_cbra_u.mat'); % Reload to get reaction names

file_loc = 'C:\Users\Sun Jin\Documents\revision\tissueGEM_revision\v2\code_to_upload\3_Flux_analysis\models';
load(strcat(file_loc,'\model_out_cbra_u.mat'));

for i = 1:length(comparison_results)
    res = comparison_results{i};
    model = models{i,1};
    fprintf('Model: %s\n', res.modelID);
    
    if res.FBA.solution.stat == 1 && res.pFBA.solution.stat == 1
        % Calculate absolute flux changes
        flux_changes = abs(res.FBA.solution.x) - abs(res.pFBA.solution.x);
        [sorted_changes, sorted_idx] = sort(flux_changes, 'descend');
        
        % Show top 10 reactions with largest flux reductions
        fprintf('  Top reactions with largest flux reductions:\n');
        for j = 1:min(10, length(sorted_idx))
            if sorted_changes(j) > 1 % Only show significant changes
                rxn_idx = sorted_idx(j);
                fprintf('    %s: FBA=%.1f, pFBA=%.1f, Change=%.1f\n', ...
                    model.rxns{rxn_idx}, res.FBA.solution.x(rxn_idx), ...
                    res.pFBA.solution.x(rxn_idx), sorted_changes(j));
            end
        end
    end
    fprintf('\n');
end

%% Save results
save(fullfile(results_dir, 'FBA_pFBA_comparison.mat'), 'comparison_results');

%% Helper function (same as script_1_v2)
function model = apply_constraints(model, model_id)
    reductionFactor = 0.5;
    scaling_factor = 1;
    
    model = changeRxnBounds(model, {'MAR03802','MAR03804'}, 0, 'u');
    model = changeRxnBounds(model, {'MAR03958','MAR00710','MAR04111','MAR04112'}, -1000, 'l');
    model = changeRxnBounds(model, {'MAR02358'}, -1000, 'l');
    model = changeRxnBounds(model, {'MAR02358'}, 0, 'u');
    model = changeRxnBounds(model, {'MAR04306'}, -1000, 'l');
    model = changeRxnBounds(model, {'MAR04306'}, 0, 'u');
    model = changeRxnBounds(model, {'MAR04332','MAR04333','MAR04335','MAR04654','MAR04655'}, 0, 'l');
    model = changeRxnBounds(model, {'MAR07702'}, 0, 'l');
    model = changeRxnBounds(model, {'MAR02185'}, 0, 'u');
    model = changeRxnBounds(model, {'MAR02185'}, -1000, 'l');
    
    model = addReaction(model, 'ATP_maintenance', ...
        'metaboliteList', {'MAM01371c[c]', 'MAM01371m[m]', 'MAM02040c[c]', 'MAM02040m[m]', ...
                           'MAM01285c[c]', 'MAM01285m[m]', 'MAM02751c[c]', 'MAM02751m[m]', ...
                           'MAM02039c[c]', 'MAM02039m[m]'}, ...
        'stoichCoeffList', [-1; -1; -1; -1; 1; 1; 1; 1; 1; 1], ...
        'reversible', false);
    
    model = addReaction(model, 'NAD_demand', ...
        'metaboliteList', {'MAM02552c[c]', 'MAM02039c[c]', 'MAM02552m[m]', 'MAM02039m[m]', ...
                           'MAM02553c[c]', 'MAM02553m[m]'}, ...
        'stoichCoeffList', [-1; -1; -1; -1; 1; 1], ...
        'reversible', false);
    
    model = addReaction(model, 'NADPH_demand', ...
        'metaboliteList', {'MAM02555c[c]', 'MAM02555m[m]',...
                            'MAM02039c[c]', 'MAM02039m[m]', 'MAM02554c[c]', 'MAM02554m[m]'}, ...
        'stoichCoeffList', [-1; -1; 1; 1; 1; 1], ...
        'reversible', false);
    
    model.lb = model.lb ./ scaling_factor;
    model.ub = model.ub ./ scaling_factor;
    
    if strcmp(model_id, 'HSD')
        % Transport reactions
        transport_rxns = {'MAR09034','MAR09426','MAR05029'};
        for j = 1:length(transport_rxns)
            rxn_idx = find(strcmp(model.rxns, transport_rxns{j}));
            if ~isempty(rxn_idx)
                model = changeRxnBounds(model, transport_rxns{j}, ...
                    model.ub(rxn_idx) * reductionFactor, 'u');
                model = changeRxnBounds(model, transport_rxns{j}, ...
                    model.lb(rxn_idx) * reductionFactor, 'l');
            end
        end
        
        gene_list = {'Tret1-1','Glut1'};
        for g = 1:length(gene_list)
            try
                rxn_struct = findRxnsFromGenes(model, gene_list{g});
                fields = fieldnames(rxn_struct);
                for f = 1:length(fields)
                    rxns = rxn_struct.(fields{f})(:,1);
                    for r = 1:length(rxns)
                        rxn_idx = find(strcmp(model.rxns, rxns{r}));
                        if ~isempty(rxn_idx)
                            model = changeRxnBounds(model, rxns{r}, ...
                                model.ub(rxn_idx) * reductionFactor, 'u');
                            model = changeRxnBounds(model, rxns{r}, ...
                                model.lb(rxn_idx) * reductionFactor, 'l');
                        end
                    end
                end
            catch
                % Skip if gene not found
            end
        end
        
        ccm_rxns = {'MAR04373'}; 
        for j = 1:length(ccm_rxns)
            rxn_idx = find(strcmp(model.rxns, ccm_rxns{j}));
            if ~isempty(rxn_idx)
                model = changeRxnBounds(model, ccm_rxns{j}, ...
                    model.ub(rxn_idx) * 0.3, 'u');
                model = changeRxnBounds(model, ccm_rxns{j}, ...
                    model.lb(rxn_idx) * 0.3, 'l');
            end
        end
        
        tca_rxns = {'MAR04410', 'MAR04652','MAR08743', 'MAR04209','MAR05297','MAR06411','MAR06413'};
        reduction_factors = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]; 
        
        for j = 1:length(tca_rxns)
            rxn_idx = find(strcmp(model.rxns, tca_rxns{j}));
            if ~isempty(rxn_idx)
                model = changeRxnBounds(model, tca_rxns{j}, ...
                    model.ub(rxn_idx) * reduction_factors(j), 'u');
                model = changeRxnBounds(model, tca_rxns{j}, ...
                    model.lb(rxn_idx) * reduction_factors(j), 'l');
            end
        end
    end
    
end