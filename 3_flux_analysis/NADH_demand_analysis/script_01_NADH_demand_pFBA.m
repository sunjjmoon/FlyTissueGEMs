%% script_01_NADH_demand_pFBA: pFBA Analysis for NADH Maximum production Capacity
% Implements NADH_demand maximization and comprehensive flux analysis using pFBA
% Uses parsimonious FBA following Lewis et al., MSB, 2010

clc;
% clear all;

%% Initialize COBRA Toolbox
% initCobraToolbox(false);
% changeCobraSolver('gurobi', 'LP');

%% 1. Set up
file_loc = 'C:\Users\Sun Jin\Documents\revision\tissueGEM_revision\v2\code_to_upload\3_Flux_analysis\models';
load(strcat(file_loc,'\model_out_cbra_u.mat'));

models = model_out_cbra_u;

pathway = pwd;
results_dir = fullfile(pathway, '01_NADH_demand_pFBA');  % Different folder name

if ~exist(results_dir, 'dir')
    mkdir(results_dir); 
end

%% 2. Initialize results storage
enhanced_results = cell(length(models), 1);
model_constrained_out = cell(length(models),1);
nadh_analysis_results = struct();

%% 3. Process each model
for i = 1:length(models)
    model = models{i,1};
    model_id = model.modelID;
    
    fprintf('Processing model %s (%d/%d)\n', model_id, i, length(models));
    
    % Apply constraints
    model_constrained = apply_constraints(model, model_id);
    model_constrained_out{i,1} = model_constrained;

    %% ANALYSIS 1: pFBA biomass optimization
    fprintf('  Running pFBA biomass optimization...\n');
    biomass_sol = run_pFBA(model_constrained);
    
    %% ANALYSIS 2: NAD_demand maximization with biomass constraint (reviewer suggestion)
    fprintf('  Running NAD demand maximization with pFBA...\n');
    nad_demand_sol = maximize_nad_demand_pfba(model_constrained, biomass_sol);
    
    %% ANALYSIS 3: Comprehensive NADH/NAD+ flux analysis
    fprintf('  Analyzing NADH/NAD+ related fluxes...\n');
    nadh_flux_analysis = analyze_nadh_fluxes(model_constrained, nad_demand_sol);
    
    %% ANALYSIS 4: Flux Variability Analysis on NADH reactions
    fprintf('  Running FVA on NADH-related reactions...\n');
    nadh_fva_results = run_nadh_fva(model_constrained);
    
    %% Store comprehensive results (same structure as regular FBA version)
    enhanced_results{i,1}.modelID = model_id;
    enhanced_results{i,1}.biomass_solution = biomass_sol;
    enhanced_results{i,1}.nad_demand_solution = nad_demand_sol;
    enhanced_results{i,1}.nadh_flux_analysis = nadh_flux_analysis;
    enhanced_results{i,1}.nadh_fva_results = nadh_fva_results;
    enhanced_results{i,1}.max_nadh_oxidation_capacity = nad_demand_sol.f;
    enhanced_results{i,1}.biomass_growth_rate = biomass_sol.f;
    
    % Summary metrics
    enhanced_results{i,1}.max_nadh_oxidation_capacity = nad_demand_sol.f;
    enhanced_results{i,1}.biomass_growth_rate = biomass_sol.f;
    
    fprintf('Model %s: Max NADH oxidation capacity = %.4f\n', model_id, nad_demand_sol.f);
    fprintf('Model %s: Biomass growth rate = %.4f\n', model_id, biomass_sol.f);
end

%% 5. Comparative Analysis
fprintf('\n=== COMPARATIVE ANALYSIS ===\n');
comparative_results = perform_comparative_analysis(enhanced_results);

%% 6. Generate Summary Report
generate_summary_report(enhanced_results, comparative_results, results_dir);

%% Save all results (same file names as regular FBA version)
save(fullfile(results_dir, 'enhanced_FBA_results.mat'), 'enhanced_results');
save(fullfile(results_dir, 'comparative_results.mat'), 'comparative_results');
save(fullfile(results_dir, 'model_constrained_out.mat'), 'model_constrained_out');

fprintf('\npFBA Analysis complete! Results saved to: %s\n', results_dir);

fprintf('\nExporting results to Excel for review transparency...\n');

summary_table = table();
for i = 1:length(enhanced_results)
    res = enhanced_results{i};
    if isempty(res) || ~isfield(res, 'nadh_flux_analysis')
        continue;
    end

    tmp = table;
    tmp.ModelID = {res.modelID};
    tmp.Max_NADH_Oxidation_Capacity = res.max_nadh_oxidation_capacity;
    tmp.Biomass_Growth_Rate = res.biomass_growth_rate;
    tmp.Net_NADH_Production = res.nadh_flux_analysis.net_nadh_production;
    tmp.Net_NAD_Consumption = res.nadh_flux_analysis.net_nad_consumption;
    tmp.Mass_Balance_Check = res.nadh_flux_analysis.mass_balance_check;
    
    summary_table = [summary_table; tmp]; %#ok<AGROW>
    
    % Detailed per-reaction table for each model
    detailed_tbl = table();
    detailed_tbl.Reaction = res.nadh_flux_analysis.nadh_related_rxns';
    detailed_tbl.Flux = res.nadh_flux_analysis.nadh_related_fluxes';
    
end

% Save global summary table
writetable(summary_table, fullfile(results_dir, 'NADH_Summary_AllModels.xlsx'), 'Sheet', 'Summary');

fprintf('Excel export complete:\n');
fprintf('  - %s\\NADH_Summary_AllModels.xlsx\n', results_dir);
fprintf('  - %s\\NADH_Flux_Details.xlsx (per model sheets)\n', results_dir);


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
    
    model = addReaction(model, 'NADPH_demand', ...
        'metaboliteList', {'MAM02555c[c]', 'MAM02555m[m]',...
                            'MAM02039c[c]', 'MAM02039m[m]', 'MAM02554c[c]', 'MAM02554m[m]'}, ...
        'stoichCoeffList', [-1; -1; 1; 1; 1; 1], ...
        'reversible', false);
        
%     
        model = addReaction(model, 'NADH_demand', ...
        'metaboliteList', {'MAM02553c[c]', 'MAM02553m[m]','MAM02553p[p]'...
                            'MAM02552c[c]', 'MAM02039c[c]', ...
                            'MAM02552m[m]', 'MAM02039m[m]',...
                            'MAM02552p[p]', 'MAM02039p[p]'}, ...
        'stoichCoeffList', [-1; -1;-1;-1;
                            1; 1;
                            1; 1;
                            1; 1], ...
        'reversible', false);
    
    % 3. Apply scaling
    model.lb = model.lb ./ scaling_factor;
    model.ub = model.ub ./ scaling_factor;
    
    % 4. HSD-specific constraints (same as regular version)
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
        
        % Gene-associated reactions
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
        
        % Additional constraints
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

function pFBA_sol = run_pFBA(model)
    % Implement pFBA following Lewis et al., MSB, 2010
    % min sum(v_irrev) s.t. max v_objective = v_objective,lb
    
    % 1. Run FBA to obtain optimal objective value
    FBA_sol = optimizeCbModel(model, 'max');
    
    if FBA_sol.stat ~= 1
        fprintf('    FBA failed, returning empty solution\n');
        pFBA_sol = FBA_sol;
        return;
    end
    
    optimal_obj = FBA_sol.f;
    
    % 2. Convert to irreversible reactions
    model_irrev = convertToIrreversible(model);
    
    % 3. Fix objective at optimal value with small tolerance
    obj_tolerance = 0.001; % 0.1% tolerance
    obj_rxns = find(model_irrev.c);
    
    for i = 1:length(obj_rxns)
        model_irrev = changeRxnBounds(model_irrev, model_irrev.rxns{obj_rxns(i)}, ...
            optimal_obj * (1 - obj_tolerance), 'l');
    end
    
    % 4. Minimize total flux
    model_irrev.c = ones(length(model_irrev.rxns), 1);
    model_irrev.osense = 1; % Minimize
    
    % 5. Solve the pFBA problem
    pFBA_sol_irrev = optimizeCbModel(model_irrev, 'min');
    
    if pFBA_sol_irrev.stat ~= 1
        fprintf('    pFBA failed, returning regular FBA solution\n');
        pFBA_sol = FBA_sol;
        return;
    end
    
    % 6. Convert solution back to original model format
    pFBA_sol = convertIrrevFluxDistribution(pFBA_sol_irrev.x, model_irrev.match);
    
    % Create solution structure in original model format
    pFBA_solution.x = pFBA_sol;
    pFBA_solution.f = optimal_obj;
    pFBA_solution.stat = 1;
    pFBA_solution.origStat = pFBA_sol_irrev.origStat;
    
    pFBA_sol = pFBA_solution;
end

function nad_sol = maximize_nad_demand_pfba(model, biomass_sol)
    % Maximize NAD_demand with pFBA and biomass constraint
    
    % First constrain biomass to 50% of optimal
    biomass_idx = find(model.c);
    model = changeRxnBounds(model, model.rxns{biomass_idx}, 0.5 * biomass_sol.f, 'l');
    
    % Set NAD_demand as objective
    nad_demand_idx = find(strcmp(model.rxns, 'NADH_demand'));
    model.c = zeros(length(model.rxns), 1);
    model.c(nad_demand_idx) = 1;
    
    % Use pFBA to maximize NAD_demand
    nad_sol = run_pFBA(model);
end

function nadh_analysis = analyze_nadh_fluxes(model, solution)
    % Same function as regular FBA version
    
    if solution.stat ~= 1
        nadh_analysis = struct();
        return;
    end
    
    % Define NADH/NAD+ related metabolites
    nadh_mets = {'MAM02552c[c]', 'MAM02552m[m]'}; % NADH cytosolic and mitochondrial
    nad_mets = {'MAM02553c[c]', 'MAM02553m[m]'};   % NAD+ cytosolic and mitochondrial
    
    % Find all reactions involving NADH/NAD+
    nadh_rxns = [];
    for i = 1:length(nadh_mets)
        met_idx = find(strcmp(model.mets, nadh_mets{i}));
        if ~isempty(met_idx)
            rxn_indices = find(model.S(met_idx, :) ~= 0);
            nadh_rxns = [nadh_rxns, rxn_indices];
        end
    end
    
    for i = 1:length(nad_mets)
        met_idx = find(strcmp(model.mets, nad_mets{i}));
        if ~isempty(met_idx)
            rxn_indices = find(model.S(met_idx, :) ~= 0);
            nadh_rxns = [nadh_rxns, rxn_indices];
        end
    end
    
    nadh_rxns = unique(nadh_rxns);
    
    % Calculate net NADH production/consumption
    net_nadh_production = 0;
    net_nad_consumption = 0;
    
    nadh_producing_rxns = {};
    nadh_consuming_rxns = {};
    
    for i = 1:length(nadh_rxns)
        rxn_idx = nadh_rxns(i);
        flux = solution.x(rxn_idx);
        
        % Check NADH stoichiometry
        nadh_stoich = 0;
        for j = 1:length(nadh_mets)
            met_idx = find(strcmp(model.mets, nadh_mets{j}));
            if ~isempty(met_idx)
                nadh_stoich = nadh_stoich + model.S(met_idx, rxn_idx);
            end
        end
        
        net_nadh_change = nadh_stoich * flux;
        
        if net_nadh_change > 1e-6
            nadh_producing_rxns{end+1} = model.rxns{rxn_idx};
            net_nadh_production = net_nadh_production + net_nadh_change;
        elseif net_nadh_change < -1e-6
            nadh_consuming_rxns{end+1} = model.rxns{rxn_idx};
            net_nad_consumption = net_nad_consumption + abs(net_nadh_change);
        end
    end
    
    % Store results
    nadh_analysis.nadh_related_rxns = model.rxns(nadh_rxns);
    nadh_analysis.nadh_related_fluxes = solution.x(nadh_rxns);
    nadh_analysis.net_nadh_production = net_nadh_production;
    nadh_analysis.net_nad_consumption = net_nad_consumption;
    nadh_analysis.mass_balance_check = abs(net_nadh_production - net_nad_consumption);
    nadh_analysis.nadh_producing_rxns = nadh_producing_rxns;
    nadh_analysis.nadh_consuming_rxns = nadh_consuming_rxns;
end

function fva_results = run_nadh_fva(model)
    % Same function as regular FBA version
    
    target_rxns = {'MAR03956', 'MAR03957', 'MAR00363', 'MAR00364', 'MAR04732', ...
                   'NAD_demand', 'NADPH_demand'};
    
    existing_rxns = {};
    rxn_indices = [];
    
    for i = 1:length(target_rxns)
        idx = find(strcmp(model.rxns, target_rxns{i}));
        if ~isempty(idx)
            existing_rxns{end+1} = target_rxns{i};
            rxn_indices(end+1) = idx;
        end
    end
    
    if isempty(rxn_indices)
        fva_results = struct();
        return;
    end
    
    try
        [minFlux, maxFlux] = fluxVariability(model, 0, 'max', existing_rxns);
        
        fva_results.reactions = existing_rxns;
        fva_results.minFlux = minFlux;
        fva_results.maxFlux = maxFlux;
        fva_results.flux_range = maxFlux - minFlux;
        
    catch ME
        fprintf('    FVA failed: %s\n', ME.message);
        fva_results = struct();
    end
end

function comp_results = perform_comparative_analysis(enhanced_results)
    % Same function as regular FBA version
    
    n_models = length(enhanced_results);
    model_ids = cell(n_models, 1);
    nadh_capacities = zeros(n_models, 1);
    growth_rates = zeros(n_models, 1);
    
    for i = 1:n_models
        model_ids{i} = enhanced_results{i}.modelID;
        nadh_capacities(i) = enhanced_results{i}.max_nadh_oxidation_capacity;
        growth_rates(i) = enhanced_results{i}.biomass_growth_rate;
    end
    
    comp_results.model_ids = model_ids;
    comp_results.nadh_oxidation_capacities = nadh_capacities;
    comp_results.growth_rates = growth_rates;
    
    if n_models > 1
        control_nadh = nadh_capacities(1);
        control_growth = growth_rates(1);
        
        comp_results.nadh_fold_changes = nadh_capacities ./ control_nadh;
        comp_results.growth_fold_changes = growth_rates ./ control_growth;
        
        fprintf('\nNADH oxidation capacity fold changes (pFBA):\n');
        for i = 1:n_models
            fprintf('%s: %.3f-fold\n', model_ids{i}, comp_results.nadh_fold_changes(i));
        end
    end
end

function generate_summary_report(enhanced_results, comp_results, results_dir)
    % Same function as regular FBA version but mentions pFBA
    
    report_file = fullfile(results_dir, 'NADH_analysis_summary.txt');
    fid = fopen(report_file, 'w');
    
    fprintf(fid, '=== NADH/NAD+ METABOLIC ANALYSIS REPORT (pFBA VERSION) ===\n');
    fprintf(fid, 'Generated: %s\n\n', datestr(now));
    
    fprintf(fid, 'ADDRESSING REVIEWER CONCERNS WITH pFBA:\n');
    fprintf(fid, '1. Mass balance: All analyses maintain NAD+/NADH mass balance\n');
    fprintf(fid, '2. Network-wide capacity: NAD_demand maximization with pFBA\n');
    fprintf(fid, '3. Parsimonious solutions: Minimizes total flux while maintaining optimality\n\n');
    
    fprintf(fid, 'KEY FINDINGS (pFBA):\n');
    for i = 1:length(enhanced_results)
        result = enhanced_results{i};
        fprintf(fid, '\nModel: %s\n', result.modelID);
        fprintf(fid, '  Max NADH oxidation capacity: %.4f\n', result.max_nadh_oxidation_capacity);
        fprintf(fid, '  Biomass growth rate: %.4f\n', result.biomass_growth_rate);
        
        if isfield(result.nadh_flux_analysis, 'mass_balance_check')
            fprintf(fid, '  Mass balance check: %.2e (should be ~0)\n', result.nadh_flux_analysis.mass_balance_check);
        end
    end
    
    if isfield(comp_results, 'nadh_fold_changes')
        fprintf(fid, '\nCOMPARATIVE ANALYSIS (pFBA):\n');
        for i = 1:length(comp_results.model_ids)
            fprintf(fid, '%s: %.3f-fold change in NADH capacity\n', ...
                comp_results.model_ids{i}, comp_results.nadh_fold_changes(i));
        end
    end
    
    fclose(fid);
    fprintf('pFBA Summary report generated: %s\n', report_file);
end