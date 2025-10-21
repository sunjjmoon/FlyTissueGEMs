%% script_01_run_pFBA
% Implements  parsimonious FBA, following the mathematical formulation
% descried in detail from Lewis et al., MSB, 2010.

clc;
% clear all;

%% Initialize COBRA Toolbox
% initCobraToolbox(false);
% changeCobraSolver('gurobi', 'LP'); % commented out like in original code


%% 1. Set up
file_loc = pwd;
load(strcat(file_loc,'\model_out_cbra_u.mat'));

models = model_out_cbra_u;

pathway = pwd;
results_dir = fullfile(pathway, '1_pFBA');

if ~exist(results_dir, 'dir')
    mkdir(results_dir); 
end


%% 2. Initialize
pFBA_results = cell(length(models), 1);
model_constrained_out = cell(length(models),1);

%% 3. Process each model
for i = 1:length(models)
    model = models{i,1};
    model_id = model.modelID;
    
    fprintf('Processing model %s (%d/%d)\n', model_id, i, length(models));
    
    % Apply constraints
    model_constrained = apply_constraints(model, model_id);
    model_constrained_out{i,1} = model_constrained;

    % Run pFBA
    fprintf('  Running proper pFBA...\n');
    pFBA_sol = run_proper_pFBA(model_constrained);
    
    % Store results
    pFBA_results{i,1}.modelID = model_id;
    pFBA_results{i,1}.solution = pFBA_sol;
    pFBA_results{i,1}.flux = pFBA_sol.x;
    pFBA_results{i,1}.objective = pFBA_sol.f;
    pFBA_results{i,1}.total_flux = sum(abs(pFBA_sol.x));
    pFBA_results{i,1}.avg_flux = mean(abs(pFBA_sol.x));
    
    fprintf('Model %s: Objective = %.4f, Total flux = %.4f\n', ...
        model_id, pFBA_sol.f, sum(abs(pFBA_sol.x)));
    fprintf('Model %s: Objective = %.4f, Avg flux = %.4f\n', ...
        model_id, pFBA_sol.f, mean(abs(pFBA_sol.x)));
end

%% Save results
save(fullfile(results_dir, 'pFBA_results.mat'), 'pFBA_results');
save(fullfile(results_dir, 'model_constrained_out.mat'), 'model_constrained_out');

% fprintf('Proper pFBA analysis completed. Results saved to pFBA_v2_results.mat\n');

%% Helper functions
function model = apply_constraints(model, model_id)
    % Apply the same constraints as in your FVA analysis
    
    % Set parameters
    reductionFactor = 0.5;
    scaling_factor = 1;
    
    % 1. NAD(P)H-dependent reactions constraints
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
    
    % 2. Add demand/maintenance reactions
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
    
    % 3. Apply scaling - No scaling is applied
    model.lb = model.lb ./ scaling_factor;
    model.ub = model.ub ./ scaling_factor;
    
    % 4. HSD-specific constraints
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
        
        % Additional constraints on central carbon metabolism
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
        
        % TCA cycle constraints
        tca_rxns = {'MAR04410', 'MAR04652','MAR08743', 'MAR04209','MAR05297','MAR06411','MAR06413'};
        reduction_factors = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]; % Fum1, SDH, SDH, OGDH, OGDH, OGDH, OGDH
        
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

function pFBA_sol = run_proper_pFBA(model)
    % Implement  pFBA following the mathematical formulation
    % min sum(v_irrev) s.t. max v_objective = v_objective,lb
    
    % 1. Run FBA to obtain optimal objective value
    FBA_sol = optimizeCbModel(model, 'max');
    
    if FBA_sol.stat ~= 1
        fprintf('    FBA failed, returning empty solution\n');
        pFBA_sol = FBA_sol;
        return;
    end
    
    optimal_obj = FBA_sol.f;
    fprintf('    Optimal objective: %.6f\n', optimal_obj);
    
    % 2. Convert to irreversible reactions
    model_irrev = convertToIrreversible(model);
%     fprintf('    Original reactions: %d, Irreversible reactions: %d\n', ...
%         length(model.rxns), length(model_irrev.rxns));
    
    % 3. Fix objective at optimal value with small tolerance
    obj_tolerance = 0.001; % 0.1% tolerance
    obj_rxns = find(model_irrev.c);
    
    for i = 1:length(obj_rxns)
        model_irrev = changeRxnBounds(model_irrev, model_irrev.rxns{obj_rxns(i)}, ...
            optimal_obj * (1 - obj_tolerance), 'l');
    end
    
    % 4. Minimize total flux    
    % Change objective to minimize sum of all fluxes
    % Since all reactions are now irreversible (non-negative), 
    % minimizing sum of fluxes = minimizing sum of absolute fluxes
    model_irrev.c = ones(length(model_irrev.rxns), 1); % Minimize sum of all fluxes
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
    pFBA_solution.f = optimal_obj; % Report the biological objective value
    pFBA_solution.stat = 1;
    pFBA_solution.origStat = pFBA_sol_irrev.origStat;
    
    pFBA_sol = pFBA_solution;
end