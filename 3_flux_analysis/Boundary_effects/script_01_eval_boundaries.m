%% Script_01_Evaluate Boundary Effects
% Goal: To evaluate how varying the boundary reactions influence the rates
% Input:
%   1. muscle-GEM (both NSD and HSD)
% Output:
% 
% Description:
% Test three boundary conditions: 1000 (default), 10000, and 50000
% Performs the same FVA and sampling analysis

clc;
% clear all;

%% 1. Set parameters
bound_conditions = [1000, 10000, 50000]; 

n_samples = 5000;        % Number of flux samples

%% 2. Load models
path_tmp = pwd;
parentDir = fileparts(path_tmp);            
model_out_cbra = load(fullfile(parentDir, 'models','model_out_cbra_u.mat'));
model_out_cbra = model_out_cbra.model_out_cbra_u;

%% 3. Test each boundary condition
for bound_idx = 1:length(bound_conditions)
    current_bound = bound_conditions(bound_idx);    
    % Create output directories
    dir_fva = sprintf('01_FBA_FVA_bounds_%d', current_bound);
    dir_bounded = sprintf('02_FVA_bounded_bounds_%d', current_bound);
    dir_sampling = sprintf('03_FSA_bounds_%d', current_bound);
    
    pathway = [pwd '/results'];
    subfolder_fva = [pathway '/' dir_fva];
    subfolder_bounded = [pathway '/' dir_bounded];
    subfolder_sampling = [pathway '/' dir_sampling];
    
    % Create directories
    if ~exist(subfolder_fva, 'dir'), mkdir(subfolder_fva); end
    if ~exist(subfolder_bounded, 'dir'), mkdir(subfolder_bounded); end
    if ~exist(subfolder_sampling, 'dir'), mkdir(subfolder_sampling); end
    
    %% 1: fva
    out_all = cell(length(model_out_cbra), 1);
    for i = 1:length(model_out_cbra)
        model = model_out_cbra{i,1};
        model_id = model.modelID;

        % Run FVA with current boundary
        [out, subSysModel] = run_fva(model, model_id, current_bound);
        
        % Store results
        out_all{i,1}.models = out;
        out_all{i,1}.subSysModel = subSysModel;
        out_all{i,1}.modelID = model_id;
        
    end
    
    % Save FVA results
    save(fullfile(subfolder_fva, 'out_all.mat'), 'out_all');
    
    %% 2. Apply fva bounds to models
    out_all_fvaBounded = cell(length(out_all), 1);
    
    for i = 1:length(out_all)
        model_fva = out_all{i,1}.subSysModel;
        model_fva = changeRxnBounds(model_fva, model_fva.rxns, out_all{i,1}.models(:,1), 'l');  % min bounds
        model_fva = changeRxnBounds(model_fva, model_fva.rxns, out_all{i,1}.models(:,2), 'u');  % max bounds
        out_all_fvaBounded{i,1} = model_fva;
    end
    
    % Save FVA-bounded models
    save(fullfile(subfolder_bounded, 'out_all_fvaBounded.mat'), 'out_all_fvaBounded');
    
    %% 3. Flux sampling
    changeCobraSolver('gurobi','all');
    solverOK = changeCobraSolver('gurobi','LP');
    
    for k = 1:length(out_all_fvaBounded)
        subSysModel = out_all_fvaBounded{k,1};
        id = subSysModel.modelID;
        
        samples = gpSampler(subSysModel, n_samples);
        
        save_name = fullfile(subfolder_sampling, [char(id), '.mat']);
        save(save_name, 'samples');
        
    end
    
    fprintf('Completed analysis for boundary condition: ±%d\n', current_bound);
end

%% 4. Compare results across Boundary Conditions
compare_boundary_effects(bound_conditions, pathway);
fprintf('\nAll boundary testing completed!\n');

function [FVAsolution, subSysModel] = run_fva(model, model_id, flux_bound)
reductionFactor = 0.5;    % Reduction factor for hsd
scaling_factor = 1;       % Scaling factor for all flux bounds
model.lb(model.lb == -1000) = -flux_bound; % Change the original flux bounds to the current flux bound
model.ub(model.ub == 1000) = flux_bound;
%  GLUD reactions 
model = changeRxnBounds(model, {'MAR03802','MAR03804'}, 0, 'u');
%  IDH reactions 
model = changeRxnBounds(model, {'MAR03958','MAR00710','MAR04111','MAR04112'}, -flux_bound, 'l');
%  Trxr-2 
model = changeRxnBounds(model, {'MAR02358'}, -flux_bound, 'l');
model = changeRxnBounds(model, {'MAR02358'}, 0, 'u');
%  Zw (G6PD)
model = changeRxnBounds(model, {'MAR04306'}, -flux_bound, 'l');
model = changeRxnBounds(model, {'MAR04306'}, 0, 'u');
% Dhfr (Dihydrofolate reductase)
model = changeRxnBounds(model, {'MAR04332','MAR04333','MAR04335','MAR04654','MAR04655'}, 0, 'l');
% CG1236/GRHPR
model = changeRxnBounds(model, {'MAR07702'}, 0, 'l');
% FASN1 (fatty acid synthase)
model = changeRxnBounds(model, {'MAR02185'}, 0, 'u');
model = changeRxnBounds(model, {'MAR02185'}, -flux_bound, 'l');

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

for i = 1:length(model)
    list = {'MAR09034','MAR09426'};  
    for j = 1:length(list)
        list_in = list(j);
        idx_in = find(ismember(model.rxns,list_in));
        if strcmp(model_id,'HSD') == 1
            fprintf('%s, %s, lb: %.2f, ub: %.2f\n', model_id, char(list_in), model.lb(idx_in), model.ub(idx_in));
            model = changeRxnBounds(model,list_in,model.ub(idx_in).*reductionFactor,'u'); 
            model = changeRxnBounds(model,list_in,model.lb(idx_in).*reductionFactor,'l'); 
            fprintf('%s_update, %s, lb: %.2f, ub: %.2f\n', model_id, char(list_in), model.lb(idx_in), model.ub(idx_in));
        end
    end    
end

% Glucose transport (extracellular -> cytosol)
for i = 1:length(model)
    list = {'MAR05029'};
    for j = 1:length(list)
        list_in = list(j);
        idx_in = find(ismember(model.rxns,list_in));
        if strcmp(model_id,'HSD') == 1
            fprintf('%s, %s, lb: %.2f, ub: %.2f\n', model_id, char(list_in), model.lb(idx_in), model.ub(idx_in));
            model = changeRxnBounds(model,list_in,model.ub(idx_in).*reductionFactor,'u'); 
            model = changeRxnBounds(model,list_in,model.lb(idx_in).*reductionFactor,'l'); 
            fprintf('%s_update, %s, lb: %.2f, ub: %.2f\n', model_id, char(list_in), model.lb(idx_in), model.ub(idx_in));
        end
    end    
end

gene_int = {'Tret1-1','Glut1'};
for i = 1:length(model)
    list_strct = findRxnsFromGenes(model,gene_int);
    fieldName = fieldnames(list_strct);
    list = [];
    for k = 1:length(fieldName)
        list_tmp = list_strct.(fieldName{k})(:,1);
        list = [list; string(list_tmp)];
    end
    for j = 1:length(list)
        list_in = list(j);
        idx_in = find(ismember(model.rxns,list_in));
        if strcmp(model_id,'HSD') == 1
            model = changeRxnBounds(model,list_in,model.ub(idx_in).*reductionFactor,'u'); 
            model = changeRxnBounds(model,list_in,model.lb(idx_in).*reductionFactor,'l'); 
        end
    end    
end

% GAPDH constraints
for i = 1:length(model)
    list = {'MAR04373'};
    for j = 1:length(list)
        list_in = list(j);
        idx_in = find(ismember(model.rxns,list_in));
        if strcmp(model_id,'HSD') == 1
            model = changeRxnBounds(model,list_in,model.ub(idx_in).*0.3,'u'); 
            model = changeRxnBounds(model,list_in,model.lb(idx_in).*0.3,'l'); 
        end
    end    
end

% TCA-related genes
gene_int = {'Ogdh','SdhA','kdn'};
for i = 1:length(model)
    list_strct = findRxnsFromGenes(model,gene_int);
    fieldName = fieldnames(list_strct);
    list = [];
    for k = 1:length(fieldName)
        list_tmp = list_strct.(fieldName{k})(:,1);
        list = [list; string(list_tmp)];
    end
    for j = 1:length(list)
        list_in = list(j);
        idx_in = find(ismember(model.rxns,list_in));
        if strcmp(model_id,'HSD') == 1
            model = changeRxnBounds(model,list_in,model.ub(idx_in).*reductionFactor,'u'); 
            model = changeRxnBounds(model,list_in,model.lb(idx_in).*reductionFactor,'l'); 
        end
    end    
end

changeCobraSolver('gurobi', 'LP');
subSysModel = model;
FBAsolution = optimizeCbModel(subSysModel, 'max');
% Run FVA (90% optimal flux range)
[minFlux, maxFlux] = fluxVariability(subSysModel, 90, 'max', subSysModel.rxns, 0, 1);
FVAsolution = [minFlux, maxFlux];

end
