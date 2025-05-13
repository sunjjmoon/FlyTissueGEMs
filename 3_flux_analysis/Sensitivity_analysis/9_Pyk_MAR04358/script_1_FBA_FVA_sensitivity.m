%% Goal:
% - To perform the FVA analysis for HSD-muscle-GEM
% Input:
%   - Cell array of models: model_out_cbra_u.mat
%   - Apply flux constraints to simulate T2D (e.g. reduced glucose uptake, TCA cycle)
% Output:
%   - Flux variability results saved in /1_FBA_FVA

clc;
% clear all;

%% Set the output directory
dir = '1_FBA_FVA';
pathway = pwd;
subfolder = [pathway '/' dir];
if ~exist(subfolder, 'dir')
    mkdir(subfolder)
end

%% Set the rate reduction parameters
reductionFactor = 0.5;    % Reduction factor for specific reactions (e.g. glucose uptake) to make HSD-muscle-GEM
scaling_factor = 1;       % Scaling factor for all flux bounds; 1 = no global scaling.
sens_reductionF = 0.8;
sens_int = {'MAR04358'};

%% Load all models
[parentFolder, ~, ~] = fileparts(pathway);
model_out_cbra = load(fullfile(parentFolder,'models', 'model_out_cbra_u.mat'));
model_out_cbra = model_out_cbra.model_out_cbra_u;

%% Perform FVA for each model
for i = 1:length(model_out_cbra)
    model = model_out_cbra{i,1};
    model_id = model.modelID;

    % Run the core FVA pipeline
    [out, subSysModel] = run_fva(model, model_id, reductionFactor, scaling_factor,sens_reductionF,sens_int);

    % Store results
    out_all{i,1}.models = out;
    out_all{i,1}.subSysModel = subSysModel;
    out_all{i,1}.modelID = model_id;

    disp([num2str(i) ' / ' num2str(length(model_out_cbra)) ' is completed']);
end

%% Save results
save('out_all.mat');
movefile('out_all.mat', fullfile(subfolder, 'out_all.mat'));


function [FVAsolution, subSysModel] = run_fva(model, model_id, reductionFactor, scaling_factor,sens_reductionF,sens_int)

%% 1. Update constraints for NAD(P)H-dependent reactions
% Based on literature evidence
% detailed description is found in the supplementary data file.

%  GLUD reactions 
model = changeRxnBounds(model, {'MAR03802','MAR03804'}, 0, 'u');

%  IDH reactions 
model = changeRxnBounds(model, {'MAR03958','MAR00710','MAR04111','MAR04112'}, -1000, 'l');

%  Trxr-2 
model = changeRxnBounds(model, {'MAR02358'}, -1000, 'l');
model = changeRxnBounds(model, {'MAR02358'}, 0, 'u');

%  Zw (G6PD)
model = changeRxnBounds(model, {'MAR04306'}, -1000, 'l');
model = changeRxnBounds(model, {'MAR04306'}, 0, 'u');

% Dhfr (Dihydrofolate reductase)
model = changeRxnBounds(model, {'MAR04332','MAR04333','MAR04335','MAR04654','MAR04655'}, 0, 'l');

% CG1236/GRHPR
model = changeRxnBounds(model, {'MAR07702'}, 0, 'l');

% FASN1 (fatty acid synthase)
model = changeRxnBounds(model, {'MAR02185'}, 0, 'u');
model = changeRxnBounds(model, {'MAR02185'}, -1000, 'l');

%% 2. Objective functions: ATP, NAD, and NADPH demand reactions
% These reactions simulate energy/redox requirements of the cell

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

%% 3. Apply global flux scaling
model.lb = model.lb ./ scaling_factor;
model.ub = model.ub ./ scaling_factor;

%% 4. Apply constraints in transport reactions 
% MAR09034 / MAR09426 = Glucose / Trehalose [e] <-> 
% MAR05029 = glucose[c] <-> glucose[e]
for i = 1:length(model)
    list = {'MAR09034','MAR09426'};  
    for j = 1:length(list)
        list_in = list(j);
        idx_in = find(ismember(model.rxns,list_in));
        if strcmp(model_id,'HSD') == 1 % HSD
            % Check the Treh gene.
            disp([model_id ', ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
            model = changeRxnBounds(model,list_in,model.ub(idx_in).*reductionFactor,'u'); 
            model = changeRxnBounds(model,list_in,model.lb(idx_in).*reductionFactor,'l'); 
            disp([model_id '_update, ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
        else
            disp([newline model_id ', ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
            disp([model_id '_update, ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
        end
    end    
end

for i = 1:length(model)
    list = {'MAR05029'};  % glucose [c] <-> glucose[e]
    for j = 1:length(list)
        list_in = list(j);
        idx_in = find(ismember(model.rxns,list_in));
        if strcmp(model_id,'HSD') == 1 % HSD
            % Check the Treh gene.
            disp([model_id ', ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
            model = changeRxnBounds(model,list_in,model.ub(idx_in).*reductionFactor,'u'); 
            model = changeRxnBounds(model,list_in,model.lb(idx_in).*reductionFactor,'l'); 
            disp([model_id '_update, ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
        else
            disp([newline model_id ', ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
            disp([model_id '_update, ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
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
        if strcmp(model_id,'HSD') == 1 % HSD
            % Check the Treh gene.
            disp([model_id ', ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
            model = changeRxnBounds(model,list_in,model.ub(idx_in).*reductionFactor,'u'); 
            model = changeRxnBounds(model,list_in,model.lb(idx_in).*reductionFactor,'l'); 
            disp([model_id '_update, ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
        else
            disp([newline model_id ', ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
            disp([model_id '_update, ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
        end
    end    
end

%% 5. Apply constraints in central carbon metabolic reactions
% TCA-related genes (Ogdh, SdhA, kdn)
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
        if strcmp(model_id,'HSD') == 1 % HSD
            % Check the Treh gene.
            disp([model_id ', ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
            model = changeRxnBounds(model,list_in,model.ub(idx_in).*reductionFactor,'u'); 
            model = changeRxnBounds(model,list_in,model.lb(idx_in).*reductionFactor,'l'); 
            disp([model_id '_update, ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
        else
            disp([newline model_id ', ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
            disp([model_id '_update, ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
        end
    end    
end

%% Perturbation
% MAR04373 - 1,3-bisphospho-D-glycerate [c] + H+ [c] + NADH [c] Ô GAP [c] + NAD+ [c] + Pi [c]
for i = 1:length(model)
    list = sens_int;  % Sens
    for j = 1:length(list)
        list_in = list(j);
        idx_in = find(ismember(model.rxns,list_in));
        if strcmp(model_id,'HSD') == 1 % HSD
            % Check the Treh gene.
            disp([model_id ', ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
            model = changeRxnBounds(model,list_in,model.ub(idx_in).*sens_reductionF,'u'); 
            model = changeRxnBounds(model,list_in,model.lb(idx_in).*sens_reductionF,'l'); 
            disp([model_id '_update, ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
        else
            disp([newline model_id ', ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
            disp([model_id '_update, ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
        end
    end    
end

% Fum1
for i = 1:length(model)
    list = {'MAR04410'};  % Fum1
    for j = 1:length(list)
        list_in = list(j);
        idx_in = find(ismember(model.rxns,list_in));
        if strcmp(model_id,'HSD') == 1 % HSD
            % Check the Treh gene.
            disp([model_id ', ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
            model = changeRxnBounds(model,list_in,model.ub(idx_in).*0.1,'u'); 
            model = changeRxnBounds(model,list_in,model.lb(idx_in).*0.1,'l'); 
            disp([model_id '_update, ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
        else
            disp([newline model_id ', ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
            disp([model_id '_update, ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
        end
    end    
end


% SDH
for i = 1:length(model)
    list = {'MAR04652','MAR08743'};  % SDH
    for j = 1:length(list)
        list_in = list(j);
        idx_in = find(ismember(model.rxns,list_in));
        if strcmp(model_id,'HSD') == 1 % HSD
            % Check the Treh gene.
            disp([model_id ', ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
            model = changeRxnBounds(model,list_in,model.ub(idx_in).*0.5,'u'); 
            model = changeRxnBounds(model,list_in,model.lb(idx_in).*0.5,'l'); 
            disp([model_id '_update, ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
        else
            disp([newline model_id ', ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
            disp([model_id '_update, ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
        end
    end    
end

% OGDH
for i = 1:length(model)
    list = {'MAR04209','MAR05297','MAR06411','MAR06413'};  % SDH
    for j = 1:length(list)
        list_in = list(j);
        idx_in = find(ismember(model.rxns,list_in));
        if strcmp(model_id,'HSD') == 1 % HSD
            % Check the Treh gene.
            disp([model_id ', ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
            model = changeRxnBounds(model,list_in,model.ub(idx_in).*0.5,'u'); 
            model = changeRxnBounds(model,list_in,model.lb(idx_in).*0.5,'l'); 
            disp([model_id '_update, ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
        else
            disp([newline model_id ', ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
            disp([model_id '_update, ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
        end
    end    
end

%% 6. Run the simulation
% changeCobraSolver('gurobi', 'LP');
subSysModel = model;

% Set composite objective: 1/3 weight to each maintenance demand
model = changeObjective(model, {'ATP_maintenance','NAD_demand','NADPH_demand'}, 1/3);

% Run FBA
FBAsolution = optimizeCbModel(model, 'max');

% Run FVA (90% optimal flux range)
[minFlux, maxFlux] = fluxVariability(subSysModel, 90, 'max', subSysModel.rxns, 0, 1);
FVAsolution = [minFlux, maxFlux];

end
