%% Goal:
% - To perform the fva analysis
% INput:
%   - 31 different models and the correspoinding results.
%% Set the directory
dir = 'A_fva_bounds_results_v2';

pathway = pwd;
subfolder = [pathway '/' dir];
if ~exist(subfolder, 'dir')
    mkdir(subfolder)
end

%% Define the decrease of the rate factor
reductionFactor = 0.5; % reduce everythin by 1/2. E.g.,) glucose uptake. Treh -> glc.
reductionF_Overall_flux = 1; % Reduce the overall flux bounds to this factor since the glucose uptake was estimated to be around ~1 mmol/g/hr
model_out_cbra = model_out_cbra_u;

% Perform flux variability analysis
for i = 1:length(model_out_cbra)
   model = model_out_cbra{i,1};
   model_id = model.modelID;
   [out,subSysModel] = run_fva(model,model_id,reductionFactor,reductionF_Overall_flux);
   out_all{i,1}.models = out;
   out_all{i,1}.subSysModel = subSysModel;
   out_all{i,1}.modelID = model.modelID;
   disp([num2str(i) '/ ' num2str(length(model_out_cbra)) ' is completed'])
end

save('out_all.mat'); movefile('out_all.mat',strcat(subfolder,'/'));
save('subSysModel.mat'); movefile('subSysModel.mat',strcat(subfolder,'/'));

function [FVAsolution,subSysModel] = run_fva(model,model_id,reductionFactor,reductionF_Overall_flux)
%% 2023. 11. 20.  -- Bring the constraints used for the previous flux analysis
% Input:
%   -   model of the interest.
% GLUD - It is expected to be generating NADPH/NADH (Lewis, ARS, 2017)
    list = {'MAR03802','MAR03804'};
    model{p} = changeRxnBounds(model{p},list,0,'u'); 
% IDH - These reactions are reversible.
    list = {'MAR03958', 'MAR00710', 'MAR04111', 'MAR04112'};
    model{p} = changeRxnBounds(model{p},list,-1000,'l'); 
% Trxr-2 - Thioredoxin reductase uses NADPH to reduce oxidized thioredoxins
    list = {'MAR02358'};
    model{p} = changeRxnBounds(model{p},list,-1000,'l'); % 
    model{p} = changeRxnBounds(model{p},list,0,'u');
% Zw - It is considered irreversible reaction (the default directionality is opposite)
    list = {'MAR04306'};
    model{p} = changeRxnBounds(model{p},list,-1000,'l'); % 
    model{p} = changeRxnBounds(model{p},list,0,'u');
% Dhfr - difydrofolate makes THF. Set the reversible reaction to 0.
    list = {'MAR04332', 'MAR04333', 'MAR04335', 'MAR04654', 'MAR04655'};
    model{p} = changeRxnBounds(model{p},list,0,'l'); % 
% CG1236 / GRHPR - This reaction uses NAD(P)H to make glycolate. Set unidirectional.
    list = {'MAR07702'};
    model{p} = changeRxnBounds(model{p},list,0,'l'); % 
% FASN2 - NADPH is used for lipid synthesis in FASN2 reaction. The default directionality is opposite.  
    list = {'MAR02185'};
    model{p} = changeRxnBounds(model{p},list,0,'u'); % 
    model{p} = changeRxnBounds(model{p},list,-1000,'l'); %  

%%  2023. 11. 20. -- Set the new objective function

            model = addReaction(model,'ATP_maintenance',...
            'metaboliteList', ...
            {'MAM01371c[c]', 'MAM01371m[m]', 'MAM02040c[c]', 'MAM02040m[m]',... %ATPc,ATPm,H2Oc,H2Om
                'MAM01285c[c]', 'MAM01285m[m]','MAM02751c[c]','MAM02751m[m]',... %ADPc,ADPm,Pic,Pim
                'MAM02039c[c]', 'MAM02039m[m]'},...%Hc,Hm,NADPc,NADPm
            'stoichCoeffList', [-1; -1; -1; -1;...
                                        1; 1; 1; 1;...
                                            1; 1],...
            'reversible',false);

        model = addReaction(model,'NAD_demand',...
            'metaboliteList', ...
            {'MAM02552c[c]', 'MAM02039c[c]', 'MAM02552m[m]', 'MAM02039m[m]',... %NADc,Hc,NADm,Hm,
                'MAM02553c[c]', 'MAM02553m[m]'},...%NADHc,NADHm
            'stoichCoeffList', [-1; -1; -1; -1;...
                                            1; 1],...
            'reversible',false);

        model  = addReaction(model,'NADPH_demand',...
            'metaboliteList', ...
            {'MAM02555c[c]', 'MAM02555m[m]',... % NADPHc, NADPHm
                'MAM02039c[c]', 'MAM02039m[m]', 'MAM02554c[c]', 'MAM02554m[m]'},...%Hc,Hm,NADPc,NADPm
            'stoichCoeffList', [-1; -1;...
                                        1; 1; 1; 1;],...
            'reversible',false);

%% Reduce the reaction bounds
model.lb = model.lb./reductionF_Overall_flux;
model.ub = model.ub./reductionF_Overall_flux;
   
%% Apply the constraints for the T2D - reduced glucose/trehalose uptake rates
for i = 1:length(model)
    list = {'MAR09034','MAR09426'};  % glucose [e], Trehalose[e] <-> 
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

% MAR05029 - GLucose [c] <-> glucose[e]. Change this bound
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

%% through gene
% tret1-2 : Transporter of trahalose and 
% glut1
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

%% GAPDH constraint
% MAR04373 - 1,3-bisphospho-D-glycerate [c] + H+ [c] + NADH [c] Ô GAP [c] + NAD+ [c] + Pi [c]
for i = 1:length(model)
    list = {'MAR04373'};  % glucose [c] <-> glucose[e]
    for j = 1:length(list)
        list_in = list(j);
        idx_in = find(ismember(model.rxns,list_in));
        if strcmp(model_id,'HSD') == 1 % HSD
            % Check the Treh gene.
            disp([model_id ', ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
            model = changeRxnBounds(model,list_in,model.ub(idx_in).*0.3,'u'); 
            model = changeRxnBounds(model,list_in,model.lb(idx_in).*0.3,'l'); 
            disp([model_id '_update, ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
        else
            disp([newline model_id ', ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
            disp([model_id '_update, ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
        end
    end    
end

%% TCA cycle constraints
% tret1-2 : Transporter of trahalose and 
% glut1
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
            model = changeRxnBounds(model,list_in,model.ub(idx_in).*0.5,'u'); 
            model = changeRxnBounds(model,list_in,model.lb(idx_in).*0.5,'l'); 
            disp([model_id '_update, ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
        else
            disp([newline model_id ', ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
            disp([model_id '_update, ' char(list_in) ', lb: ' num2str(model.lb(idx_in)) ', ub:' num2str(model.ub(idx_in)) ])
        end
    end    
end

%% Set the solver
changeCobraSolver('gurobi','LP');
subSysModel = model;
model = changeObjective(model,{'ATP_maintenance','NAD_demand','NADPH_demand'},1/3); %ATP synthase       
%     model_in = changeObjective(model{i},'MAR00021'); %set the objective function
FBAsolution  = optimizeCbModel(model,'max'); %FBA analysis
FBAsolution_out.sol = FBAsolution;
%% bondary flux change
% Run the FBA
FVAsolution = FBAsolution; 
[minFlux,maxFlux] = fluxVariability(subSysModel,90,'max',subSysModel.rxns,0,1);
FVAsolution = [minFlux,maxFlux];

end


