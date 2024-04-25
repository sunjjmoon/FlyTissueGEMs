%% Goal:
% - To perform the fva analysis
% INput:
%   - 32 different models and the correspoinding results.
%% Set the directory
dir = 'A_fva_bounds_results_v2';
pathway = pwd;
subfolder = [pathway '\' dir];
if ~exist(subfolder, 'dir')
    mkdir(subfolder)
end

load(strcat(pwd,'\models\','model_out_cbra.mat')) 

% Perform flux variability analysis
for i = 1:length(model_out_cbra)
   model = model_out_cbra{i,1};
   [out,subSysModel] = run_fva(model);
   out_all{i,1}.models = out;
   out_all{i,1}.subSysModel = subSysModel;
   out_all{i,1}.modelID = model.modelID;
end

save('out_all.mat'); movefile('out_all.mat',strcat(subfolder,'\'));
save('subSysModel.mat'); movefile('subSysModel.mat',strcat(subfolder,'\'));

function [FVAsolution,subSysModel] = run_fva(model)
%% Change the reaction bounds (v2)
% Input:
%   -   model of the interest.
% GLUD - It is expected to be generating NADPH/NADH (Lewis, ARS, 2017)
    list = {'MAR03802','MAR03804'};
    model = changeRxnBounds(model,list,0,'u'); 
% IDH - These reactions are reversible.
    list = {'MAR03958', 'MAR00710', 'MAR04111', 'MAR04112'};
    model = changeRxnBounds(model,list,-1000,'l'); 
% Trxr-2 - Thioredoxin reductase uses NADPH to reduce oxidized thioredoxins
    list = {'MAR02358'};
    model = changeRxnBounds(model,list,-1000,'l'); 
    model = changeRxnBounds(model,list,0,'u'); 
% Zw - It is considered irreversible reaction (the default directionality is opposite)
    list = {'MAR04306'};
    model = changeRxnBounds(model,list,-1000,'l'); 
    model = changeRxnBounds(model,list,0,'u'); 
% Dhfr - difydrofolate makes THF. Set the reversible reaction to 0.
    list = {'MAR04332', 'MAR04333', 'MAR04335', 'MAR04654', 'MAR04655'};
    model = changeRxnBounds(model,list,0,'l'); 
% CG1236 / GRHPR - This reaction uses NAD(P)H to make glycolate. Set unidirectional.
    list = {'MAR07702'};
    model = changeRxnBounds(model,list,0,'l'); 
% FASN2 - NADPH is used for lipid synthesis in FASN2 reaction. The default directionality is opposite.  
    list = {'MAR02185'};
    model = changeRxnBounds(model,list,0,'u'); 
    model = changeRxnBounds(model,list,-1000,'l'); 

%% Set the new objective function
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

%% Set the solver
changeCobraSolver('gurobi','LP');
subSysModel = model;

model_in = changeObjective(model,{'ATP_maintenance','NAD_demand','NADPH_demand'},1/3); %ATP synthase       
FBAsolution  = optimizeCbModel(model_in,'max'); %FBA analysis
FBAsolution_out.sol = FBAsolution;

%% bondary flux change
% Run the FBA
FVAsolution = FBAsolution; 
[minFlux,maxFlux] = fluxVariability(subSysModel,90,'max',subSysModel.rxns,0,1);
FVAsolution = [minFlux,maxFlux];

end


