% To update the reaction bounds obtained from FVA analysis onto model.

%% Define the save folder
output_folder = '2_FBA_FVA_rxn_bound_updated';  % Output directory

pathway = pwd;
subfolder_out = strcat(pathway,'\', output_folder);
if ~exist(subfolder_out, 'dir')
    mkdir(subfolder_out)
end

%% Load the FVA results
load(strcat(pathway,'\', '1_FBA_FVA','\out_all.mat'));

%% Apply FVA-bounded fluxes to models
out_all_fvaBounded = cell(length(out_all), 1);

for i = 1:length(out_all)
    model_fva = out_all{i,1}.subSysModel;
    model_fva = changeRxnBounds(model_fva, model_fva.rxns, out_all{i,1}.models(:,1), 'l');  % min bounds
    model_fva = changeRxnBounds(model_fva, model_fva.rxns, out_all{i,1}.models(:,2), 'u');  % max bounds
    out_all_fvaBounded{i,1} = model_fva;
end

%% Save updated model set
output_folder = '2_FBA_FVA_rxn_bound_updated';  % Output directory

save(strcat(pathway,'\', output_folder, '\out_all_fvaBounded.mat'), 'out_all_fvaBounded');
