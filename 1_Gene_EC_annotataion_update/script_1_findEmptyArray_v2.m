%% Identify Reactions Without EC Numbers and Extract Associated KEGG IDs

clc;
clear;

%% Define Save Folder
pathway = pwd;
save_dir = '1_rxn_woEC';
subfolder = fullfile(pathway, save_dir);
if ~exist(subfolder, 'dir')
    mkdir(subfolder);
end

%% Load GEM Model
load(fullfile(pathway, 'models', 'fruitfly1.mat'));  % Assumes model saved as 'fruitflyGEM'
GEM_model = fruitflyGEM;

%% Run Analysis
[out_empty, uniqKEGGid, T] = run_findEmptyArray_v2(GEM_model,save_dir);

% Save outputs
movefile('out_empty.mat', subfolder);
writetable(T, fullfile(subfolder, 'raw_rxns_wo_ec.xlsx'));

%% Function: Identify Reactions Missing EC Numbers
function [out_empty, uniqKEGGid, T] = run_findEmptyArray_v2(GEM_model,save_dir)

out_empty = [];

for i = 1:length(GEM_model.eccodes)
    if isempty(GEM_model.eccodes{i})
        idx_empty = i;
        reaction = GEM_model.rxns{idx_empty};
        keggID   = GEM_model.rxnKEGGID{idx_empty};
        grRule   = GEM_model.grRules{idx_empty};
        out_empty = [out_empty; {idx_empty, reaction, keggID, grRule}];
    end
end

out_empty = string(out_empty);
fprintf('\n# of reactions without EC numbers: %d\n', size(out_empty,1));
fprintf('Total number of reactions        : %d\n', length(GEM_model.rxns));
fprintf('Percentage without EC            : %.1f%%\n', size(out_empty,1) / length(GEM_model.rxns) * 100);

% Create table for export
T = table(out_empty);
T = splitvars(T);
T.Properties.VariableNames = {'idx_empty', 'Reaction', 'KEGG_ID', 'grRules'};

% Extract KEGG IDs
keggIDs_all = out_empty(:,3);
keggIDs_nonempty = keggIDs_all(~cellfun('isempty', cellstr(keggIDs_all)));
uniqKEGGid = unique(keggIDs_nonempty, 'stable');
uniqKEGGid_table = table(uniqKEGGid, 'VariableNames', {'KEGG_ID'});
writetable(uniqKEGGid_table, strcat(save_dir,'\uniqKEGGid_list.csv'));

fprintf('\n# of reactions with KEGG ID      : %d\n', length(keggIDs_nonempty));
fprintf('# of reactions without KEGG ID   : %d\n', length(keggIDs_all) - length(keggIDs_nonempty));
fprintf('%% of EC-missing reactions with KEGG ID: %.1f%%\n', length(keggIDs_nonempty) / length(keggIDs_all) * 100);

% Save for future use
save('out_empty.mat', 'out_empty');

end
