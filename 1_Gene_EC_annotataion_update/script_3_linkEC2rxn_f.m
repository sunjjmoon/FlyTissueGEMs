%% KEGG EC Linking to GEM Reactions
% ------------------------------------------------------------
% Goal:
%   - Link KEGG-derived EC numbers to reactions missing EC annotations.
%   - Update GEM model accordingly.
% ------------------------------------------------------------

clc;

%% Define Save Folder
load(fullfile(pwd,'1_rxn_woEC','out_empty.mat'))
pathway = pwd;
save_dir = '3_linkEC2rxn';
subfolder = fullfile(pathway, save_dir);
if ~exist(subfolder, 'dir')
    mkdir(subfolder);
end

%% Load Required Files
fileIn = {fullfile(pathway, '2_parseKEGGID', 'out_kgID2ec_1to80'), ...
          fullfile(pathway, '2_parseKEGGID', 'out_kgID2ec_81toRest')};

load(fullfile(pathway, 'models', 'fruitfly2.mat'));  % loads model_uuu
gem = model_uuu;

%% Link EC to Reactions
[out_empty2, gem_u] = run_linkEC2rxn(out_empty, fileIn, gem,save_dir);

% Save outputs
movefile('out_empty2.mat', subfolder);
movefile('fruitfly2_ECupd.mat', subfolder);


%% Function: Link EC Numbers to GEM
function [out_empty2, gem_u] = run_linkEC2rxn(out_empty, fileIn, gem,save_dir)

% Load KEGG API outputs
df = [];
for i = 1:length(fileIn)
    file_path = strcat(fileIn{i}, '.txt');
    df_tmp = readtable(file_path, 'Delimiter', '\t', 'ReadVariableNames', false);
    df = [df; df_tmp];
end

df = table2cell(df);
df(:,1) = strrep(df(:,1), 'rn:', '');
df(:,2) = strrep(df(:,2), 'ec:', '');

fprintf('\n# of EC entries                  : %d\n', size(df,1));
fprintf('# of unique query KEGG IDs      : %d\n', length(unique(out_empty(:,3))));
fprintf('%% of matched EC entries         : %.1f%%\n', ...
    size(df,1) / length(unique(out_empty(:,3))) * 100);

% Concatenate ECs for KEGG IDs with multiple ECs
[unique_kgid, ia] = unique(df(:,1), 'stable');
df2 = df(ia,:);
for i = 1:length(unique_kgid)
    matching = ismember(df(:,1), unique_kgid{i});
    if sum(matching) > 1
        ec_combined = strjoin(string(df(matching,2)), ';');
        df2{i,2} = ec_combined;
    end
end

fprintf('# of KEGG IDs with ECs (including multiple): %d\n', size(df2,1));

% Match and assign ECs to out_empty
ec_col = strings(size(out_empty,1),1);
for i = 1:size(out_empty,1)
    idx = ismember(df2(:,1), out_empty(i,3));
    if any(idx)
        ec_col(i) = df2{idx,2};
    end
end

% Combine updated ECs into output table
out_empty2 = [out_empty, ec_col];
T = splitvars(table(out_empty2));
T.Properties.VariableNames = {'idx_empty', 'Reaction', 'KEGG_ID', 'grRules', 'EC'};
writetable(T, strcat(save_dir,'\raw_rxns_wo_ec_upd.xlsx'));

% Display summary
num_updated = sum(~cellfun('isempty', cellstr(out_empty2(:,5))));
fprintf('\n# of reactions updated with ECs : %d / %d (%.1f%%)\n', ...
    num_updated, size(out_empty2,1), num_updated / size(out_empty2,1) * 100);

% Update GEM model
non_empty_logic = ~cellfun('isempty', cellstr(out_empty2(:,5)));
update_idx = str2double(out_empty2(non_empty_logic,1));
ec_to_add = cellstr(out_empty2(non_empty_logic,5));

gem_u = gem;
for i = 1:length(update_idx)
    gem_u.eccodes(update_idx(i)) = ec_to_add(i);
end

gem_u.id = 'fruitfly2_ECupd';
save([gem_u.id, '.mat'], 'gem_u');
save('out_empty2.mat', 'out_empty2');

end
