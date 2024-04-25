%% Define the save folder
pathway = pwd;
save_dir = 'C';
subfolder = [pathway '\' save_dir];
if ~exist(subfolder, 'dir')
    mkdir(subfolder)
end

%% Load the KeggID - EC numbers
fileIn = {strcat(pathway,'\B\','out_kgID2ec_1to80'),...
            strcat(pathway,'\B\','out_kgID2ec_81toRest')};

%% Once load, deactivate
load(strcat(pathway,'\models\','fruitfly2.mat'));
gem = model_uuu;
%%      
[out_empty2,gem_u] = run_linkEC2rxn(out_empty,fileIn,gem);
movefile('out_empty2.mat',subfolder);
movefile('fruitfly2_ECupd.mat',subfolder);

function [out_empty2,gem_u] = run_linkEC2rxn(out_empty,fileIn,gem)
% Goal:
%   - To link the EC obtained from the KEGG API to rxns.
% Input:
%   -1. out_empty : Generated from the run_findEmptyArray_v2 function. 
%   2. fileIn: This will be files obtained from KEGG
%       e.g.,  fileIn = {'out_kgID2ec_1to80', 'out_kgID2ec_81toRest'};
%   3. gem: model that will be updated
%       e.g,) gem = fruitfly-gem

%% Load the txt files from KEGG API
df = []; % will be concetanated by row
for i = 1:length(fileIn)
   file_tmp = strcat(fileIn{i}, '.txt');
   df_tmp = readtable(string(file_tmp),'Delimiter', '\t', 'ReadVariableNames', false);
   df = [df; df_tmp];
end

df = table2cell(df);
df(:,1) = strrep(df(:,1),'rn:',''); % remove 'ec' from 'Rxxxec'.
df(:,2) = strrep(df(:,2),'ec:',''); % remove 'ec' from 'Rxxxec'.
disp(['# of EC: ' num2str(length(df)), newline,...
        '# of unique query rxns: ' num2str(length(unique(out_empty(:,3)))), newline,...
        '% of Ec from  unique query rxns: ' num2str(length(df)./length(unique(out_empty(:,3)))*100)]);
%% Make one array for kegg-ID that has multiple ECs.
[unique_kgid, ia] = unique(df(:,1),'stable');
df2 = df(ia,:);
for i = 1:length(unique_kgid)
    % Check which array matches to the kegg-id that has EC number.
   [Lia,~] = ismember(df(:,1),unique_kgid(i,1));
   if sum(Lia) ~= 1 % meaning there are multiple ECs for unique kgid
       ec_tmp = df(Lia,2);
       ec_tmp = strjoin(string(ec_tmp),';'); % join the EC with ';'
       df2(i,2) = cellstr(ec_tmp); % update the ec
   end % no need else. Because, unique match has already linked to EC.
end

disp(['# of KeggID with EC concetanated (some have multiple EC): ' num2str(length(df2))]);
    
%% link the EC to reactions
ec_col = strings(length(out_empty),1);
for i = 1:length(out_empty) % loop the original data
    % Check which array matches to the kegg-id that has EC number.
   [Lia,~] = ismember(df2(:,1),out_empty(i,3));
   % update the 
   if sum(Lia) ~= 0 % that means, there are KeggID that needs to be updated
       ec_tmp = df2(Lia,2); % extract the corresponding EC.
       ec_col(i) = ec_tmp; % update the ec_col
   end
end

out_empty2 = [out_empty, ec_col];
T = table(out_empty2); 
T = splitvars(T);
T.Properties.VariableNames = {'idx_empty', 'reaction', 'keggID', 'grRules', 'EC'};
filename = 'raw_rxns_wo_ec_upd.xlsx';
% write(T,filename);

%% check how many reactions will be updated
tmp = cellstr(out_empty2(:,5));
disp(['# of reactions that will have new EC from none: ',...
        num2str(sum(~cellfun('isempty',tmp))), newline,...
        'out of initial no-EC linked ', num2str(length(out_empty2)), ' rxns', newline,...
        '% of reactions updated for EC from blank: ' num2str(sum(~cellfun('isempty',tmp))./length(out_empty2)*100)])

%% update the EC to the model
non_empty_logic = ~cellfun('isempty',cellstr(out_empty2(:,5))); % 5th col is EC info col.

non_empty_idx = str2double(out_empty2(non_empty_logic,1));% 1st is idx of the original gem associated with reactions.
ec_to_add = cellstr(out_empty2(non_empty_logic,5));

gem_u = gem;
for i = 1:length(ec_to_add)
    gem_u.eccodes(non_empty_idx(i),1) = ec_to_add(i);
end
gem_u.id = 'fruitfly2_ECupd';
save([strcat(gem_u.id,'.mat')],'gem_u');
save('out_empty2.mat','out_empty2')
end