%% Set the save and load folder
pathway = pwd;
loadfolder = [pathway '\' 'B_fva_update_v2'];

load(strcat(loadfolder,'\','out_all_fvaBounded.mat'));
load(strcat(pwd,'\files\','subsysAll_fruitflyGEM.mat'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter setting
subfolder = [pathway '\' 'D_1'];
if ~exist(subfolder, 'dir')
    mkdir(subfolder)
end

subSysThreshold = 200; % flux difference is greathan 2 times

modelID = graph_box_max_min_metabolicSubsys(out_all_fvaBounded,subfolder,subsysAll_fruitflyGEM,subSysThreshold);

function modelID = graph_box_max_min_metabolicSubsys(out_all_fvaBounded,subfolder,subsysAll_fruitflyGEM,subSysThreshold)
close all


matrix_row_subSys_col_tissue = [];
for j = 1:length(subsysAll_fruitflyGEM)
%% Extract the reaction of interest lb and ub
    lb = zeros(length(out_all_fvaBounded),1);
    ub = zeros(length(out_all_fvaBounded),1);
    pfv = zeros(length(out_all_fvaBounded),1); % store the range of Vmax-Vmin in the pathway
    
    for i = 1:length(out_all_fvaBounded) % tissues
        
        model = out_all_fvaBounded{i,1};
        modelID{i,1} = model.modelID;
        subsys_model = model.subSystems;
    
        for k = 1:length(subsys_model)
            if numel(subsys_model{k,1}) ~= 1  % if there is more than 1 cell in the subsystem, use the first one
               subsys_model{k,1} = subsys_model{k,1}{1,1};
            end
        end
        subsys_model = string(subsys_model);
        [Lia,~] = ismember(subsys_model,subsysAll_fruitflyGEM(j));
        idxTmp_orig = find(Lia==1);
        idxTmp =find(ismember(subsys_model,subsysAll_fruitflyGEM(j)));
        
       if isempty(idxTmp) ~= 1 % If the reaction exists
            lb_temp = model.lb(idxTmp);
            ub_temp = model.ub(idxTmp);        
            pfv(i,1) = length(idxTmp).*sum(ub_temp-lb_temp);% Find the median of the range    
       else
            pfv(i,1) = 0;% Find the median            
       end 
    end
matrix_row_subSys_col_tissue{j,1} = pfv';

end
% row represents the subsystem and the columns are tissues
matrix_row_subSys_col_tissue = cell2mat(matrix_row_subSys_col_tissue);

% Calculate z-score of summed flux variability
fvs_zscore = zscore(matrix_row_subSys_col_tissue,0,2);

%% Initialize
% Create a cell array to store subsystem names
subsystemNames = cell(1, length(subsysAll_fruitflyGEM));
% Create a cell array to store tissue name
subsystemData = cell(length(modelID), length(subsysAll_fruitflyGEM));
% Create a cell array to store score data
scoreData = cell(length(modelID), length(subsysAll_fruitflyGEM));
% Create a zero array to store the sum of the score data 
scoreData_sum = zeros(length(modelID),1);
% Create a cell array to store flux data
fluxData = zeros(length(modelID), length(subsysAll_fruitflyGEM));

%% Sort by the highest to lowest
for i = 1:length(matrix_row_subSys_col_tissue)
   subsys_tmp = subsysAll_fruitflyGEM(i);
   val = matrix_row_subSys_col_tissue(i,:)';
   [B,I_tmp] = sort(val,'descend');
   tissue_sorted_noDistinguish4sameVal = modelID(I_tmp);
% Get the unique values and their indices
[uniqueValues, ~, indexInOriginal] = unique(val);
% Sort the unique values in descending order
sortedUniqueValues = sort(uniqueValues, 'descend');

    % Assign scores based on the order of sorted unique values
    scores = zeros(size(val));
    for k = 1:length(sortedUniqueValues)
        indices = find(indexInOriginal == k);
        scores(indices) = length(sortedUniqueValues) + 1 - k;
    end
    scoreData_sum = scoreData_sum +scores;
    scores = scores(I_tmp); % this is going to be reported
%% Save the data
% Store subsystem name
subsystemNames{i} = subsys_tmp;
% Store sorted tissues in subsystemData
subsystemData(:, i) =tissue_sorted_noDistinguish4sameVal;
% Store scores in scoreData
scoreData(:, i) = num2cell(scores);
fluxData(:,i) = val;

end
[scoreData_sum_sorted, idx_scoreSum] = sort(scoreData_sum,'ascend');
tissue_name_sorted_f = modelID(idx_scoreSum);

% Replace subsystem names with systematic names
for i = 1:length(subsystemNames)
    newSubsystemName = sprintf('subsys%d', i);
    subsystemNames{i} = newSubsystemName;
end

%% Report to the excel file
% Create a table for subsystem data
subsystemTable = cell2table(subsystemData,'VariableNames',string(subsystemNames));
% Create a table for score data
scoreTable = cell2table(scoreData,'VariableNames',string(subsystemNames));
% Create a table for score data
fluxTable = table(fluxData);
fluxTable = splitvars(fluxTable);
fluxTable.Properties.VariableNames = string(subsystemNames);
% Create a table for fva data
fvs_zscoreDataTable = table(fvs_zscore');
fvs_zscoreDataTable = splitvars(fvs_zscoreDataTable);
fvs_zscoreDataTable.Properties.VariableNames = string(subsystemNames);
% Create a tissue table
tissue_table = table(modelID);

% Creat a table for the summed score data
scoreSumTable = table(tissue_name_sorted_f,scoreData_sum_sorted);
% Create a table for the original subsystem name
subsysOrig = table(subsysAll_fruitflyGEM);

% Write tables to an Excel file with two sheets
writetable(subsystemTable, strcat(subfolder,'\','output.xlsx'), 'Sheet', 'Subsystem', 'WriteRowNames', true);
writetable(scoreTable, strcat(subfolder,'\','output.xlsx'), 'Sheet', 'Score', 'WriteRowNames', true);
writetable(fluxTable, strcat(subfolder,'\','output.xlsx'), 'Sheet', 'flux', 'WriteRowNames', true);
writetable(scoreSumTable, strcat(subfolder,'\','output.xlsx'), 'Sheet', 'scoreSum', 'WriteRowNames', true);
writetable(subsysOrig,strcat(subfolder,'\','output.xlsx'),'Sheet','SubSysOrig','WriteRowNames',true);
writetable(fvs_zscoreDataTable,strcat(subfolder,'\','output.xlsx'),'Sheet','FVS','WriteRowNames',true);
writetable(tissue_table,strcat(subfolder,'\','output.xlsx'),'Sheet','tissues','WriteRowNames',true);

end
