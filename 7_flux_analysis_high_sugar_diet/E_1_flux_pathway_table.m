%% Set the save and load folder
pathway = pwd;
loadfolder = [pathway '\' 'B_fva_update_v2'];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model_sampling_in_name = {'NSD','HSD'};
model_sampling_in = [];
for i = 1:length(model_sampling_in_name)
   fileName = strcat(pathway,'\C_sampling\',string(model_sampling_in_name(i)));
   model_sampling_in{i} = load(fileName);
    
end

%% Once loaded, deactivate;
load(strcat(loadfolder,'\','out_all_fvaBounded.mat'));
load(strcat(pwd,'\files\','subsysAll_fruitflyGEM.mat'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter setting
subfolder = [pathway '\' 'E_PFI_table'];
if ~exist(subfolder, 'dir')
    mkdir(subfolder)
end

subSysThreshold = 200; % flux difference is greathan 2 times

modelID = graph_box_max_min_metabolicSubsys(out_all_fvaBounded,subfolder,subsysAll_fruitflyGEM,subSysThreshold,model_sampling_in);

function modelID = graph_box_max_min_metabolicSubsys(out_all_fvaBounded,subfolder,subsysAll_fruitflyGEM,subSysThreshold,model_sampling_in)
close all
% Extract the mean values from the sampling analysis and normalize to the
% range.
matrix_row_subSys_col_tissue = [];
for j = 1:length(subsysAll_fruitflyGEM)
%     j=3
%% Extract the reaction of interest lb and ub
    lb = zeros(length(out_all_fvaBounded),1);
    ub = zeros(length(out_all_fvaBounded),1);
    pfi = zeros(length(out_all_fvaBounded),1); % store the range of Vmax-Vmin in the pathway
    
    for i = 1:length(out_all_fvaBounded)
        
        %% Find the index of specific subsystem from the flux_bound updated model
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
           
           %Extract the reaction name
            rxnName = model.rxns(idxTmp);
            % Match the reaction in the sampling model and extract the mean
            % values of the sampling.
            mean_val = zeros(length(rxnName),1);
            iqr_val = zeros(length(rxnName),1);
            for p = 1:length(rxnName)
                model_tmp = model_sampling_in{1,i}.samples;
                idxTmpTmp =find(ismember(model_tmp.rxns,rxnName(p)));
                mean_val(p,1) = mean(model_tmp.points(idxTmpTmp,:));
                iqr_val(p,1) = iqr(model_tmp.points(idxTmpTmp,:)); %interquartile range of data set
            end
        
            lb_temp = model.lb(idxTmp);
            ub_temp = model.ub(idxTmp); 
 
            pfi(i,1) = length(idxTmp).*sum(abs(ub_temp-lb_temp).*abs(mean_val),'omitnan');% Find the median of the range    
       else
            pfi(i,1) = 0;% Find the median            
       end 
    end
matrix_row_subSys_col_tissue{j,1} = pfi';

end
% row represents the subsystem and the columns are tissues
matrix_row_subSys_col_tissue = cell2mat(matrix_row_subSys_col_tissue);

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
fluxData = fluxData';
fluxTable = table(subsysAll_fruitflyGEM,fluxData);
fluxTable = splitvars(fluxTable);
% fluxTable.Properties.VariableNames = string(subsystemNames);

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

end
