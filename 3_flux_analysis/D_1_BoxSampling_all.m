current_path = pwd;
model_sampling_in_name = {'NSD','HSD'};
%% load
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model_sampling_in = [];
for i = 1:length(model_sampling_in_name)
   fileName = strcat(current_path,'\C_sampling\',string(model_sampling_in_name(i)));
   model_sampling_in{i} = load(fileName);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter setting
graph = 0; % 1 is to graph all
x_range = [-10 50];
y_range = [0 0.5];
y_range_box = [-1000 1000];
scale = 'log'; %'linear' or 'lin
legName = model_sampling_in_name;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Set the inputs information %%%%%%%%%
file_name_tmp = 'Reactions for subsystem glycolysis_gluconeogenesis';
pathway_int = 'glycolysis';

% file_name_tmp = 'Reactions for subsystem tricarboxylic_acid_cycle_and_glyoxylate_dicarboxylate_metabolism';
% pathway_int = 'tca';

% file_name_tmp = 'Reactions for metabolite MAM02552c_nadC';
% pathway_int = 'nadC';



%% Load tsv files - for gene-rxns
load_folder = [current_path,'\rxn_path_files\'];
tmp1 = tdfread(strcat(load_folder,file_name_tmp,'.tsv'));
% Extract the equations
opts = delimitedTextImportOptions('Delimiter','\t','Encoding','UTF-8');
tmp1_table = readtable(strcat(load_folder,file_name_tmp,'.tsv'),opts);
tmp1_table.Properties.VariableNames=tmp1_table{1,:};
tmp1_table(1, :) = [];  % Remove the first row after assigning column names
reactions = string(tmp1.Reaction_ID);
reactions_ID = string(tmp1.Genes);
trimmedStr = strtrim(reactions_ID);

%% Load the files - to extract the equation
eqn = string(tmp1.Equation);
%% List of reactions
list = reactions;
%% Make graphs
[dataOut,subfolder] = run_graphBoxSampling_eachRxns(model_sampling_in,graph,x_range,y_range,legName,y_range_box,list,pathway_int,scale,trimmedStr,tmp1_table);

%% Extract the mean, min, max
t = table(dataOut.mean,dataOut.min,dataOut.max);
t = splitvars(t);
t(:,[7,8,9,10,13,14,15]) = []; % remove the redundant columns
t2 = t;
t2(:,8) = t(:,9); % make the min and max for CTR, instead of min-ctr and min-exp.
t2(:,9) = t(:,8); % similarly, change the order.
t2.Properties.VariableNames = {'Rxns','Gene','Eqn','ctr_mean','exp_mean','relChg','ctr_min','ctr_max','exp_min','exp_max'};

%% Calculate relChange for min and max
min_ctr_tmp = t2.ctr_min;
emptyCells = cellfun(@isempty, min_ctr_tmp); % Find empty cells and replace them with 0
min_ctr_tmp(emptyCells) = {0};
min_ctr_tmp = cell2mat(min_ctr_tmp); % Convert the cell array to a double array

min_exp_tmp = t2.exp_min;
emptyCells = cellfun(@isempty, min_exp_tmp); % Find empty cells and replace them with 0
min_exp_tmp(emptyCells) = {0};
min_exp_tmp = cell2mat(min_exp_tmp); % Convert the cell array to a double array

min_relchg = min_exp_tmp./min_ctr_tmp;

%% Calculate relChange for  max
max_ctr_tmp = t2.ctr_max;
emptyCells = cellfun(@isempty, max_ctr_tmp); % Find empty cells and replace them with 0
max_ctr_tmp(emptyCells) = {0};
max_ctr_tmp = cell2mat(max_ctr_tmp); % Convert the cell array to a double array

max_exp_tmp = t2.exp_max;
emptyCells = cellfun(@isempty, max_exp_tmp); % Find empty cells and replace them with 0
max_exp_tmp(emptyCells) = {0};
max_exp_tmp = cell2mat(max_exp_tmp); % Convert the cell array to a double array

max_relchg = max_exp_tmp./max_ctr_tmp;
t2_tmp = table(min_relchg,max_relchg);
t3 = [t2, t2_tmp];

% Replace empty elements in the 'Gene' variable with a placeholder value
placeholderValue = 'NA';  % You can use any value you prefer
isEmptyCellArray = cellfun(@isempty, t3.Gene);
isNestedEmptyCellArray = cellfun(@(x) all(cellfun(@isempty, x)), t3.Gene);
t3.Gene(isEmptyCellArray | isNestedEmptyCellArray) = {placeholderValue};
writetable(t3,strcat(subfolder,'\',pathway_int,'_summary.xlsx'));

%% Save the structure folder
save(fullfile(subfolder,strcat('dataOut_',pathway_int,'.mat')),'dataOut');
% Save the structure file into excel
excelFilename = strcat(subfolder,'\',pathway_int,'.xlsx');
structFieldnames = fieldnames(dataOut); % <--- where myStruct is your struct of data

for k = 1:length(structFieldnames)
    fieldname = structFieldnames{k};
   % Write header row using writetable
   table_tmp = table(dataOut.(fieldname));
   table_tmp = splitvars(table_tmp);
   if k == 1 % for gene cells
       % Flatten the nested cell array into a string
       table_tmp.Properties.VariableNames = {'rxn','gene','eqn'};
   elseif k == 2 || k == 3 % mean or median 
       table_tmp.Properties.VariableNames = {'rxn','gene','eqn','ctrVal','expVal','RelChg','RelChgLog2'};
   else
       table_tmp.Properties.VariableNames = {'rxn','gene','eqn','ctrVal','expVal'};
   end
    C = table2cell(table_tmp);
    varName = table_tmp.Properties.VariableNames;
    C = [varName ; C];
    writecell(C, excelFilename, 'Sheet', fieldname);

end
file_info = ['1st col is WT, 2nd col is yki' newline '1st and 2nd set are dataset of Ctr and Yki' newline ...
                '3rd coln is relative change, 4th column is the log change - order of magnitude change'];
writematrix(file_info,excelFilename,'Sheet','info');
writematrix(trimmedStr,excelFilename,'Sheet','gene_info');
writetable(tmp1_table,excelFilename,'Sheet','subSysInfo');


function [dataOut,subfolder] = run_graphBoxSampling_eachRxns(model_sampling_in,graph,x_range,y_range,legName,y_range_box,list,list_gene,scale,trimmedStr,tmp1_table)
close all
path_c = pwd;
subfolder = [path_c '\D_flux_analysis\' scale '\' list_gene];
if ~exist(subfolder, 'dir')
    mkdir(subfolder)
end


%% 
% Extract the model id
modelID = {};
for i = 1:length(model_sampling_in)
%     modelID{i} = string(model_sampling_in{1,i}.SamplingResults.modelID);
    modelID{i} = string(legName(i));
end

dataOut = struct(); % Initialize the data structure

% Initialize data storage for Excel export
NSD_data = [];
HSD_data = [];
list_labels = [];
list_gene = [];

for k = 1:length(list)
    dataAll_sample =[];
    g = [];
    stats = [];
    c_vec =[];
    
    dataOut.genes{k,1} = list(k);
    dataOut.mean{k,1} = list(k);
    dataOut.med{k,1} = list(k);
    dataOut.mode{k,1} = list(k);
    dataOut.min{k,1} = list(k);
    dataOut.max{k,1} = list(k); 
    dataOut.std{k,1} = list(k); 
    dataOut.range{k,1} = list(k); 
    dataOut.iqr{k,1} = list(k); 
    dataOut.twenty5{k,1} = list(k); 
    dataOut.seventy5{k,1} = list(k); 

    dataOut.genes{k,2} = trimmedStr(k);
    dataOut.mean{k,2} = trimmedStr(k);
    dataOut.med{k,2} = trimmedStr(k);
    dataOut.mode{k,2} = trimmedStr(k);
    dataOut.min{k,2} = trimmedStr(k);
    dataOut.max{k,2} = trimmedStr(k);    
    dataOut.std{k,2} = trimmedStr(k); 
    dataOut.range{k,2} = trimmedStr(k); 
    dataOut.iqr{k,2} = trimmedStr(k); 
    dataOut.twenty5{k,2} = trimmedStr(k);
    dataOut.seventy5{k,2} = trimmedStr(k);

    dataOut.genes{k,3} = string(tmp1_table.Equation(k));
    dataOut.mean{k,3} = string(tmp1_table.Equation(k));
    dataOut.med{k,3} = string(tmp1_table.Equation(k));
    dataOut.mode{k,3} = string(tmp1_table.Equation(k));
    dataOut.min{k,3} = string(tmp1_table.Equation(k));
    dataOut.max{k,3} = string(tmp1_table.Equation(k));        
    dataOut.std{k,3} = string(tmp1_table.Equation(k));    
    dataOut.range{k,3} = string(tmp1_table.Equation(k));    
    dataOut.iqr{k,3} = string(tmp1_table.Equation(k));    
    dataOut.twenty5{k,3} = string(tmp1_table.Equation(k));
    dataOut.seventy5{k,3} = string(tmp1_table.Equation(k));

c_color = flip(jet(length(model_sampling_in)));
c_color = [0 0 0; 1 0 0];

data_ctr = [];
for i = 1:length(model_sampling_in)

    idxTmp = find(contains(model_sampling_in{1,i}.samples.rxns,list(k)));
    dataTmp = model_sampling_in{1,i}.samples.points(idxTmp,:)';
    dataTmp = reshape(dataTmp,[size(dataTmp,2)*length(dataTmp),1]);
    
    stats_tmp = [mean(dataTmp),median(dataTmp),mode(dataTmp)];
    stats_min = min(dataTmp); 
    stats_max = max(dataTmp);
    stats = [stats; stats_tmp];
        
    % 2024. 01. 02. Extract quartile values directly from the data
    lowerQuartile = quantile(dataTmp, 0.25);
    upperQuartile = quantile(dataTmp, 0.75);
 
    dataAll_sample = [dataAll_sample; dataTmp];
    g_tmp = repmat(string(modelID(i)),size(dataTmp,2).*length(dataTmp),1);
    g = [g; g_tmp];
    
    data(i).val = dataTmp;
    data(i).g = g_tmp;
    
    dataOut.genes{k,1} = list(k);
    dataOut.mean{k,i+3} = mean(dataTmp);
    dataOut.med{k,i+3} = median(dataTmp);
    dataOut.mode{k,i+3} = mode(dataTmp);
    dataOut.min{k,i+3} = min(dataTmp);
    dataOut.max{k,i+3} = max(dataTmp); 
    dataOut.std{k,i+3} = std(dataTmp);
    dataOut.range{k,i+3} = range(dataTmp);
    dataOut.iqr{k,i+3} = iqr(dataTmp);
    dataOut.twenty5{k,i+3} = lowerQuartile;
    dataOut.seventy5{k,i+3} = upperQuartile;
    c_vec_tmp = ones(length(dataTmp),1).*c_color(i,:);
    c_vec = [c_vec; c_vec_tmp];
    
    if i ==1 %store the mean of the control
        data_ctr = dataTmp;
        NSD_data = [NSD_data,dataTmp];
        if ~isempty(dataTmp)
         list_labels = [list_labels, list(k)];
         list_gene = [list_gene, trimmedStr(k)];
        end
    end    
    if i ==2 % calculate the scaledChange 
        mean_of_ctr = mean(data_ctr);
        mean_of_exp = mean(dataTmp);
       rel_change = (dataTmp-mean_of_ctr)./mean_of_ctr*100;
       rel_change = (dataTmp)./mean_of_ctr;

       rel_change_ctr = (data_ctr-mean_of_ctr)./mean_of_ctr*100;
       % Calculate the minimum and maximum of the relative change
        min_rel_change = min(rel_change);
        max_rel_change = max(rel_change);
        
        HSD_data = [HSD_data,dataTmp];  
    end
    
end


if isempty(g) ~= 1
g_cat = categorical(g);
g_cat = reordercats(g_cat,{'NSD';'HSD'});

%% Calculate the mean
mean_of_exp = mean(dataTmp);
mean_all = [mean_of_ctr, mean_of_exp];

%% Boxplot graph
if graph == 1
    figure(1)
    bp2 = boxchart(g_cat,dataAll_sample,...
        'MarkerStyle','none',... %'MarkerSize','3',
        'BoxFaceColor','k','WhiskerLineStyle','-');
    bp2.JitterOutliers ='off';
    bp2.BoxFaceAlpha = 0.2;
    bp2.LineWidth = 2;
    bp2.MarkerColor = 'k';
    hold on
    ylabel('Flux (mmol g-DCW^{-1} hr^{-1})')
    ax = gca;
    ax.YGrid = 'on';
    yline(0,'--k','LineWidth',1.5)
    title(tmp1_table.Equation(k),'FontSize',10)
    ylim(y_range_box)
    set(gca,'yscale',scale,'fontsize',10)
    set(gca,'TickLabelInterpreter','tex');
    grid minor
    hold off   
end

filename_jpeg = char(strcat(subfolder,'\',list(k),'.jpg'));

    if graph ==1
        
        print(gcf, filename_jpeg, '-djpeg', '-r300')
        saveas(gcf, char(strcat(subfolder,'\',list(k),'.fig')))
        
  
    end
end
end

% Define the Excel filename
excel_filename = fullfile(subfolder, 'flux_data_raw.xlsx');

% Find columns in NSD_data that are **not all zeros**
nonzero_cols_NSD = any(NSD_data ~= 0, 1);
NSD_data_filtered = NSD_data(:, nonzero_cols_NSD);
list_labels_NSD = list_labels(nonzero_cols_NSD); % Adjust list labels accordingly
list_gene_NSD = list_gene(nonzero_cols_NSD);

% Find columns in HSD_data that are **not all zeros**
nonzero_cols_HSD = any(HSD_data ~= 0, 1);
HSD_data_filtered = HSD_data(:, nonzero_cols_HSD);
list_labels_HSD = list_labels(nonzero_cols_HSD); % Adjust list labels accordingly
list_gene_HSD = list_gene(nonzero_cols_HSD);

% Replace empty entries in list_gene_NSD and list_gene_HSD with "NA"
list_gene_NSD(ismissing(list_gene_NSD) | list_gene_NSD == "") = "NA";
list_gene_HSD(ismissing(list_gene_HSD) | list_gene_HSD == "") = "NA";

% Convert filtered data to cell format and add list_gene as the first row
NSD_combined = [cellstr(list_gene_NSD); num2cell(NSD_data_filtered)];
HSD_combined = [cellstr(list_gene_HSD); num2cell(HSD_data_filtered)];

% Convert to table
NSD_table = cell2table(NSD_combined, 'VariableNames', list_labels_NSD);
HSD_table = cell2table(HSD_combined, 'VariableNames', list_labels_HSD);

% Write to Excel
writetable(NSD_table, excel_filename, 'Sheet', 'NSD', 'WriteVariableNames', true);
writetable(HSD_table, excel_filename, 'Sheet', 'HSD', 'WriteVariableNames', true);

disp(['Excel file saved: ' excel_filename]);

    % Relative change
    tmp = cell2mat(dataOut.mean(:,end-1:end));
    relChange = (tmp(:,end)-tmp(:,end-1))./tmp(:,end-1)*100;
    relChange = (tmp(:,end))./tmp(:,end-1);
    dataOut.mean=[dataOut.mean, num2cell(relChange)];
    
    log_change = log10(tmp(:,end)./tmp(:,end-1));
    dataOut.mean=[dataOut.mean, num2cell(log_change)];
    
    % Relative change for median
    tmp = cell2mat(dataOut.med(:,end-1:end));
    relChange = (tmp(:,end)-tmp(:,end-1))./tmp(:,end-1)*100;
    relChange = (tmp(:,end))./tmp(:,end-1);
    dataOut.med=[dataOut.med, num2cell(relChange)];
    
    log_change = log10(tmp(:,end)./tmp(:,end-1));
    dataOut.med=[dataOut.med, num2cell(log_change)];

%% Check normality
% figure;
% subplot(1,2,1);
% histogram(NSD_data_filtered(:), 'Normalization', 'probability');
% hold on;
% histogram(HSD_data_filtered(:), 'Normalization', 'probability');
% legend('NSD', 'HSD');
% title('Flux Distribution');
% xlabel('Flux');
% ylabel('Probability');
% 
% subplot(1,2,2);
% boxplot([NSD_data_filtered(:), HSD_data_filtered(:)], {'NSD', 'HSD'});
% title('Boxplot of Flux Distributions');
% saveas(gcf, char(strcat(subfolder,'\normality.fig')))
% saveas(gcf, char(strcat(subfolder,'\normality.png')))

    
% eqn_filtered = tmp1_table.Equation(nonzero_cols_NSD);
% Initialize eqn_filtered as a cell array
eqn_filtered = cell(size(list_labels_NSD));

% Loop through each reaction ID in list_labels_NSD
for i = 1:length(list_labels_NSD)
    % Find the index where Reaction ID in tmp1_table matches list_labels_NSD
    match_idx = find(strcmp(tmp1_table.("Reaction ID"), list_labels_NSD(i)));
    
    % If a match is found, extract the corresponding Equation
    if ~isempty(match_idx)
        eqn_filtered{i} = tmp1_table.Equation{match_idx}; % Extract equation
    else
        eqn_filtered{i} = 'No Match Found'; % If no match, assign placeholder
    end
end

% Convert to a column cell array if needed
eqn_filtered = eqn_filtered(:);

% Display the first few equations
disp(eqn_filtered(1:min(10, length(eqn_filtered))));

% Initialize arrays
num_rxns = size(NSD_data_filtered, 2);
p_values = nan(num_rxns, 1);
median_NSD = nan(num_rxns, 1);
median_HSD = nan(num_rxns, 1);
effect_size = nan(num_rxns, 1);
relChange = nan(num_rxns, 1);
relChange_mean = nan(num_rxns, 1);
relChange_sd = nan(num_rxns, 1);
relChange_sem = nan(num_rxns, 1);
relChange_CI_lower = nan(num_rxns, 1);
relChange_CI_upper = nan(num_rxns, 1);
mode_NSD = nan(num_rxns, 1);
mode_HSD = nan(num_rxns, 1);

% Confidence level
confidence_level = 0.95;
alpha = 1 - confidence_level;

% Loop through each reaction
for i = 1:num_rxns
    NSD_flux = NSD_data_filtered(:, i);  % 10,000 x 1
    HSD_flux = HSD_data_filtered(:, i);  % 10,000 x 1

    % Mann-Whitney U test (ranksum test)
    p_values(i) = ranksum(NSD_flux, HSD_flux);
    
    % Compute median fluxes
    median_NSD(i) = median(NSD_flux);
    median_HSD(i) = median(HSD_flux);

    rounded_flux_nsd = round(NSD_flux, 3);  % round to 3 decimals
    [unique_vals, ~, idx_u] = unique(rounded_flux_nsd, 'stable');  % preserve order
    counts = accumarray(idx_u, 1);  % count frequencies
    max_count = max(counts);
    modes = unique_vals(counts == max_count);  % all tied modes
    mode_NSD(i) = modes(1);  % pick the first one
                    
    rounded_flux_hsd = round(HSD_flux, 3);  % round to 3 decimals
    [unique_vals, ~, idx_u] = unique(rounded_flux_hsd, 'stable');  % preserve order
    counts = accumarray(idx_u, 1);  % count frequencies
    max_count = max(counts);
    modes = unique_vals(counts == max_count);  % all tied modes
    mode_HSD(i) = modes(1);  % pick the first one
                        
    % Compute effect size (median difference)
    effect_size(i) = median_HSD(i) - median_NSD(i);
    
    % Compute log2 fold-change (avoid division by zero)
    if median_NSD(i) ~= 0
        relChange(i) = median_HSD(i) / median_NSD(i);
    else
        relChange(i) = NaN; % Undefined if median NSD is zero
    end
    
    % Compute individual relChange values per sample
    if median_NSD(i) ~= 0
        relChange_values = HSD_flux / median_NSD(i);
        relChange_mean(i) = mean(relChange_values);
        relChange_sd(i) = std(relChange_values);
        relChange_sem(i) = relChange_sd(i) / sqrt(length(HSD_flux)); % SEM = SD / sqrt(n)
        
        % Compute 95% confidence interval
        n = length(HSD_flux);
        if n > 1
            t_value = tinv(1 - alpha/2, n-1); % t-score for 95% CI
            relChange_CI_lower(i) = relChange_mean(i) - t_value * relChange_sem(i);
            relChange_CI_upper(i) = relChange_mean(i) + t_value * relChange_sem(i);
        else
            relChange_CI_lower(i) = NaN;
            relChange_CI_upper(i) = NaN;
        end
    else
        relChange_mean(i) = NaN;
        relChange_sd(i) = NaN;
        relChange_sem(i) = NaN;
        relChange_CI_lower(i) = NaN;
        relChange_CI_upper(i) = NaN;
    end
end

% Store results in a table
results_table = table(list_labels_NSD', list_gene_NSD', eqn_filtered, ...
                      median_NSD, median_HSD, effect_size, relChange, ...
                      relChange_mean, relChange_sd, relChange_sem, ...
                      relChange_CI_lower, relChange_CI_upper, p_values, ...
                      mode_NSD,mode_HSD,...
                      'VariableNames', {'Reaction', 'Gene', 'Equation', ...
                                        'Median_NSD', 'Median_HSD', ...
                                        'MeanDiff', 'RelChange', ...
                                        'Mean_RelChange', 'SD_RelChange', 'SEM_RelChange', ...
                                        'CI_Lower_95%', 'CI_Upper_95%', 'p_value',...
                                        'Mode_NSD','Mode_HSD'});

% Define Excel filename
excel_filename_pvals = fullfile(subfolder, 'Flux_Stats.xlsx');

% Save results to an Excel file
writetable(results_table, excel_filename_pvals, 'Sheet', 'Flux_Analysis');

disp(['Excel file saved with p-values: ', excel_filename_pvals]);

close all
end
