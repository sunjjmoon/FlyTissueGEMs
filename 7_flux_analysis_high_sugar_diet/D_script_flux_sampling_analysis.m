
pathway = pwd;
save_dirFlux = 'D';
subfolder = [pathway '\' save_dirFlux];
if ~exist(subfolder, 'dir')
    mkdir(subfolder)
end

model_sampling_in_name = {'NSD','HSD'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model_sampling_in = [];
for i = 1:length(model_sampling_in_name)
   fileName = strcat(pathway,'\C_sampling\',string(model_sampling_in_name(i)));
   model_sampling_in{i} = load(fileName);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter setting
graph = 0;
% x_range = [10^-1 10^2];
x_range = [-10 50];

y_range = [0 0.5];
y_range_box = [-1000 1000];
scale = 'lin'; %'linear'
scale = 'log'; %'linear'

legName = {'Ctr','Yki^{act}'};
legName = model_sampling_in_name;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Set the inputs information %%%%%%%%%
file_name_tmp = 'Reactions for subsystem glycolysis_gluconeogenesis';
list_gene = 'glycolysis';

%% Load tsv files - for gene-rxns
load_folder = [pathway,'\rxn_path_files\'];
tmp1 = tdfread(strcat(load_folder,file_name_tmp,'.tsv'));
tmp1 = tdfread(strcat(load_folder,file_name_tmp,'.tsv'));

% Extract the equations
opts = delimitedTextImportOptions('Delimiter','\t','Encoding','UTF-8');
tmp1_table = readtable(strcat(load_folder,file_name_tmp,'.tsv'),opts);
tmp1_table.Properties.VariableNames=tmp1_table{1,:};
tmp1_table(1, :) = [];  % Remove the first row after assigning column names

reactions = string(tmp1.Reaction_ID);
reactions_ID = string(tmp1.Genes);
trimmedStr = strtrim(reactions_ID);

%% extract the equation
eqn = string(tmp1.Equation);
%% List of reactions
list = reactions;
%% Make graphs
[dataOut,subfolder] = run_graphBoxSampling_eachRxns(model_sampling_in,graph,x_range,y_range,legName,y_range_box,list,list_gene,scale,trimmedStr,tmp1_table);

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

writetable(t3,strcat(subfolder,'\',list_gene,'_summary.xlsx'));

%% Save the structure folder
save(fullfile(subfolder,strcat('dataOut_',list_gene,'.mat')),'dataOut');
% Save the structure file into excel
excelFilename = strcat(subfolder,'\',list_gene,'.xlsx');
structFieldnames = fieldnames(dataOut); % <--- where myStruct is your struct of data


for k = 1:length(structFieldnames)
    fieldname = structFieldnames{k};
       % Write header row using writetable
   table_tmp = table(dataOut.(fieldname));
   table_tmp = splitvars(table_tmp);
%    table_tmp.gene = cellfun(@(c) strjoin(c(~cellfun('isempty', c)), ', '), table_tmp.gene, 'UniformOutput', false);
   
   if k == 1 % for gene cells
       % Flatten the nested cell array into a string
%        table_tmp.gene = cellfun(@(c) strjoin(c, ', '), table_tmp.gene, 'UniformOutput', false);
       table_tmp.Properties.VariableNames = {'rxn','gene','eqn'};
   elseif k == 2 || k == 3 % mean or median 
       table_tmp.Properties.VariableNames = {'rxn','gene','eqn','ctrVal','expVal','RelChg','RelChgLog2'};
   else
       table_tmp.Properties.VariableNames = {'rxn','gene','eqn','ctrVal','expVal'};
   end
%     writetable(table_tmp, excelFilename, 'Sheet', (fieldname));
    C = table2cell(table_tmp);
    varName = table_tmp.Properties.VariableNames;
    C = [varName ; C];
%     writecell(dataOut.(fieldname), excelFilename, 'Sheet', fieldname);
    writecell(C, excelFilename, 'Sheet', fieldname);

end
file_info = ['1st col is WT, 2nd col is yki' newline '1st and 2nd set are dataset of Ctr and Yki' newline ...
                '3rd coln is relative change, 4th column is the log change - order of magnitude change'];
writematrix(file_info,excelFilename,'Sheet','info');
writematrix(trimmedStr,excelFilename,'Sheet','gene_info');
writetable(tmp1_table,excelFilename,'Sheet','subSysInfo');

%%
function [dataOut,subfolder] = run_graphBoxSampling_eachRxns(model_sampling_in,graph,x_range,y_range,legName,y_range_box,list,list_gene,scale,trimmedStr,tmp1_table)
close all
path_c = pwd;
subfolder = [path_c '\D\swarm\' scale '\' list_gene];
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

for k = 1:length(list)
    dataAll_sample =[];
    g = [];
    % data = struct([]);
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
%         dataOut.scaledChange{k,1} = list(k); 

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
%         dataOut.scaledChange{k,2} = trimmedStr(k); 

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

%         dataOut.scaledChange{k,3} = string(tmp1_table.Equation(k));    


% 2023. 11. 21. -- define the color for jet. not the manual.
c_color = flip(jet(length(model_sampling_in)));
c_color = [0 0 0; 1 0 0];

% c_color = [0 0.4470 0.7410; 0.8500 0.3250 0.0980];
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
%     dataOut.data{k,i} = dataTmp;
%     dataOut.zScore{k,i+3} = (dataTmp-mean(dataTmp))./std(dataTmp);

    c_vec_tmp = ones(length(dataTmp),1).*c_color(i,:);
    c_vec = [c_vec; c_vec_tmp];
    
    if i ==1 %store the mean of the control
        data_ctr = dataTmp;
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
        % Define the desired range for scaling (e.g., -1 to 1)
%         scaled_min = -1; scaled_max = 1;

%     % Scale the relative change to the desired range
%     scaled_change = scaled_min + ((rel_change - min_rel_change) / (max_rel_change - min_rel_change)) * (scaled_max - scaled_min);
%     scaled_change_ctr = scaled_min + ((rel_change_ctr - min_rel_change) / (max_rel_change - min_rel_change)) * (scaled_max - scaled_min);
    % Calculate the maximum absolute change
%     max_absolute_change = max(abs(rel_change));
%     symmetric_percent_change = rel_change / max_absolute_change;
% 
% 
%     dataOut.scaledChange{k,i+3} = symmetric_percent_change;
%     dataOut.scaledChange{k,i+2} = rel_change_ctr./max_absolute_change;

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
%     plot(mean_all,'-o','LineWidth',1.5)
    ylabel('Flux (mmol g-DCW^{-1} hr^{-1})')
    ax = gca;
    ax.YGrid = 'on';
    yline(0,'--k','LineWidth',1.5)
%     title(trimmedStr(k))
    title(tmp1_table.Equation(k),'FontSize',10)
    ylim(y_range_box)
    set(gca,'yscale',scale,'fontsize',10)
    set(gca,'TickLabelInterpreter','tex');
    grid minor
    hold off      
end

filename_jpeg = char(strcat(subfolder,'\',list(k),'.jpg'));

    if graph ==1
%         print(gcf, filename_tiff,'-dtiffn','-r300')
%         saveas(gcf,filename_tiff_fig)
        
        print(gcf, filename_jpeg, '-djpeg', '-r300')
        saveas(gcf, char(strcat(subfolder,'\',list(k),'.fig')))
        
  
    end
end
end
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

close all
end
