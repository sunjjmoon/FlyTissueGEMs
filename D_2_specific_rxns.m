close all
pathway = pwd;
model_sampling_in_name = {'NSD','HSD'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model_sampling_in = [];
% for i = 1:length(model_sampling_in_name)
%    fileName = strcat(pathway,'\C_sampling\',string(model_sampling_in_name(i)));
%    model_sampling_in{i} = load(fileName);
%     
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter setting
graph = 1;
% x_range = [10^-1 10^2];
x_range = [-10 50];

y_range = [0 0.5];
y_range_box = [-1000 1000];
scale = 'lin'; %'linear', 'log'
scale = 'log'; %'linear', 'log'

legName = {'Ctr','Yki^{act}'};
legName = model_sampling_in_name;
legName = {'NSD', 'HSD'};
oppositeDir = 1; % -1 make the negative reaction to positivie


%% Specifiy the interest reaction
reactions = "MAR04388"; % LDH
% reactions = "MAR00479"; % GPDH - nadh
% reactions = "MAR04683"; % aldh

%%%%%% Set the inputs information %%%%%%%%%
file_name_tmp = 'Reactions for metabolite MAM02552c_nadC';
folder_int = 'nadC';
reactions = "MAR04388"; % LDH
reactions = "MAR00479"; % GPDH - nadh
reactions = "MAR04683"; % aldh
reactions = "MAR08507"; % T3dh


%% Load tsv files - for gene-rxns
load_folder = [pathway,'\rxn_path_files\'];
tmp1 = tdfread(strcat(load_folder,file_name_tmp,'.tsv'));

% Extract the equations
opts = delimitedTextImportOptions('Delimiter','\t','Encoding','UTF-8');
tmp1_table = readtable(strcat(load_folder,file_name_tmp,'.tsv'),opts);
tmp1_table.Properties.VariableNames=tmp1_table{1,:};
tmp1_table(1, :) = [];  % Remove the first row after assigning column names

%% 2023.12.11 - Find the specific reaction info
rxns_tmp = string(table2cell(tmp1_table(:,1)));
idx = find(ismember(rxns_tmp,(reactions)));
tmp1_table = tmp1_table(idx,:);

tmp1_table.Equation 
if oppositeDir == -1
% tmp1_table.Equation =  {'glucose [e] ⇔ glucose [c]'};
    tmp1_table.Equation =  strrep(tmp1_table.Equation,'⇒','⇔');
% Split the string into reactants and products
    splitReaction = strsplit(string(tmp1_table.Equation) , '⇔');
% Swap the reactants and products
    tmp1_table.Equation = [strtrim(splitReaction{2}), ' ⇔ ', strtrim(splitReaction{1})];
    tmp1_table.Equation = cellstr(tmp1_table.Equation);
end
% reactions = string(tmp1.Reaction_ID);
reactions_ID = string(tmp1.Genes);
reactions_ID = reactions_ID(idx);
trimmedStr = strtrim(reactions_ID);

%% Load the files (2023.12.11) - to extract the equation
eqn = string(tmp1.Equation);
%% List of reactions
list = reactions;
%% Make graphs
[dataOut,subfolder] = run_graphBoxSampling_eachRxns(model_sampling_in,graph,x_range,y_range,legName,...
    y_range_box,list,folder_int,scale,trimmedStr,tmp1_table,...
    oppositeDir);

%% Save the structure folder
save(fullfile(subfolder,strcat('dataOut_',folder_int,'.mat')),'dataOut');
% Save the structure file into excel
excelFilename = strcat(subfolder,'\',folder_int,'.xlsx');
structFieldnames = fieldnames(dataOut); % <--- where myStruct is your struct of data

for k = 1:length(structFieldnames)
    fieldname = structFieldnames{k};
    writecell(dataOut.(fieldname), excelFilename, 'Sheet', fieldname);
end
file_info = ['1st col is WT, 2nd col is yki' newline '1st and 2nd set are dataset of Ctr and Yki' newline ...
                '3rd coln is relative change, 4th column is the log change - order of magnitude change'];
writematrix(file_info,excelFilename,'Sheet','info');
writematrix(trimmedStr,excelFilename,'Sheet','gene_info');
writetable(tmp1_table,excelFilename,'Sheet','subSysInfo');

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

writetable(t3,strcat(subfolder,'\',folder_int,'_summary.xlsx'));

%%
function [dataOut,subfolder] = run_graphBoxSampling_eachRxns(model_sampling_in,graph,x_range,y_range,legName,...
    y_range_box,list,list_gene,scale,trimmedStr,tmp1_table,...
    oppositeDir)
close all
path_c = pwd;
subfolder = [path_c '\D_flux_analysis\specific\' list_gene];
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

    dataOut.genes{k,2} = trimmedStr(k);
    dataOut.mean{k,2} = trimmedStr(k);
    dataOut.med{k,2} = trimmedStr(k);
    dataOut.mode{k,2} = trimmedStr(k);
    dataOut.min{k,2} = trimmedStr(k);
    dataOut.max{k,2} = trimmedStr(k);    
    dataOut.std{k,2} = trimmedStr(k); 
    dataOut.range{k,2} = trimmedStr(k); 

    dataOut.genes{k,3} = string(tmp1_table.Equation(k));
    dataOut.mean{k,3} = string(tmp1_table.Equation(k));
    dataOut.med{k,3} = string(tmp1_table.Equation(k));
    dataOut.mode{k,3} = string(tmp1_table.Equation(k));
    dataOut.min{k,3} = string(tmp1_table.Equation(k));
    dataOut.max{k,3} = string(tmp1_table.Equation(k));        
    dataOut.std{k,3} = string(tmp1_table.Equation(k));    
    dataOut.range{k,3} = string(tmp1_table.Equation(k));    

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
    c_vec_tmp = ones(length(dataTmp),1).*c_color(i,:);
    c_vec = [c_vec; c_vec_tmp];
    
    
    if i ==1 %store the mean of the control
        data_ctr = dataTmp;
    end    
    if i ==2 % calculate the scaledChange 
        mean_of_ctr = mean(data_ctr);
        mean_of_exp = mean(dataTmp);
       rel_change = (dataTmp-mean_of_ctr)./mean_of_ctr*100;
       rel_change_ctr = (data_ctr-mean_of_ctr)./mean_of_ctr*100;
       % Calculate the minimum and maximum of the relative change
        min_rel_change = min(rel_change);
        max_rel_change = max(rel_change);
    end
    
end

if isempty(g) ~= 1
g = strrep(g, 'NSD','CTR');
g_cat = categorical(g);
g_cat = reordercats(g_cat,{'CTR';'HSD'});

%% Calculate the mean
mean_of_exp = mean(dataTmp);
mean_all = [mean_of_ctr.*oppositeDir, mean_of_exp.*oppositeDir];

%% Boxplot graph
if graph == 1
    figure(1)
    bp2 = boxchart(g_cat,dataAll_sample.*oppositeDir,...
        'MarkerStyle','none',... %'MarkerSize','3',
        'BoxFaceColor','k','WhiskerLineStyle','-');
    bp2.JitterOutliers ='off';
    bp2.BoxFaceAlpha = 0.3;
    bp2.LineWidth = 2;
    bp2.MarkerColor = 'k';
    hold on
    ylabel('Flux (mmol g-DCW^{-1} hr^{-1})')
    ax = gca;
    ax.YGrid = 'on';
    yline(0,'--k','LineWidth',1.5)
    title(tmp1_table.Equation(k),'FontSize',10)
    ylim(y_range_box)
    set(gca, 'yscale', scale, 'fontsize', 12, 'FontWeight', 'bold');
    set(gca,'TickLabelInterpreter','tex');
    grid minor
    hold off   
    
% Create a table with g_cat and dataAll_sample
export_table = table(string(g_cat), dataAll_sample .* oppositeDir, ...
    'VariableNames', {'Group', 'FluxValues'});

% Define the output filename
excelFilename = fullfile(subfolder, strcat(list(k), '_boxplot_data.xlsx'));

% Save the table to an Excel file
writetable(export_table, excelFilename, 'Sheet', 'Boxplot_Data');

disp(['Boxplot data saved to: ' excelFilename]);
end

filename_jpeg = char(strcat(subfolder,'\',list(k),'.jpg'));

    if graph ==1    
        print(gcf, filename_jpeg, '-djpeg', '-r300')
        saveas(gcf, char(strcat(subfolder,'\',list(k),'.fig')))
    end
end

%% Check effect size : If 1, they are 0.1 - 0.3 → Small effect
%0.3 - 0.5 → Medium effect
%> 0.5 → Large effect. Little overlap. most of them are all lower or higher
U = ranksum(dataAll_sample(g_cat == "CTR"), dataAll_sample(g_cat == "HSD"));
n1 = sum(g_cat == "CTR");
n2 = sum(g_cat == "HSD");
effect_size = 1 - (2 * U) / (n1 * n2);
fprintf('Effect Size: %.3f\n', effect_size);


%% Statistics
[p_mwu, h_mwu] = ranksum(dataAll_sample(g_cat == "CTR"), dataAll_sample(g_cat == "HSD"));
disp(['Mann-Whitney U Test p-value: ', num2str(p_mwu)]);

figure;
bp = boxchart(g_cat, dataAll_sample .* oppositeDir, ...
    'MarkerStyle', 'none', 'BoxFaceColor', 'k', 'WhiskerLineStyle', '-');
bp.BoxFaceAlpha = 0.3;
bp.LineWidth = 2;
hold on;
yline(0, '--k', 'LineWidth', 1.5);
ylabel('Flux (mmol g-DCW^{-1} hr^{-1})');
ax = gca;
ax.YGrid = 'on';
set(gca, 'yscale', scale, 'fontsize', 10);
grid minor;
hold off;

% Add P-Value Annotation
text(1.5, max(dataAll_sample) * 0.9, ...
    ['p = ', num2str(p_mwu, '%.2g')], ...
    'FontSize', 12, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');

% Save Visualization
output_folder = fullfile(subfolder, 'figures');
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end
saveas(gcf, char(strcat(output_folder,'\',list(k), '_boxplot_stat_comparison.svg')));
saveas(gcf, char(strcat(output_folder,'\',list(k), '_boxplot_stat_comparison.png')));
disp('✅ Boxplot with p-value saved!');

end

    % Relative change
    tmp = cell2mat(dataOut.mean(:,end-1:end));
    relChange = (tmp(:,end)-tmp(:,end-1))./tmp(:,end-1)*100;
    dataOut.mean=[dataOut.mean, num2cell(relChange)];
    
    log_change = log10(tmp(:,end)./tmp(:,end-1));
    dataOut.mean=[dataOut.mean, num2cell(log_change)];
    
    % Relative change for median
    tmp = cell2mat(dataOut.med(:,end-1:end));
    relChange = (tmp(:,end)-tmp(:,end-1))./tmp(:,end-1)*100;
    dataOut.med=[dataOut.med, num2cell(relChange)];
    
    log_change = log10(tmp(:,end)./tmp(:,end-1));
    dataOut.med=[dataOut.med, num2cell(log_change)];

%% Histogram
bWidth = 3e-2;
NumBins = 20;
%% Frequency Histogram

stats_store = [];
leg = [];
min_max = [];

colors = {'k', 'r'};  % Define colors for NSD and HSD
alphas = [0.5, 0.5];  % Define transparency for NSD and HSD
alphas = [0.6, 0.8];  % Define transparency for NSD and HSD

for i = 1:length(model_sampling_in)
    
    data_tmp = data(i).val.*oppositeDir;
    idx_neg = find(data_tmp<0);
    data_tmp(idx_neg) = [];
    
    data_log = log10(data_tmp);
    
    mean_ctr = mean(data_log);
    med_ctr = median(data_log);
    mod_ctr = mode(data_log);
    
    stats_store(i,:) = [mean_ctr, med_ctr, mod_ctr];
    leg{i} = data(i).g(1);
   

    figure(2)
    h =  histogram(data_log,'BinWidth',bWidth, 'Normalization', 'probability', 'FaceColor', colors{i}, 'FaceAlpha', alphas(i));
    min_max(i) = [max(h.Values)];
    hold on
    
    [counts,binLocations]=hist(data_log);
    [~,I] = max(counts);
    freq_max(i) = counts(I)./sum(counts);
    min_x(i) = binLocations(1);
    max_x(i) = binLocations(end);

end

xlabel('log_{10} Flux (mmol g-DCW^{-1} hr^{-1})','FontWeight','bold')
ylabel('Frequency','FontWeight','bold')
title(tmp1_table.Equation(k),'FontSize',8)
legend(legName,'fontsize',10,'Location','northwest');
grid on
grid minor
legend boxoff
set(gca, 'fontsize',9)

filename_jpeg = char(strcat(subfolder,'\',list(k),'_hist.jpg'));

    if graph ==1
        
        print(gcf, filename_jpeg, '-djpeg', '-r300')
        saveas(gcf, char(strcat(subfolder,'\',list(k),'_hist.fig')))
        saveas(gcf, char(strcat(subfolder,'\',list(k),'_hist.svg')))
     
        set(gcf,'Units','inches');
        screenposition = get(gcf,'Position');
        set(gcf,...
            'PaperPosition',[0 0 screenposition(3:4)],...
            'PaperSize',[screenposition(3:4)]);
        print('-dpdf', '-painters', char(strcat(subfolder,'\',list(k),'_hist.pdf')));

    end
     
%% **📌 Save P-Value, Effect Size, and Statistics to Excel**
excel_filename = strcat(subfolder,'\',list(k), '_pvalue_statistics.xlsx');

% Create summary statistics table
stat_table = table(["CTR"; "HSD"], stats_store(:,1), stats_store(:,2), stats_store(:,3), ...
    'VariableNames', {'Group', 'Mean_log10Flux', 'Median_log10Flux', 'Mode_log10Flux'});

% Create p-value and effect size table
pvalue_table = table("Mann-Whitney U Test", p_mwu, effect_size, ...
    'VariableNames', {'Statistical_Test', 'P_Value', 'Effect_Size'});

% Write to Excel
writetable(stat_table, excel_filename, 'Sheet', 'Statistics');
writetable(pvalue_table, excel_filename, 'Sheet', 'P_Value');

disp(['✅ P-value, effect size, and statistics saved to: ', excel_filename]);
   
end
