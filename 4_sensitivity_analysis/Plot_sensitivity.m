clc
clear all

%% Set the save and load folder
pathway = pwd;

% Define save folder
subfolder = [pathway '\' 'plot_sensitivity'];
% subfolder = [pathway '\' 'plot_sensitivity_v2'];

if ~exist(subfolder, 'dir')
    mkdir(subfolder)
end

% Get all directory contents
dirContents = dir(pathway);
folderNames = {};

% Extract valid subfolder names
for i = 1:numel(dirContents)
    if dirContents(i).isdir && ~strcmp(dirContents(i).name, '.') && ~strcmp(dirContents(i).name, '..')
        folderNames{end+1} = dirContents(i).name;
    end
end

% Filter folder names that start with digits + underscore
filteredFolderNames = folderNames(startsWith(folderNames, ...
    {'1_','2_','3_','4_','5_','6_','7_','8_', '9_','10_','11_','12','13','14'}));

%% Extract gene names from folder names
extractedWords = cell(size(filteredFolderNames));
for i = 1:numel(filteredFolderNames)
    tokens = regexp(filteredFolderNames{i}, '_([^_]+)_', 'tokens');
    if ~isempty(tokens)
        extractedWords{i} = tokens{1}{1};
    else
        extractedWords{i} = '';
    end
end

disp('Extracted words:');
disp(extractedWords');

%% Load sensitivity data
sens = zeros(length(filteredFolderNames),1);
gene_names = [];

for i = 1:length(filteredFolderNames)
    name = filteredFolderNames{i};
    loadfolder = fullfile(pathway, name);
    df = readtable(fullfile(loadfolder, '3_sensitivity', 'summary_data_v3_avg.xlsx'), 'Sheet', 'Summary');
    sens(i,1) = df.sensitivity;
    gene = extractBetween(name, '_', '_');
    gene_names{i,1} = gene;
end

%% Sort sensitivity by descending order
[sens_ordered, I] = sort(sens, 'descend');
gene_names_ordered = string(gene_names(I));
gene_names_ordered_cat = categorical(gene_names_ordered);
gene_names_ordered_cat_reorder = reordercats(gene_names_ordered_cat, gene_names_ordered);

%% Plot
xlabText = '\Sensitivity to glycolysis flux';

figure(1)
h = barh(gene_names_ordered_cat_reorder, sens_ordered);
title(xlabText);
set(gca, 'FontSize', 12);
grid on
grid minor



%% Export results
t = table(gene_names_ordered_cat_reorder, sens_ordered);

writetable(t, fullfile(subfolder,'Sensitivity_out_v3_avg.xlsx'));
saveas(gcf, fullfile(subfolder,'perturbation_effect_plot_v3_avg.png'));
saveas(gcf, fullfile(subfolder,'perturbation_effect_plot_v3_avg.svg'));

% writetable(t, fullfile(subfolder,'Sensitivity_out_v3_avg_v2.xlsx'));
% saveas(gcf, fullfile(subfolder,'perturbation_effect_plot_v3_avg_v2.png'));
