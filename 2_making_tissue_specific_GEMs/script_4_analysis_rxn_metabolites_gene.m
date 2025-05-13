%% Script: Analyze Reactions, Metabolites, and Genes Across Models

% clc;
% clear;
% close all;

%% Create Output Folder
outputFolder = '4_analysis_rxn_metabolites_gene';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% Initialize Summary Variables
nModels = length(model_all);
id           = cell(nModels, 1);
reactions    = zeros(nModels, 1);
metabolites  = zeros(nModels, 1);
genes        = zeros(nModels, 1);

%% Extract Model Data
for i = 1:nModels
   reactions(i)   = length(model_all{i}.rxns);
   metabolites(i) = length(model_all{i}.mets);
   genes(i)       = length(model_all{i}.genes);
   id{i}          = model_all{i}.id; 
end

%% Save Model Summary Table
model_summary = table(id, reactions, metabolites, genes, ...
                      'VariableNames', {'ID', 'Reactions', 'Metabolites', 'Genes'});

writetable(model_summary, fullfile(outputFolder, 'model_summary.xlsx'));
disp('Model summary saved.');

%% Prepare Data for Visualization
bar_data      = [reactions, metabolites, genes];
bar_data_norm = bar_data ./ sum(bar_data, 2);  % Normalize per model

%% 1. Grouped Bar Chart
figure;
bar(bar_data);
xticks(1:nModels);
xticklabels(id);
legend({'Reactions', 'Metabolites', 'Genes'}, 'Location', 'Best');
xlabel('Model ID');
ylabel('Count');
title('Comparison of Reactions, Metabolites, and Genes');
grid on;
saveas(gcf, fullfile(outputFolder, 'Grouped_Bar.png'));
saveas(gcf, fullfile(outputFolder, 'Grouped_Bar.svg'));

%% 2. Stacked Bar Chart (Normalized)
figure;
bar(bar_data_norm, 'stacked');
xticks(1:nModels);
xticklabels(id);
legend({'Reactions', 'Metabolites', 'Genes'}, 'Location', 'Best');
xlabel('Model ID');
ylabel('Proportion');
title('Proportional Comparison of Reactions, Metabolites, and Genes');
grid on;
saveas(gcf, fullfile(outputFolder, 'Stacked_Bar.png'));
saveas(gcf, fullfile(outputFolder, 'Stacked_Bar.svg'));

%% 3. Correlation & Linear Fit Analysis

% Spearman Correlation
[rho_rm, pval_rm] = corr(reactions, metabolites, 'Type', 'Spearman');
[rho_rg, pval_rg] = corr(reactions, genes,      'Type', 'Spearman');

% Linear Regression & R²: Reactions vs Metabolites
p_rm       = polyfit(reactions, metabolites, 1);
y_pred_rm  = polyval(p_rm, reactions);
R2_rm      = 1 - sum((metabolites - y_pred_rm).^2) / sum((metabolites - mean(metabolites)).^2);

% Linear Regression & R²: Reactions vs Genes
p_rg       = polyfit(reactions, genes, 1);
y_pred_rg  = polyval(p_rg, reactions);
R2_rg      = 1 - sum((genes - y_pred_rg).^2) / sum((genes - mean(genes)).^2);

% Save Correlation Table
corr_results = table( ...
    {'Reactions vs Metabolites'; 'Reactions vs Genes'}, ...
    [rho_rm; rho_rg], [pval_rm; pval_rg], ...
    [R2_rm; R2_rg], ...
    'VariableNames', {'Comparison', 'Spearman_Correlation', 'P_Value', 'R2'});

writetable(corr_results, fullfile(outputFolder, 'Spearman_R2_Correlation.xlsx'));
disp('Correlation analysis saved.');

%% 4. Scatter Plot: Reactions vs Metabolites
figure;
scatter(reactions, metabolites, 100, 'filled');
hold on;
x_rm = linspace(min(reactions), max(reactions), 100);
y_rm = polyval(p_rm, x_rm);
plot(x_rm, y_rm, 'k--', 'LineWidth', 1.5);
xlabel('Number of Reactions');
ylabel('Number of Metabolites');
title(sprintf('Reactions vs Metabolites (Spearman ρ = %.2f, R² = %.2f)', rho_rm, R2_rm));
grid on;
saveas(gcf, fullfile(outputFolder, 'Scatter_Reactions_vs_Metabolites.png'));
saveas(gcf, fullfile(outputFolder, 'Scatter_Reactions_vs_Metabolites.svg'));

%% 5. Scatter Plot: Reactions vs Genes
figure;
scatter(reactions, genes, 100, 'filled');
hold on;
x_rg = linspace(min(reactions), max(reactions), 100);
y_rg = polyval(p_rg, x_rg);
plot(x_rg, y_rg, 'k--', 'LineWidth', 1.5);
xlabel('Number of Reactions');
ylabel('Number of Genes');
title(sprintf('Reactions vs Genes (Spearman ρ = %.2f, R² = %.2f)', rho_rg, R2_rg));
grid on;
saveas(gcf, fullfile(outputFolder, 'Scatter_Reactions_vs_Genes.png'));
saveas(gcf, fullfile(outputFolder, 'Scatter_Reactions_vs_Genes.svg'));

disp('All figures and data saved successfully.');
