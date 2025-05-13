%% Functional Comparison Across  Tissue-Specific GEMs
clc;
close all;

%% Define Save Folder and Load Tissue Info
save_dir = '5_analysis_fn_metabolic_task';
outputFolder = fullfile(pwd, save_dir);
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Load list of selected tissues
tissueInfoPath = fullfile(pwd, 'files', 'tissue_info.xlsx');
tissue_info_tbl = readtable(tissueInfoPath, 'ReadVariableNames', true);
tissue_info = table2cell(tissue_info_tbl(:, 1));  % Extract only the first column

%% Run Functional Comparison
save_results = 1;  % Set to 1 to save outputs
out = fn_run_compareGEM_fns_all(model_all, tissue_info, save_results, save_dir);

%% Compute and Plot Pass Rate
task_number = fullfile(pwd, 'files', 'metabolic_task_Full.xlsx');
task_number = readtable(task_number);

task_total = length(unique(task_number.ID(~isnan(task_number.ID))));

resultPath = strcat(pwd, '/5_analysis_fn_metabolic_task', '/availabel_fns_all.xlsx');
df = readtable(resultPath);


pass_rate = df.y2 ./ task_total;

df_plot = table(df.Tissue, pass_rate, 'VariableNames', {'Tissue', 'PassRate'});
writetable(df_plot, fullfile(save_dir, 'pass_rate.csv'));

% Plot
figure;
x = categorical(df_plot.Tissue);
x = reordercats(x, df_plot.Tissue);
bar(x, df_plot.PassRate);
grid on;
ylabel(['Pass Rate' newline '(Passed Tasks / Total Tasks)']);

% Save the figure
saveas(gcf, fullfile(save_dir, 'PassRate_BarPlot.png'));
saveas(gcf, fullfile(save_dir, 'PassRate_BarPlot.svg'));
