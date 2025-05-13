function out = fn_run_compareGEM_fns_all(model_all, tissue_info, save_in, save_dir)
% RUN_COMPAREGEM_FNS_ALL
% Compare the functionality (metabolic task performance) across selected GEMs.
%
% INPUTS:
%   model_all    - Cell array of tissue-specific GEMs
%   tissue_info  - Cell array of tissue names (in order)
%   save_in      - Boolean (1 = save results)
%   save_dir     - Output directory name (relative path)
%
% OUTPUT:
%   out          - Plot handle for spy plot (optional)

%% Set Up Save Directory
subfolder = fullfile(pwd, save_dir);
if ~exist(subfolder, 'dir')
    mkdir(subfolder);
end

%% Load Metabolic Task Definitions
taskFilePath = fullfile(pwd, 'files', 'metabolic_task_Full.xlsx');
res_func = compareMultipleModels(model_all, false, false, [], true, taskFilePath);

%% Assign Tissue Names
useModels = tissue_info;
res_func.modelIDs = useModels;

%% Identify Differentiated Tasks
isDiff = ~all(res_func.funcComp.matrix == 0, 2) & ~all(res_func.funcComp.matrix == 1, 2);
diffTasks = res_func.funcComp.tasks(isDiff);

%% Count and Sort Passed Tasks
pass_counts = sum(res_func.funcComp.matrix(isDiff, :), 1)';
[sorted_counts, sort_idx] = sort(pass_counts, 'descend');
tissue_sorted = useModels(sort_idx);
zscore_counts = normalize(sorted_counts);

% Save Task Summary Table
T = table(tissue_sorted, sorted_counts, zscore_counts, ...
    'VariableNames', {'Tissue', 'y2', 'zscore'});
writetable(T, fullfile(subfolder, 'availabel_fns_all.xlsx'));
save(fullfile(subfolder, 'T.mat'), 'T');

%% Plot Binary Task Completion Matrix
out = figure;
spy(res_func.funcComp.matrix(isDiff, :), 30);
set(gca, ...
    'XTick', 1:numel(useModels), ...
    'XTickLabel', useModels, ...
    'XTickLabelRotation', 90, ...
    'YTick', 1:numel(diffTasks), ...
    'YTickLabel', diffTasks, ...
    'YAxisLocation', 'right');
xlabel('');
title('Differentially Performed Metabolic Tasks');
legend('Pass', 'Location', 'southeast');
legend('boxoff');

if save_in == 1
    saveas(gcf, fullfile(subfolder, 'task.png'));
    saveas(gcf, fullfile(subfolder, 'task.fig'));
end

%% Save Matrix to CSV
matrix = res_func.funcComp.matrix(isDiff, :);
output_matrix = array2table(matrix, ...
    'VariableNames', useModels, ...
    'RowNames', diffTasks);

writetable(output_matrix, fullfile(subfolder, 'functional_comparison_matrix.csv'), ...
    'WriteRowNames', true);
disp(['Saved functional matrix to: ' fullfile(subfolder, 'functional_comparison_matrix.csv')]);

end
