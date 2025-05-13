%% Fisher's Exact Test for Functional Task Enrichment in Tissue-Specific GEMs

clc;
clear;

%% Define Save Directory
savePath = '6_analysis_fisher_test';
saveFolder = fullfile(pwd, savePath);
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end

%% Load Input Data

% Reference: SYSTEM & SUBSYSTEM info
refFile = fullfile(pwd, 'files', 'metabolic_task_Full.xlsx');
df_ref = readtable(refFile, 'ReadVariableNames', true);

% Functional task completion matrix
funcMatrixFile = fullfile(pwd, '5_analysis_fn_metabolic_task', 'functional_comparison_matrix.csv');
df = readtable(funcMatrixFile, 'ReadVariableNames', true);
df.Properties.VariableNames{1} = 'Description';

% Initialize SYSTEM and SUBSYSTEM columns
df.SYSTEM = repmat({''}, height(df), 1);
df.SUBSYSTEM = repmat({''}, height(df), 1);

% Map SYSTEM and SUBSYSTEM to main table
for i = 1:height(df)
    desc = df.Description{i};
    match_idx = find(strcmpi(df_ref.DESCRIPTION, desc));
    
    if ~isempty(match_idx)
        df.SYSTEM{i}    = df_ref.SYSTEM{match_idx};
        df.SUBSYSTEM{i} = df_ref.SUBSYSTEM{match_idx};
    else
        df.SYSTEM{i}    = 'Not Found';
        df.SUBSYSTEM{i} = 'Not Found';
    end
end

% Reorder SYSTEM and SUBSYSTEM to appear first
column_order = [{'SYSTEM', 'SUBSYSTEM'}, df.Properties.VariableNames(1:end-2)];
df = df(:, column_order);

% Save updated task matrix with SYSTEM & SUBSYSTEM
writetable(df, fullfile(saveFolder, 'functional_comparison_matrix_u.csv'));
disp('✅ Functional matrix with SYSTEM/SUBSYSTEM info saved.');

%% Task Count Summary by SYSTEM and SUBSYSTEM

tissue_columns = df.Properties.VariableNames(4:end);

% SYSTEM task counts
system_summary = groupsummary(df, 'SYSTEM', @(x) sum(x == 1), tissue_columns);
system_summary.Properties.VariableNames(3:end) = tissue_columns;
writetable(system_summary, fullfile(saveFolder, 'tasks_per_system.csv'));

% SUBSYSTEM task counts
subsystem_summary = groupsummary(df, 'SUBSYSTEM', @(x) sum(x == 1), tissue_columns);
subsystem_summary.Properties.VariableNames(3:end) = tissue_columns;
writetable(subsystem_summary, fullfile(saveFolder, 'tasks_per_subsystem.csv'));

%% Reference Summary Counts

% SYSTEM reference summary
unique_systems = unique(df_ref.SYSTEM);
unique_systems = unique_systems(~cellfun('isempty', unique_systems));
ref_counts = arrayfun(@(s) sum(strcmp(df_ref.SYSTEM, s)), unique_systems);
ref_summary = table(unique_systems, ref_counts, 'VariableNames', {'SYSTEM', 'TotalTasks'});
writetable(ref_summary, fullfile(saveFolder, 'tasks_per_system_reference.csv'));

% SUBSYSTEM reference summary
unique_subsystems = unique(df_ref.SUBSYSTEM);
unique_subsystems = unique_subsystems(~cellfun('isempty', unique_subsystems));
ref_counts_sub = arrayfun(@(s) sum(strcmp(df_ref.SUBSYSTEM, s)), unique_subsystems);
ref_summary_subsystem = table(unique_subsystems, ref_counts_sub, 'VariableNames', {'SUBSYSTEM', 'TotalTasks'});
writetable(ref_summary_subsystem, fullfile(saveFolder, 'tasks_per_subsystem_reference.csv'));

%% Fisher’s Exact Test (SYSTEM)

fisher_results = table(system_summary.SYSTEM, 'VariableNames', {'SYSTEM'});
odd_ratio_results = table(system_summary.SYSTEM, 'VariableNames', {'SYSTEM'});

for col = 1:numel(tissue_columns)
    tissue = tissue_columns{col};
    p_values = zeros(height(system_summary), 1);
    odds_ratios = zeros(height(system_summary), 1);

    for i = 1:height(system_summary)
        tissue_count = system_summary{i, tissue};
        ref_count = ref_summary.TotalTasks(strcmp(ref_summary.SYSTEM, system_summary.SYSTEM{i}));
        total_in_tissue = sum(system_summary{:, tissue});
        total_in_ref = sum(ref_summary.TotalTasks);

        % Build contingency table
        contingency = [
            tissue_count,                 ref_count - tissue_count;
            total_in_tissue - tissue_count, total_in_ref - ref_count - (total_in_tissue - tissue_count)
        ];

        [~, p_values(i)] = fishertest(contingency);
        odds_ratios(i) = (contingency(1,1) * contingency(2,2)) / max((contingency(1,2) * contingency(2,1)), 1); % Avoid divide-by-zero
    end

    fisher_results.(tissue) = p_values;
    odd_ratio_results.(tissue) = odds_ratios;
end

writetable(fisher_results, fullfile(saveFolder, 'fisher_test_results.csv'));
writetable(odd_ratio_results, fullfile(saveFolder, 'fisher_test_results_oddR.csv'));

%% Adjust SYSTEM p-values (Benjamini-Hochberg FDR)
adjusted_p = zeros(height(fisher_results), numel(tissue_columns));
for col = 1:numel(tissue_columns)
    adjusted_p(:, col) = mafdr(fisher_results{:, tissue_columns{col}}, 'BHFDR', true);
end

adjusted_table = array2table(adjusted_p, 'VariableNames', tissue_columns);
adjusted_table.SYSTEM = fisher_results.SYSTEM;
adjusted_table = [adjusted_table(:, end), adjusted_table(:, 1:end-1)];
writetable(adjusted_table, fullfile(saveFolder, 'fisher_test_results_adjusted.csv'));

%% Fisher’s Exact Test (SUBSYSTEM)

fisher_sub = table(subsystem_summary.SUBSYSTEM, 'VariableNames', {'SUBSYSTEM'});
odds_sub = table(subsystem_summary.SUBSYSTEM, 'VariableNames', {'SUBSYSTEM'});

for col = 1:numel(tissue_columns)
    tissue = tissue_columns{col};
    p_vals = zeros(height(subsystem_summary), 1);
    odds_vals = zeros(height(subsystem_summary), 1);

    for i = 1:height(subsystem_summary)
        tissue_count = subsystem_summary{i, tissue};
        ref_count = ref_summary_subsystem.TotalTasks(strcmp(ref_summary_subsystem.SUBSYSTEM, subsystem_summary.SUBSYSTEM{i}));
        total_tissue = sum(subsystem_summary{:, tissue});
        total_ref = sum(ref_summary_subsystem.TotalTasks);

        contingency = [
            tissue_count, ref_count - tissue_count;
            total_tissue - tissue_count, total_ref - ref_count
        ];

        [~, p_vals(i)] = fishertest(contingency);
        odds_vals(i) = (contingency(1,1) * contingency(2,2)) / max((contingency(1,2) * contingency(2,1)), 1); % Avoid divide-by-zero
    end

    fisher_sub.(tissue) = p_vals;
    odds_sub.(tissue) = odds_vals;
end

writetable(fisher_sub, fullfile(saveFolder, 'fisher_subsystem.csv'));
writetable(odds_sub, fullfile(saveFolder, 'fisher_subsystem_oddR.csv'));

disp('✅ All enrichment results saved.');
