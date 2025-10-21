%% script_04_compare_all
% Extract all up/down regulated pathways from pFBA and FVA-sampling

clc;
% clear all;

%% 1. Set up
pathway = pwd;
results_dir = fullfile(pathway, '4_compare_all');

if ~exist(results_dir, 'dir')
    mkdir(results_dir); 
end

%% 2. Load all three method results
fprintf('=== LOADING ALL THREE METHOD RESULTS ===\n');

% Load pFBA results
pFBA_file = fullfile(pathway, '2_pFBA_PFI', 'pathway_flux_analysis.xlsx');
pFBA_data = readtable(pFBA_file, 'Sheet', 'Mean_HSD_vs_NSD', 'ReadRowNames', true);

% Load FVA-sampling results
FVA_file = fullfile(pathway, '3_sampling_PFI', 'pathway_flux_analysis.xlsx');
FVA_data = readtable(FVA_file, 'Sheet', 'Mean_HSD_vs_NSD', 'ReadRowNames', true);

fprintf('pFBA data: %d pathways loaded\n', height(pFBA_data));
fprintf('FVA-sampling data: %d pathways loaded\n', height(FVA_data));

%% 3. Define thresholds
up_threshold = 0;    % Log2FC > 0 (any increase)
down_threshold = 0;  % Log2FC < 0 (any decrease)

fprintf('\nThresholds:\n');
fprintf('  Up-regulated: Log2FC > %.3f\n', up_threshold);
fprintf('  Down-regulated: Log2FC < %.3f\n', down_threshold);

%% 4. Extract up/down pathways from each method
fprintf('\n=== EXTRACTING UP/DOWN PATHWAYS ===\n');

% pFBA
pFBA_pathways = pFBA_data.Properties.RowNames;
pFBA_log2fc = pFBA_data.Log2FC;
pFBA_up = pFBA_pathways(pFBA_log2fc > up_threshold & ~isnan(pFBA_log2fc) & ~isinf(pFBA_log2fc));
pFBA_down = pFBA_pathways(pFBA_log2fc < down_threshold & ~isnan(pFBA_log2fc) & ~isinf(pFBA_log2fc));

fprintf('pFBA:\n');
fprintf('  Up-regulated: %d pathways\n', length(pFBA_up));
fprintf('  Down-regulated: %d pathways\n', length(pFBA_down));

% FVA-sampling
FVA_pathways = FVA_data.Properties.RowNames;
FVA_log2fc = FVA_data.Log2FC;
FVA_up = FVA_pathways(FVA_log2fc > up_threshold & ~isnan(FVA_log2fc) & ~isinf(FVA_log2fc));
FVA_down = FVA_pathways(FVA_log2fc < down_threshold & ~isnan(FVA_log2fc) & ~isinf(FVA_log2fc));

fprintf('\nFVA-sampling:\n');
fprintf('  Up-regulated: %d pathways\n', length(FVA_up));
fprintf('  Down-regulated: %d pathways\n', length(FVA_down));

%% 5. Get all unique pathways for each direction
all_up_pathways = unique([pFBA_up; FVA_up]);
all_down_pathways = unique([pFBA_down; FVA_down]);

fprintf('\n=== TOTAL UNIQUE PATHWAYS ===\n');
fprintf('All up-regulated pathways: %d\n', length(all_up_pathways));
fprintf('All down-regulated pathways: %d\n', length(all_down_pathways));

%% 6. Create binary matrices for UpSet plot
% Up-regulated pathways
upset_up = table();
upset_up.Pathway = all_up_pathways;
upset_up.pFBA = double(ismember(all_up_pathways, pFBA_up));
upset_up.FVA_sampling = double(ismember(all_up_pathways, FVA_up));

% Down-regulated pathways
upset_down = table();
upset_down.Pathway = all_down_pathways;
upset_down.pFBA = double(ismember(all_down_pathways, pFBA_down));
upset_down.FVA_sampling = double(ismember(all_down_pathways, FVA_down));

%% 7. Export to Excel
fprintf('\n=== EXPORTING TO EXCEL ===\n');

excel_filename = fullfile(results_dir, 'UpSet_pathway_data.xlsx');

% Export binary matrices for UpSet plots
writetable(upset_up, excel_filename, 'Sheet', 'Up_Regulated', 'WriteRowNames', false);
writetable(upset_down, excel_filename, 'Sheet', 'Down_Regulated', 'WriteRowNames', false);


% Also export the pathway lists separately
max_len_up = max([length(pFBA_up), length(FVA_up)]);
max_len_down = max([length(pFBA_down), length(FVA_down)]);

% Pad lists to same length
up_lists = table();
up_lists.pFBA_Up = [pFBA_up; repmat({''}, max_len_up - length(pFBA_up), 1)];
up_lists.FVA_Up = [FVA_up; repmat({''}, max_len_up - length(FVA_up), 1)];

down_lists = table();
down_lists.pFBA_Down = [pFBA_down; repmat({''}, max_len_down - length(pFBA_down), 1)];
down_lists.FVA_Down = [FVA_down; repmat({''}, max_len_down - length(FVA_down), 1)];

writetable(up_lists, excel_filename, 'Sheet', 'Up_Lists', 'WriteRowNames', false);
writetable(down_lists, excel_filename, 'Sheet', 'Down_Lists', 'WriteRowNames', false);

% Save MATLAB workspace
save(fullfile(results_dir, 'upset_data.mat'), ...
     'upset_up', 'upset_down', 'pFBA_up', 'pFBA_down', ...
     'FVA_up', 'FVA_down', ...
     'all_up_pathways', 'all_down_pathways');

fprintf('\n✅ Data exported to: %s\n', excel_filename);
fprintf('   - Up_Regulated sheet: binary matrix for UpSet plot\n');
fprintf('   - Down_Regulated sheet: binary matrix for UpSet plot\n');
fprintf('   - Up_Lists sheet: pathway lists by method\n');
fprintf('   - Down_Lists sheet: pathway lists by method\n');
fprintf('\n=== EXPORT COMPLETE ===\n');