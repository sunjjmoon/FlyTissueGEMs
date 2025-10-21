%% script_02_PFI_pFBA

clc;
close all
% clear all;

%% Initialize COBRA Toolbox
% initCobraToolbox(false);
% changeCobraSolver('gurobi', 'LP'); % commented out like in original code


%% 1. Set up
pathway = pwd;
results_dir = fullfile(pathway, '2_pFBA_PFI');

if ~exist(results_dir, 'dir')
    mkdir(results_dir); 
end

%% 2. Load files
% Load the constrained models (has pathway info in subSystems)
load(fullfile(pathway, '01_results', 'model_constrained_out.mat'));

% Load the pFBA results (has the flux values)
load(fullfile(pathway, '01_results', 'pFBA_results.mat'));

fprintf('Loaded %d models from constrained models\n', length(model_constrained_out));
fprintf('Loaded %d pFBA results\n', length(pFBA_results));

% Quick check that we have matching data
for i = 1:length(model_constrained_out)
    model_id_constrained = model_constrained_out{i,1}.modelID;
    model_id_pFBA = pFBA_results{i,1}.modelID;
    
    if ~strcmp(model_id_constrained, model_id_pFBA)
        error('Model IDs do not match at index %d: %s vs %s', ...
            i, model_id_constrained, model_id_pFBA);
    end
end

fprintf('All model IDs match between constrained models and pFBA results\n');

%% 3-4. Extract pathway-level flux summaries and create comparison matrices
pathway_flux_summary = cell(length(model_constrained_out), 1);

all_pathways = {};
for i = 1:length(model_constrained_out)
    model = model_constrained_out{i,1};
    unique_pathways = unique(string(model.subSystems));
    unique_pathways = unique_pathways(~cellfun(@isempty, unique_pathways));
    all_pathways = [all_pathways; cellstr(unique_pathways)];
end
all_pathways = unique(all_pathways);
fprintf('Found %d unique pathways across all models\n', length(all_pathways));

% Initialize comparison matrices
model_ids = cell(length(model_constrained_out), 1);
mean_flux_matrix = zeros(length(all_pathways), length(model_constrained_out));
median_flux_matrix = zeros(length(all_pathways), length(model_constrained_out));
q25_flux_matrix = zeros(length(all_pathways), length(model_constrained_out));
q75_flux_matrix = zeros(length(all_pathways), length(model_constrained_out));
std_flux_matrix = zeros(length(all_pathways), length(model_constrained_out));
max_flux_matrix = zeros(length(all_pathways), length(model_constrained_out));
min_flux_matrix = zeros(length(all_pathways), length(model_constrained_out));
total_flux_matrix = zeros(length(all_pathways), length(model_constrained_out));
reaction_count_matrix = zeros(length(all_pathways), length(model_constrained_out));

% Second pass: extract statistics for each model
fprintf('\n=== EXTRACTING PATHWAY STATISTICS ===\n');
for i = 1:length(model_constrained_out)
    model = model_constrained_out{i,1};
    fluxes = pFBA_results{i,1}.flux;
    model_id = model.modelID;
    model_ids{i} = model_id;
    
    fprintf('Processing model %s (%d/%d)\n', model_id, i, length(model_constrained_out));
    
    % Initialize results for this model
    pathway_results = struct();
    pathway_results.modelID = model_id;
    pathway_results.pathway_names = all_pathways;
    pathway_results.stats = struct();
    
    % Calculate statistics for each pathway
    subsys_tmp = string(model.subSystems);
    
    for p = 1:length(all_pathways)
        pathway_name = all_pathways{p};
        
        % Find reactions in this pathway
        rxn_indices = find(strcmp(subsys_tmp, pathway_name));
        
        if ~isempty(rxn_indices)
            % Get absolute fluxes for reactions in this pathway
            pathway_fluxes = abs(fluxes(rxn_indices));
            
            % Calculate comprehensive statistics
            mean_flux_matrix(p, i) = mean(pathway_fluxes);
            median_flux_matrix(p, i) = median(pathway_fluxes);
            std_flux_matrix(p, i) = std(pathway_fluxes);
            max_flux_matrix(p, i) = max(pathway_fluxes);
            min_flux_matrix(p, i) = min(pathway_fluxes);
            total_flux_matrix(p, i) = sum(pathway_fluxes);
            reaction_count_matrix(p, i) = length(rxn_indices);
            
            % Percentiles
            q25_flux_matrix(p, i) = prctile(pathway_fluxes, 25);
            q75_flux_matrix(p, i) = prctile(pathway_fluxes, 75);
            
            % Store in pathway_results structure
            pathway_results.stats.(matlab.lang.makeValidName(pathway_name)) = struct(...
                'mean', mean_flux_matrix(p, i), ...
                'median', median_flux_matrix(p, i), ...
                'std', std_flux_matrix(p, i), ...
                'min', min_flux_matrix(p, i), ...
                'max', max_flux_matrix(p, i), ...
                'q25', q25_flux_matrix(p, i), ...
                'q75', q75_flux_matrix(p, i), ...
                'total', total_flux_matrix(p, i), ...
                'reaction_count', reaction_count_matrix(p, i));
        end
        % If pathway not in model, matrices remain 0 (already initialized)
    end
    
    % Store results
    pathway_flux_summary{i,1} = pathway_results;
end

% Find HSD and NSD model indices
hsd_idx = find(strcmp(model_ids, 'HSD'));
nsd_idx = find(strcmp(model_ids, 'NSD'));

if isempty(hsd_idx) || isempty(nsd_idx)
    error('Could not find both HSD and NSD models');
end


% Get reaction counts for NSD and HSD
nsd_rxn_count = reaction_count_matrix(:, nsd_idx);
hsd_rxn_count = reaction_count_matrix(:, hsd_idx);

% MEAN comparisons (filtered)
mean_nsd = mean_flux_matrix(:, nsd_idx);
mean_hsd = mean_flux_matrix(:, hsd_idx);
% Filter: only include pathways where at least one condition has non-zero flux
mean_valid_idx = (mean_nsd > 0) | (mean_hsd > 0);
mean_nsd_filtered = mean_nsd;
mean_hsd_filtered = mean_hsd;
mean_nsd_filtered(~mean_valid_idx) = NaN;
mean_hsd_filtered(~mean_valid_idx) = NaN;
mean_ratio = mean_hsd_filtered ./ mean_nsd_filtered;
mean_log2fc = log2(mean_ratio);

% MEDIAN comparisons (filtered)
median_nsd = median_flux_matrix(:, nsd_idx);
median_hsd = median_flux_matrix(:, hsd_idx);
median_valid_idx = (median_nsd > 0) | (median_hsd > 0);
median_nsd_filtered = median_nsd;
median_hsd_filtered = median_hsd;
median_nsd_filtered(~median_valid_idx) = NaN;
median_hsd_filtered(~median_valid_idx) = NaN;
median_ratio = median_hsd_filtered ./ median_nsd_filtered;
median_log2fc = log2(median_ratio);

% MODE comparisons (filtered, using Q25 as proxy)
mode_nsd = q25_flux_matrix(:, nsd_idx);
mode_hsd = q25_flux_matrix(:, hsd_idx);
mode_valid_idx = (mode_nsd > 0) | (mode_hsd > 0);
mode_nsd_filtered = mode_nsd;
mode_hsd_filtered = mode_hsd;
mode_nsd_filtered(~mode_valid_idx) = NaN;
mode_hsd_filtered(~mode_valid_idx) = NaN;
mode_ratio = mode_hsd_filtered ./ mode_nsd_filtered;
mode_log2fc = log2(mode_ratio);

% Show filtering results
fprintf('\n=== FILTERING RESULTS ===\n');
fprintf('Mean   : %d pathways remain after filtering\n', sum(mean_valid_idx));
fprintf('Median : %d pathways remain after filtering\n', sum(median_valid_idx));
fprintf('Mode   : %d pathways remain after filtering\n', sum(mode_valid_idx));


% 
% % MEAN comparisons
% mean_nsd = mean_flux_matrix(:, nsd_idx);
% mean_hsd = mean_flux_matrix(:, hsd_idx);
% mean_ratio = mean_hsd ./ (mean_nsd );
% mean_log2fc = log2(mean_ratio);
% 
% % MEDIAN comparisons  
% median_nsd = median_flux_matrix(:, nsd_idx);
% median_hsd = median_flux_matrix(:, hsd_idx);
% median_ratio = median_hsd ./ (median_nsd);
% median_log2fc = log2(median_ratio);
% 
% % MODE comparisons (using Q25 as proxy for mode since mode might not be meaningful for continuous flux)
% mode_nsd = q25_flux_matrix(:, nsd_idx);
% mode_hsd = q25_flux_matrix(:, hsd_idx);
% mode_ratio = mode_hsd ./ (mode_nsd);
% mode_log2fc = log2(mode_ratio);

% Create comparison tables
% Create comparison tables WITH REACTION COUNTS
mean_comparison_table = table(mean_nsd, mean_hsd, mean_ratio, mean_log2fc, nsd_rxn_count, hsd_rxn_count, ...
    'VariableNames', {'NSD', 'HSD', 'HSD_NSD_Ratio', 'Log2FC', 'NSD_n', 'HSD_n'}, ...
    'RowNames', all_pathways);

median_comparison_table = table(median_nsd, median_hsd, median_ratio, median_log2fc, nsd_rxn_count, hsd_rxn_count, ...
    'VariableNames', {'NSD', 'HSD', 'HSD_NSD_Ratio', 'Log2FC', 'NSD_n', 'HSD_n'}, ...
    'RowNames', all_pathways);

mode_comparison_table = table(mode_nsd, mode_hsd, mode_ratio, mode_log2fc, nsd_rxn_count, hsd_rxn_count, ...
    'VariableNames', {'NSD', 'HSD', 'HSD_NSD_Ratio', 'Log2FC', 'NSD_n', 'HSD_n'}, ...
    'RowNames', all_pathways);

% Create readable tables (original tables)
mean_flux_table = array2table(mean_flux_matrix, 'VariableNames', model_ids, 'RowNames', all_pathways);
median_flux_table = array2table(median_flux_matrix, 'VariableNames', model_ids, 'RowNames', all_pathways);


%% Save results to Excel
fprintf('\n=== SAVING RESULTS TO EXCEL ===\n');
excel_filename = fullfile(results_dir, 'pathway_flux_analysis.xlsx');


%% Save results to Excel
fprintf('\n=== SAVING RESULTS TO EXCEL ===\n');
excel_filename = fullfile(results_dir, 'pathway_flux_analysis.xlsx');

% Get reaction counts for NSD and HSD
nsd_rxn_count = reaction_count_matrix(:, nsd_idx);
hsd_rxn_count = reaction_count_matrix(:, hsd_idx);

% Save original data tables
writetable(mean_flux_table, excel_filename, 'Sheet', 'Mean_Flux_Raw', 'WriteRowNames', true);
writetable(median_flux_table, excel_filename, 'Sheet', 'Median_Flux_Raw', 'WriteRowNames', true);

% Save individual comparison vectors as tables WITH REACTION COUNTS
mean_components = table(mean_nsd, mean_hsd, mean_ratio, mean_log2fc, nsd_rxn_count, hsd_rxn_count, ...
    'VariableNames', {'NSD', 'HSD', 'HSD_NSD_Ratio', 'Log2FC', 'NSD_n', 'HSD_n'}, ...
    'RowNames', all_pathways);
writetable(mean_components, excel_filename, 'Sheet', 'Mean_All_Pathways', 'WriteRowNames', true);

median_components = table(median_nsd, median_hsd, median_ratio, median_log2fc, nsd_rxn_count, hsd_rxn_count, ...
    'VariableNames', {'NSD', 'HSD', 'HSD_NSD_Ratio', 'Log2FC', 'NSD_n', 'HSD_n'}, ...
    'RowNames', all_pathways);
writetable(median_components, excel_filename, 'Sheet', 'Median_All_Pathways', 'WriteRowNames', true);

mode_components = table(mode_nsd, mode_hsd, mode_ratio, mode_log2fc, nsd_rxn_count, hsd_rxn_count, ...
    'VariableNames', {'NSD', 'HSD', 'HSD_NSD_Ratio', 'Log2FC', 'NSD_n', 'HSD_n'}, ...
    'RowNames', all_pathways);
writetable(mode_components, excel_filename, 'Sheet', 'Mode_All_Pathways', 'WriteRowNames', true);

% Save MATLAB results
% save(fullfile(results_dir, 'pathway_flux_analysis.mat'), ...
%      'mean_comparison_table', 'median_comparison_table', 'mode_comparison_table', ...
%      'mean_flux_table', 'median_flux_table', 'all_pathways', 'model_ids');

 % Save MATLAB results
% Save MATLAB results
save(fullfile(results_dir, 'pathway_flux_analysis.mat'), ...
     'mean_comparison_sorted', 'median_comparison_sorted', 'mode_comparison_sorted', ...
     'mean_comparison_table', 'median_comparison_table', 'mode_comparison_table', ...
     'mean_nsd', 'mean_hsd', 'mean_ratio', 'mean_log2fc', ...
     'median_nsd', 'median_hsd', 'median_ratio', 'median_log2fc', ...
     'mode_nsd', 'mode_hsd', 'mode_ratio', 'mode_log2fc', ...
     'mean_flux_table', 'median_flux_table', ...
     'all_pathways', 'model_ids');

% Mean
valid_mean = sum(~isnan(mean_log2fc) & ~isinf(mean_log2fc) & mean_log2fc ~= 0);
fprintf('Mean Log2FC: %d valid entries\n', valid_mean);

% Median
valid_median = sum(~isnan(median_log2fc) & ~isinf(median_log2fc) & median_log2fc ~= 0);
fprintf('Median Log2FC: %d valid entries\n', valid_median);

% Mode (Q25 proxy)
valid_mode = sum(~isnan(mode_log2fc) & ~isinf(mode_log2fc) & mode_log2fc ~= 0);
fprintf('Mode Log2FC: %d valid entries\n', valid_mode);

total_paths = numel(all_pathways);

fprintf('\n=== VALID Log2FC COUNTS ===\n');
fprintf('Mean   : %d / %d pathways\n', valid_mean, total_paths);
fprintf('Median : %d / %d pathways\n', valid_median, total_paths);
fprintf('Mode   : %d / %d pathways\n', valid_mode, total_paths);





%% Create bar graph for top 5 and bottom 5 Log2FC pathways
fprintf('\n=== CREATING BAR GRAPH ===\n');

% Get valid mean data (sorted by Log2FC, not absolute)
valid_idx = ~isnan(mean_log2fc) & ~isinf(mean_log2fc);
valid_pathways = all_pathways(valid_idx);
valid_log2fc = mean_log2fc(valid_idx);

% Sort by actual Log2FC (highest to lowest)
[sorted_log2fc, sort_idx] = sort(valid_log2fc, 'descend');
sorted_pathways = valid_pathways(sort_idx);

% Get top 5 and bottom 5
top5_pathways = sorted_pathways(1:5);
top5_log2fc = sorted_log2fc(1:5);

bottom5_pathways = sorted_pathways(end-4:end);
bottom5_log2fc = sorted_log2fc(end-4:end);

% Combine for plotting
plot_pathways = [top5_pathways; bottom5_pathways];
plot_log2fc = [top5_log2fc; bottom5_log2fc];

% Create figure
figure(1)
set(gcf, 'WindowState', 'maximized');   % 

% Create bar plot
b = bar(plot_log2fc, 'FaceColor', 'flat');

% Color bars: red for positive (upregulated), blue for negative (downregulated)
colors = [plot_log2fc > 0] * [1 0.2 0.2] + [plot_log2fc <= 0] * [0.2 0.4 0.8];
b.CData = colors;

% Customize plot
% xlabel('Pathways', 'FontSize', 12, 'FontWeight', 'bold');
ylabel(['Log2 FC (HSD/NSD)'], 'FontSize', 12, 'FontWeight', 'bold');
% title('Top 5 Up-regulated and Top 5 Down-regulated Pathways', 'FontSize', 14, 'FontWeight', 'bold');

% Add pathway names as x-axis labels (rotated)
xticks(1:length(plot_pathways));
xticklabels(plot_pathways);
xtickangle(45);

% Add horizontal line at y=0
hold on;
yline(0, 'k--', 'LineWidth', 1);

% Add value labels on bars
for i = 1:length(plot_log2fc)
    if plot_log2fc(i) >= 0
        text(i, plot_log2fc(i) + 0.1, sprintf('%.2f', plot_log2fc(i)), ...
             'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
    else
        text(i, plot_log2fc(i) - 0.2, sprintf('%.2f', plot_log2fc(i)), ...
             'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
    end
end

% Improve layout
grid on;
grid minor;
set(gca, 'FontSize', 11);

% Save figure
saveas(gcf, fullfile(results_dir, 'Top_Bottom_Pathways_Log2FC.png'));
saveas(gcf, fullfile(results_dir, 'Top_Bottom_Pathways_Log2FC.fig'));

fprintf('Bar graph saved: Top_Bottom_Pathways_Log2FC.png/.fig\n');


%% 7. Analyze reactions within top 5 and bottom 5 pathways
fprintf('\n=== ANALYZING REACTIONS IN TOP/BOTTOM PATHWAYS ===\n');

% Get the pathway names we're interested in
pathways_of_interest = [top5_pathways; bottom5_pathways];
fprintf('Analyzing reactions in %d pathways\n', length(pathways_of_interest));

% Initialize results
reaction_analysis = {};
row_count = 1;

% Add headers
reaction_analysis{row_count, 1} = 'Pathway';
reaction_analysis{row_count, 2} = 'Reaction_ID';
reaction_analysis{row_count, 3} = 'Reaction_Name';
reaction_analysis{row_count, 4} = 'NSD_Flux';
reaction_analysis{row_count, 5} = 'HSD_Flux';
reaction_analysis{row_count, 6} = 'HSD_NSD_Ratio';
reaction_analysis{row_count, 7} = 'Log2FC';
row_count = row_count + 1;

% Analyze each pathway
for p = 1:length(pathways_of_interest)
    pathway_name = pathways_of_interest{p};
    fprintf('Processing pathway: %s\n', pathway_name);
    
    % Find reactions in this pathway (use HSD model as reference)
    model = model_constrained_out{hsd_idx, 1};
    subsys_tmp = string(model.subSystems);
    rxn_indices = find(strcmp(subsys_tmp, pathway_name));
    
    % Get fluxes for both models
    nsd_fluxes = pFBA_results{nsd_idx, 1}.flux;
    hsd_fluxes = pFBA_results{hsd_idx, 1}.flux;
    
    for r = 1:length(rxn_indices)
        rxn_idx = rxn_indices(r);
        
        % Get reaction info
        rxn_id = model.rxns{rxn_idx};
        rxn_name = model.rxnNames{rxn_idx};
        
        % Get fluxes (use absolute values)
        nsd_flux = abs(nsd_fluxes(rxn_idx));
        hsd_flux = abs(hsd_fluxes(rxn_idx));
        
        % Calculate ratio and log2fc
        flux_ratio = hsd_flux / (nsd_flux + eps);
        flux_log2fc = log2(flux_ratio);
        
        % Store results
        reaction_analysis{row_count, 1} = pathway_name;
        reaction_analysis{row_count, 2} = rxn_id;
        reaction_analysis{row_count, 3} = rxn_name;
        reaction_analysis{row_count, 4} = nsd_flux;
        reaction_analysis{row_count, 5} = hsd_flux;
        reaction_analysis{row_count, 6} = flux_ratio;
        reaction_analysis{row_count, 7} = flux_log2fc;
        row_count = row_count + 1;
    end
end

% Convert to table
reaction_table = cell2table(reaction_analysis(2:end, :), 'VariableNames', reaction_analysis(1, :));

% Remove rows with Inf or -Inf in Log2FC
% reaction_table_noInf = reaction_table(~isinf(reaction_table.Log2FC), :);

% reaction_table.Log2FC(isinf(reaction_table.Log2FC)) = NaN;

% Now sort by Pathway (ascending) and Log2FC (descending)
% reaction_table_sorted = sortrows(reaction_table_noInf, ...
%                                  {'Pathway','Log2FC'}, ...
%                                  {'ascend','descend'});
% 
                             
% Sort by absolute Log2FC within each pathway
% reaction_table_sorted = sortrows(reaction_table, {'Pathway', 'Log2FC'}, {'ascend', 'descend'});

% Save to Excel
reaction_excel_filename = fullfile(results_dir, 'reaction_analysis_top_bottom_pathways.xlsx');

writetable(reaction_table, reaction_excel_filename, 'Sheet', 'Reaction_Analysis', 'WriteRowNames', false);

% Save to MATLAB
save(fullfile(results_dir, 'pathway_flux_analysis.mat'), 'reaction_table', '-append');

fprintf('Reaction analysis completed and saved to Excel sheet: Reaction_Analysis\n');
fprintf('Total reactions analyzed: %d\n', height(reaction_table));


%% 8. Create grouped bar chart with error bars for top/bottom pathways

% Get std values for the top/bottom pathways
top5_nsd_mean = mean_flux_matrix(sort_idx(1:5), nsd_idx);
top5_hsd_mean = mean_flux_matrix(sort_idx(1:5), hsd_idx);
top5_nsd_std = std_flux_matrix(sort_idx(1:5), nsd_idx);
top5_hsd_std = std_flux_matrix(sort_idx(1:5), hsd_idx);

bottom5_nsd_mean = mean_flux_matrix(sort_idx(end-4:end), nsd_idx);
bottom5_hsd_mean = mean_flux_matrix(sort_idx(end-4:end), hsd_idx);
bottom5_nsd_std = std_flux_matrix(sort_idx(end-4:end), nsd_idx);
bottom5_hsd_std = std_flux_matrix(sort_idx(end-4:end), hsd_idx);

% Combine data for plotting
plot_nsd_mean = [top5_nsd_mean; bottom5_nsd_mean];
plot_hsd_mean = [top5_hsd_mean; bottom5_hsd_mean];
plot_nsd_std = [top5_nsd_std; bottom5_nsd_std];
plot_hsd_std = [top5_hsd_std; bottom5_hsd_std];

% Create grouped bar data
bar_data = [plot_nsd_mean, plot_hsd_mean];
error_data = [plot_nsd_std, plot_hsd_std];

% Create figure
figure(2);
set(gcf, 'WindowState', 'maximized');   % works in R2018b+

% figure('Position', [150, 150, 900, 600]);

% Create grouped bar plot with error bars
b = bar(bar_data, 'grouped');
hold on;

% Add error bars
ngroups = size(bar_data, 1);
nbars = size(bar_data, 2);
groupwidth = min(0.8, nbars/(nbars + 1.5));

for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    errorbar(x, bar_data(:,i), error_data(:,i), 'k.', 'LineWidth', 1.5);
end

% Customize colors
% b(1).FaceColor = [0.3 0.6 0.9]; % Blue for NSD
% b(2).FaceColor = [0.9 0.4 0.3]; % Red for HSD

% Customize plot
% xlabel('Pathways', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Mean flux (a.u.)', 'FontSize', 14, 'FontWeight', 'bold');
% title('Pathway Activity Comparison: HSD vs NSD (Top 5 Up & Top 5 Down)', 'FontSize', 14, 'FontWeight', 'bold');

% Add pathway names as x-axis labels
xticks(1:length(plot_pathways));
xticklabels(plot_pathways);
xtickangle(45);

% Add legend
legend({'NSD', 'HSD'}, 'Location', 'northeast', 'FontSize', 11);

% Improve layout
grid on;
grid minor;
% set(gca, 'FontSize', 10);

% Save figure
saveas(gcf, fullfile(results_dir, 'Grouped_Bar_Mean_SD_Pathways.png'));
saveas(gcf, fullfile(results_dir, 'Grouped_Bar_Mean_SD_Pathways.fig'));

