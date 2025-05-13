%% =================== Extended PFI Visualization Script (v3, Fixed) ===================
% Purpose: Visualize Pathway Flux Index (PFI) with SEM error bars, removing the top artifact.
% ================================================================

clearvars; close all; clc;

%% ------------------ Create Output Folder ------------------
% Define main folder paths
pathway = pwd;
output_folder = fullfile(pathway, 'G_PFVI_visualize');

% Create folder if it doesn't exist
if ~exist(output_folder, 'dir')
    mkdir(output_folder);
end

%% ------------------ Load PFI Data ------------------
% Load mean PFI values and SEM values from F_PFI
mean_data = readtable(fullfile(pathway, 'F_PFVI_mode', 'PFVI_results.xlsx'));  
sem_data = readtable(fullfile(pathway, 'F_PFVI_mode', 'PFVI_results_sem.xlsx'));  
count_data = readtable(fullfile(pathway, 'F_PFVI_mode', 'PFVI_results_count.xlsx'));  

% Extract pathway names
pathway_names = mean_data{:, 1}; % First column is pathways

% Extract Mean PFI values for NSD and HSD
mean_nsd = mean_data{:, 2};  
mean_hsd = mean_data{:, 3};  

% Extract SEM values for NSD and HSD
sem_nsd = sem_data{:, 2};  
sem_hsd = sem_data{:, 3};  

% Load count values
count_nsd = count_data{:, 2};
count_hsd = count_data{:, 3};
valid_count_idx = (count_nsd > 2) & (count_hsd > 2);

%% ------------------ Compute Log2 Fold Change & SEM Propagation ------------------
% Compute Log2 Fold Change (HSD/NSD)
log2_relative_change = log2(mean_hsd ./ mean_nsd);

% Compute SEM propagation for log2(HSD/NSD)
sem_ratio = sqrt((sem_hsd ./ mean_hsd).^2 + (sem_nsd ./ mean_nsd).^2);
sem_log2 = sem_ratio ./ log(2);  % Convert SEM to log2 scale

% Filter all arrays accordingly
pathway_names = pathway_names(valid_count_idx);
mean_nsd = mean_nsd(valid_count_idx);
mean_hsd = mean_hsd(valid_count_idx);
sem_nsd = sem_nsd(valid_count_idx);
sem_hsd = sem_hsd(valid_count_idx);
log2_relative_change = log2_relative_change(valid_count_idx);
sem_log2 = sem_log2(valid_count_idx);

%% ------------------ Extract Top 10 Increased & Decreased Pathways ------------------
num_top = 10; % Number of pathways to highlight

% Get indices for top 10 most increased (positive Log2 change)
[~, idx_increased] = maxk(log2_relative_change, num_top);
top_increased = log2_relative_change(idx_increased);
top_increased_pathways = pathway_names(idx_increased);
top_increased_sem = sem_log2(idx_increased); % Extract SEM

% Get indices for top 10 most decreased (negative Log2 change)
[~, idx_decreased] = mink(log2_relative_change, num_top);
top_decreased = log2_relative_change(idx_decreased);
top_decreased_pathways = pathway_names(idx_decreased);
top_decreased_sem = sem_log2(idx_decreased); % Extract SEM

% Define colors: Red for increased, Blue for decreased
colors_increased = repmat([0.8 0.2 0.2], num_top, 1); % Red
colors_decreased = repmat([0.2 0.2 0.8], num_top, 1); % Blue

%% ------------------ Generate Improved Top 10 Most Increased Pathways Plot with SEM ------------------
figure;
hold on;

% Horizontal bar plot
barh(flip(top_increased), 'FaceColor', 'flat'); 

% Apply colors
for i = 1:num_top
    barh(i, top_increased(num_top + 1 - i), 'FaceColor', colors_increased(i, :), 'EdgeColor', 'none');
end

% Adjust error bar positioning
errorbar(flip(top_increased), 1:num_top, flip(top_increased_sem), 'k.', 'horizontal', 'LineWidth', 1, 'CapSize', 6); 

% Adjust Y-axis labels and appearance
set(gca, 'YTickLabel', flip(top_increased_pathways), 'YTick', 1:num_top, 'FontSize', 10);
xlabel('Log2 Fold Change', 'FontSize', 12, 'FontWeight', 'bold');
title('Top 10 Most Increased rPFVI', 'FontSize', 14, 'FontWeight', 'bold');
set(gca,'FontName','Arial')
grid on;
xlim([min(log2_relative_change)-0.5, max(log2_relative_change)+0.5]); 
xlim([-0.8, 11]); 

% Save figure
saveas(gcf, fullfile(output_folder, 'PFI_Top10_Increased_v3_Fixed.png'));
saveas(gcf, fullfile(output_folder, 'PFI_Top10_Increased_v3_Fixed.svg'));

fprintf("✅ Top 10 most increased pathways plot saved correctly (with SEM)!\n");

%% ------------------ Generate Improved Top 10 Most Decreased Pathways Plot with SEM ------------------
figure;
hold on;

% Horizontal bar plot
barh(top_decreased, 'FaceColor', 'flat');

% Apply colors
for i = 1:num_top
    barh(i, top_decreased(i), 'FaceColor', colors_decreased(i, :), 'EdgeColor', 'none');
end

% Adjust error bar positioning
errorbar(top_decreased, 1:num_top, top_decreased_sem, 'k.', 'horizontal', 'LineWidth', 1, 'CapSize', 6); 



% Adjust Y-axis labels and appearance
set(gca, 'YTickLabel', top_decreased_pathways, 'YTick', 1:num_top, 'FontSize', 10);
xlabel('Log2 Fold Change', 'FontSize', 12, 'FontWeight', 'bold');
title('Top 10 Most Decreased rPFVI', 'FontSize', 14, 'FontWeight', 'bold');
set(gca,'FontName','Arial')

grid on;
xlim([-11, 0.8]); 

% Save figure
saveas(gcf, fullfile(output_folder, 'PFI_Top10_Decreased_v3_Fixed.png'));
saveas(gcf, fullfile(output_folder, 'PFI_Top10_Decreased_v3_Fixed.svg'));

fprintf("✅ Top 10 most decreased pathways plot saved correctly (with SEM)!\n");

%% ------------------ Save Updated Sorted Pathways to Excel ------------------
combined_results_table = table([top_decreased_pathways; top_increased_pathways], ...
                               [top_decreased; top_increased], ...
                               [top_decreased_sem; top_increased_sem], ...
    'VariableNames', {'Pathway', 'Log2FoldChange', 'SEM'});

% Save to Excel file
writetable(combined_results_table, fullfile(output_folder, 'PFI_Top10_Increased_Decreased_Corrected_v3_Fixed.xlsx'));

fprintf("✅ Updated results (with SEM) saved in 'PFI_Top10_Increased_Decreased_Corrected_v3_Fixed.xlsx'.\n");

%% ------------------ Generate Full Sorted Horizontal Pathway Change Plot (Log2) with SEM ------------------
%% ------------------ Remove NaN Values ------------------
valid_idx = ~isnan(log2_relative_change); % Find non-NaN indices
filtered_change = log2_relative_change(valid_idx);
filtered_pathways = pathway_names(valid_idx);
filtered_sem = sem_log2(valid_idx);

%% ------------------ Sort Pathways by Log2 Fold Change ------------------
% Sort so that the highest increased (red) is first, lowest decreased (blue) is last
[sorted_change, sort_idx] = sort(filtered_change, 'descend'); 
sorted_pathways = filtered_pathways(sort_idx);
sorted_sem = filtered_sem(sort_idx);

% ** FLIP ORDER ** So that RED (increased) is on top and BLUE (decreased) is at bottom
sorted_change = flipud(sorted_change);
sorted_pathways = flipud(sorted_pathways);
sorted_sem = flipud(sorted_sem);

%% ------------------ Generate Corrected Sorted Horizontal Pathway Change Plot (Log2) with SEM ------------------
figure;
hold on;

% Define colors: Red for increased (now at the top), Blue for decreased (now at the bottom)
bar_colors = zeros(length(sorted_change), 3);
bar_colors(sorted_change > 0, :) = repmat([0.8 0.2 0.2], sum(sorted_change > 0), 1); % Red (most increased)
bar_colors(sorted_change < 0, :) = repmat([0.2 0.2 0.8], sum(sorted_change < 0), 1); % Blue (most decreased)

% Create horizontal bar plot
barh(sorted_change, 'FaceColor', 'flat');

% Apply colors to bars
for i = 1:length(sorted_change)
    barh(i, sorted_change(i), 'FaceColor', bar_colors(i, :), 'EdgeColor', 'none');
end

% Add error bars (SEM)
errorbar(sorted_change, 1:length(sorted_change), sorted_sem, 'k.', 'horizontal', 'LineWidth', 0.5, 'CapSize', 3);

a = get(gca,'XTickLabel');
% set(gca,'XTickLabel',a,'FontName','Times','fontsize',18)

% Adjust Y-axis labels and appearance
set(gca, 'YTickLabel', sorted_pathways, 'YTick', 1:length(sorted_pathways), 'FontSize', 3.5);
% set(gca,'XTickLabel',a,'fontsize',10)
xlabel('Log2 Fold Change', 'FontSize', 10, 'FontWeight', 'bold');
title('Relative PFVI (HSD/NSD)', 'FontSize', 12, 'FontWeight', 'bold');
set(gca,'FontName','Arial')
grid on;

% Save the figure
saveas(gcf, fullfile(output_folder, 'PFI_Sorted_Pathways_Log2_Horizontal_Flipped_v3.png'));
saveas(gcf, fullfile(output_folder, 'PFI_Sorted_Pathways_Log2_Horizontal_Flipped_v3.pdf'));
saveas(gcf, fullfile(output_folder, 'PFI_Sorted_Pathways_Log2_Horizontal_Flipped_v3.svg'));

fprintf("✅ Fixed sorted pathway change plot with SEM (flipped order) saved!\n");

%% ------------------ Save Updated Sorted Pathways to Excel ------------------
sorted_change_table = table(sorted_pathways, sorted_change, sorted_sem, ...
    'VariableNames', {'Pathway', 'Log2FoldChange', 'SEM'});

writetable(sorted_change_table, fullfile(output_folder, 'PFI_Sorted_Pathways_Log2_Flipped_v3.csv'));

fprintf("✅ Updated sorted pathways (flipped order) saved in 'PFI_Sorted_Pathways_Log2_Flipped_v3.csv'.\n");
