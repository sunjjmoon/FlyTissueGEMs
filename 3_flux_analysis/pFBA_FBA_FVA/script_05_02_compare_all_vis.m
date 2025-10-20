%% Visualization Script: FBA vs pFBA vs sampling Comparison

clc;
clear;
close all;

%% Load data
current_path = pwd;
results_dir = fullfile(current_path, '04_sampling_comp');

load(fullfile(results_dir, 'fba_pfba_sampling_comparison.mat'));
load(fullfile(current_path,'01_results','model_constrained_out.mat'));

%% Setup visualization directory
vis_dir = fullfile(current_path, '05_visualizations');
if ~exist(vis_dir, 'dir')
    mkdir(vis_dir);
end


%% Separate Panel Figures: NSD (150-152) and HSD (160-162)

% Common settings
methods = {'FBA', 'pFBA', 'FVA-sampling'};
threshold = 0;

%% Figure 180: Zero-flux Reactions (Both Models)
fig180 = figure(180);
clf;

% Prepare data
nsd_data = comparison_data.NSD;
hsd_data = comparison_data.HSD;

% NSD
nsd_fba = abs(nsd_data.fba_fluxes);
nsd_pfba = abs(nsd_data.pfba_fluxes);
nsd_samp = abs(nsd_data.sampling_median);
nsd_total = length(nsd_fba);
nsd_zeros = [sum(nsd_fba <= threshold)/nsd_total, sum(nsd_pfba <= threshold)/nsd_total, sum(nsd_samp <= threshold)/nsd_total] * 100;

% HSD
hsd_fba = abs(hsd_data.fba_fluxes);
hsd_pfba = abs(hsd_data.pfba_fluxes);
hsd_samp = abs(hsd_data.sampling_median);
hsd_total = length(hsd_fba);
hsd_zeros = [sum(hsd_fba <= threshold)/hsd_total, sum(hsd_pfba <= threshold)/hsd_total, sum(hsd_samp <= threshold)/hsd_total] * 100;

% Create grouped bar
data_matrix = [nsd_zeros; hsd_zeros];
b = bar(data_matrix);
set(gca, 'XTickLabel', {'NSD', 'HSD'});
ylabel('% zero-flux reactions', 'FontSize', 14, 'FontWeight', 'bold');
% title('Inactive Reactions', 'FontSize', 16, 'FontWeight', 'bold');
legend({'FBA', 'pFBA', 'FVA-sampling'}, 'Location', 'best');
grid on;
ylim([0 100]);

set(gcf, 'Position', [100, 100, 600, 450]);

savefig(fig180, fullfile(vis_dir, '180_Grouped_Inactive_Reactions.fig'));
saveas(fig180, fullfile(vis_dir, '180_Grouped_Inactive_Reactions.png'), 'png');
set(fig180, 'Units', 'inches');
screenposition = get(fig180, 'Position');
set(fig180, 'PaperPosition', [0 0 screenposition(3:4)], ...
         'PaperSize', [screenposition(3:4)], ...
         'PaperUnits', 'inches', ...
         'PaperOrientation', 'portrait');
print(fig180, fullfile(vis_dir, '180_Grouped_Inactive_Reactions.pdf'), '-dpdf', '-painters');

%% Figure 181: Active Reaction Distribution (Both Models)
fig181 = figure(181);
clf;

% NSD
nsd_fba_nonzero = nsd_fba(nsd_fba >= threshold);
nsd_pfba_nonzero = nsd_pfba(nsd_pfba >= threshold);
nsd_samp_nonzero = nsd_samp(nsd_samp >= threshold);

subplot(1,2,1);
[fba_f, fba_x] = ecdf(nsd_fba_nonzero);
[pfba_f, pfba_x] = ecdf(nsd_pfba_nonzero);
[samp_f, samp_x] = ecdf(nsd_samp_nonzero);

semilogx(fba_x, fba_f, 'LineWidth', 2.5, 'DisplayName', 'FBA');
hold on;
semilogx(pfba_x, pfba_f, 'LineWidth', 2.5, 'DisplayName', 'pFBA');
semilogx(samp_x, samp_f, 'LineWidth', 2.5, 'DisplayName', 'FVA-sampling');

xlabel('Flux (a.u.)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel(['Cumulative probability' newline 'of non-zero reactions'], 'FontSize', 14, 'FontWeight', 'bold');
title('NSD', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'southeast', 'FontSize', 12);
xlim([1e-3 1e5]);
grid on;

% HSD
hsd_fba_nonzero = hsd_fba(hsd_fba >= threshold);
hsd_pfba_nonzero = hsd_pfba(hsd_pfba >= threshold);
hsd_samp_nonzero = hsd_samp(hsd_samp >= threshold);

subplot(1,2,2);
[fba_f, fba_x] = ecdf(hsd_fba_nonzero);
[pfba_f, pfba_x] = ecdf(hsd_pfba_nonzero);
[samp_f, samp_x] = ecdf(hsd_samp_nonzero);

semilogx(fba_x, fba_f, 'LineWidth', 2.5, 'DisplayName', 'FBA');
hold on;
semilogx(pfba_x, pfba_f, 'LineWidth', 2.5, 'DisplayName', 'pFBA');
semilogx(samp_x, samp_f, 'LineWidth', 2.5, 'DisplayName', 'FVA-sampling');

xlabel('Flux (a.u.)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel(['Cumulative probability' newline 'of non-zero reactions'], 'FontSize', 14, 'FontWeight', 'bold');
title('HSD', 'FontSize', 16, 'FontWeight', 'bold');
legend('Location', 'southeast', 'FontSize', 12);
xlim([1e-3 1e5]);
grid on;

set(gcf, 'Position', [100, 100, 1000, 450]);

savefig(fig181, fullfile(vis_dir, '181_Grouped_Active_Distribution.fig'));
saveas(fig181, fullfile(vis_dir, '181_Grouped_Active_Distribution.png'), 'png');
set(fig181, 'Units', 'inches');
screenposition = get(fig181, 'Position');
set(fig181, 'PaperPosition', [0 0 screenposition(3:4)], ...
         'PaperSize', [screenposition(3:4)], ...
         'PaperUnits', 'inches', ...
         'PaperOrientation', 'portrait');
print(fig181, fullfile(vis_dir, '181_Grouped_Active_Distribution.pdf'), '-dpdf', '-painters');

%% Figure 182: Total Flux (Both Models)
fig182 = figure(182);
clf;

% Calculate total fluxes
nsd_totals = [sum(nsd_fba), sum(nsd_pfba), sum(nsd_samp)];
hsd_totals = [sum(hsd_fba), sum(hsd_pfba), sum(hsd_samp)];

% Create grouped bar
data_matrix = [nsd_totals; hsd_totals];
b = bar(data_matrix);
set(gca, 'XTickLabel', {'NSD', 'HSD'},'fontsize',16);
ylabel('Total flux (a.u.)', 'FontSize', 18, 'FontWeight', 'bold');
% title('Total Metabolic Activity', 'FontSize', 16, 'FontWeight', 'bold');
legend({'FBA', 'pFBA', 'FVA-sampling'}, 'Location', 'northeast','fontsize',12);
grid on;

set(gcf, 'Position', [100, 100, 600, 450]);

savefig(fig182, fullfile(vis_dir, '182_Grouped_Total_Flux.fig'));
saveas(fig182, fullfile(vis_dir, '182_Grouped_Total_Flux.png'), 'png');
set(fig182, 'Units', 'inches');
screenposition = get(fig182, 'Position');
set(fig182, 'PaperPosition', [0 0 screenposition(3:4)], ...
         'PaperSize', [screenposition(3:4)], ...
         'PaperUnits', 'inches', ...
         'PaperOrientation', 'portrait');
print(fig182, fullfile(vis_dir, '182_Grouped_Total_Flux.pdf'), '-dpdf', '-painters');

%% Figure 190: Correlation Analysis (Both Models in Subplots)
fig190 = figure(190);
clf;

% NSD Correlation
model_data = comparison_data.NSD;
fba_fluxes = abs(model_data.fba_fluxes);
pfba_fluxes = abs(model_data.pfba_fluxes);
sampling_fluxes = abs(model_data.sampling_median);

corr_fba_pfba = corr(fba_fluxes, pfba_fluxes);
corr_fba_samp = corr(fba_fluxes, sampling_fluxes);
corr_pfba_samp = corr(pfba_fluxes, sampling_fluxes);

corr_matrix_nsd = [1, corr_fba_pfba, corr_fba_samp;
                   corr_fba_pfba, 1, corr_pfba_samp;
                   corr_fba_samp, corr_pfba_samp, 1];

subplot(1,2,1);
imagesc(corr_matrix_nsd);
colorbar;
caxis([0 1]);
colormap('gray');
set(gca, 'XTick', 1:3, 'XTickLabel', {'FBA', 'pFBA', 'FVA-sampling'});
set(gca, 'YTick', 1:3, 'YTickLabel', {'FBA', 'pFBA', 'FVA-sampling'});
title('NSD', 'FontSize', 14, 'FontWeight', 'bold');

for i = 1:3
    for j = 1:3
        text(j, i, sprintf('%.3f', corr_matrix_nsd(i,j)), ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 12, 'FontWeight', 'bold', 'Color', 'white');
    end
end

% HSD Correlation
model_data = comparison_data.HSD;
fba_fluxes = abs(model_data.fba_fluxes);
pfba_fluxes = abs(model_data.pfba_fluxes);
sampling_fluxes = abs(model_data.sampling_median);

corr_fba_pfba = corr(fba_fluxes, pfba_fluxes);
corr_fba_samp = corr(fba_fluxes, sampling_fluxes);
corr_pfba_samp = corr(pfba_fluxes, sampling_fluxes);

corr_matrix_hsd = [1, corr_fba_pfba, corr_fba_samp;
                   corr_fba_pfba, 1, corr_pfba_samp;
                   corr_fba_samp, corr_pfba_samp, 1];

subplot(1,2,2);
imagesc(corr_matrix_hsd);
colorbar;
caxis([0 1]);
colormap('gray');
set(gca, 'XTick', 1:3, 'XTickLabel', {'FBA', 'pFBA', 'FVA-sampling'});
set(gca, 'YTick', 1:3, 'YTickLabel', {'FBA', 'pFBA', 'FVA-sampling'});
title('HSD', 'FontSize', 14, 'FontWeight', 'bold');

for i = 1:3
    for j = 1:3
        text(j, i, sprintf('%.3f', corr_matrix_hsd(i,j)), ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 12, 'FontWeight', 'bold', 'Color', 'white');
    end
end

sgtitle('Flux Correlation', 'FontSize', 16, 'FontWeight', 'bold');

set(gcf, 'Position', [100, 100, 1000, 450]);

savefig(fig190, fullfile(vis_dir, '190_Grouped_Correlation.fig'));
saveas(fig190, fullfile(vis_dir, '190_Grouped_Correlation.png'), 'png');
set(fig190, 'Units', 'inches');
screenposition = get(fig190, 'Position');
set(fig190, 'PaperPosition', [0 0 screenposition(3:4)], ...
         'PaperSize', [screenposition(3:4)], ...
         'PaperUnits', 'inches', ...
         'PaperOrientation', 'portrait');
print(fig190, fullfile(vis_dir, '190_Grouped_Correlation.pdf'), '-dpdf', '-painters');

%% Export results for grouped figures
excel_filename = fullfile(vis_dir, 'Grouped_Figures_Results.xlsx');
fprintf('\nExporting grouped figure results to: %s\n', excel_filename);

% Correlation data (same as before)
sheet_name = 'Correlation_Matrix';
corr_data = {
    'Model', 'FBA_vs_pFBA', 'FBA_vs_sampling', 'pFBA_vs_sampling';
    'NSD', corr_matrix_nsd(1,2), corr_matrix_nsd(1,3), corr_matrix_nsd(2,3);
    'HSD', corr_matrix_hsd(1,2), corr_matrix_hsd(1,3), corr_matrix_hsd(2,3)
};
writecell(corr_data, excel_filename, 'Sheet', sheet_name);

fprintf('\n✓ Grouped figures complete!\n');
fprintf('Figures: 180 (Inactive), 181 (Distribution), 182 (Total Flux), 190 (Correlation)\n');
fprintf('Results exported to: %s\n', excel_filename);






%% Export Results for Figures 150-162 to Excel
% Add this code after generating figures 150-162

excel_filename = fullfile(vis_dir, 'Figures_Results.xlsx');
fprintf('\nExporting figure results to: %s\n', excel_filename);

%% Prepare common data
methods = {'FBA', 'pFBA', 'sampling'};
models = {'NSD', 'HSD'};
threshold = 0;

%% Sheet 1: Zero-flux Reactions Summary (Figures 150, 160)
sheet_name = '01_Zero_Flux_Summary';
fprintf('Writing sheet: %s\n', sheet_name);

zero_data = {};
row_idx = 1;
zero_data{row_idx, 1} = 'Model';
zero_data{row_idx, 2} = 'Method';
zero_data{row_idx, 3} = 'Zero_Flux_Count';
zero_data{row_idx, 4} = 'Total_Reactions';
zero_data{row_idx, 5} = 'Zero_Flux_Percent';
row_idx = row_idx + 1;

for model_idx = 1:2
    if model_idx == 1
        model_data = comparison_data.NSD;
        model_name = 'NSD';
    else
        model_data = comparison_data.HSD;
        model_name = 'HSD';
    end
    
    fba_fluxes = abs(model_data.fba_fluxes);
    pfba_fluxes = abs(model_data.pfba_fluxes);
    sampling_fluxes = abs(model_data.sampling_median);
    
    total_reactions = length(fba_fluxes);
    fba_zeros = sum(fba_fluxes <= threshold);
    pfba_zeros = sum(pfba_fluxes <= threshold);
    sampling_zeros = sum(sampling_fluxes <= threshold);
    
    % FBA
    zero_data{row_idx, 1} = model_name;
    zero_data{row_idx, 2} = 'FBA';
    zero_data{row_idx, 3} = fba_zeros;
    zero_data{row_idx, 4} = total_reactions;
    zero_data{row_idx, 5} = (fba_zeros/total_reactions) * 100;
    row_idx = row_idx + 1;
    
    % pFBA
    zero_data{row_idx, 1} = model_name;
    zero_data{row_idx, 2} = 'pFBA';
    zero_data{row_idx, 3} = pfba_zeros;
    zero_data{row_idx, 4} = total_reactions;
    zero_data{row_idx, 5} = (pfba_zeros/total_reactions) * 100;
    row_idx = row_idx + 1;
    
    % sampling
    zero_data{row_idx, 1} = model_name;
    zero_data{row_idx, 2} = 'sampling';
    zero_data{row_idx, 3} = sampling_zeros;
    zero_data{row_idx, 4} = total_reactions;
    zero_data{row_idx, 5} = (sampling_zeros/total_reactions) * 100;
    row_idx = row_idx + 1;
end

zero_table = cell2table(zero_data(2:end, :), 'VariableNames', zero_data(1, :));
writetable(zero_table, excel_filename, 'Sheet', sheet_name);

%% Sheet 2: Active Reaction Distribution Statistics (Figures 151, 161)
sheet_name = '02_Active_Distribution_Stats';
fprintf('Writing sheet: %s\n', sheet_name);

dist_data = {};
row_idx = 1;
dist_data{row_idx, 1} = 'Model';
dist_data{row_idx, 2} = 'Method';
dist_data{row_idx, 3} = 'N_Active';
dist_data{row_idx, 4} = 'Min';
dist_data{row_idx, 5} = 'Q25';
dist_data{row_idx, 6} = 'Median';
dist_data{row_idx, 7} = 'Mean';
dist_data{row_idx, 8} = 'Q75';
dist_data{row_idx, 9} = 'Max';
dist_data{row_idx, 10} = 'Std';
row_idx = row_idx + 1;

for model_idx = 1:2
    if model_idx == 1
        model_data = comparison_data.NSD;
        model_name = 'NSD';
    else
        model_data = comparison_data.HSD;
        model_name = 'HSD';
    end
    
    fba_fluxes = abs(model_data.fba_fluxes);
    pfba_fluxes = abs(model_data.pfba_fluxes);
    sampling_fluxes = abs(model_data.sampling_median);
    
    fba_nonzero = fba_fluxes(fba_fluxes > threshold);
    pfba_nonzero = pfba_fluxes(pfba_fluxes > threshold);
    sampling_nonzero = sampling_fluxes(sampling_fluxes > threshold);
    
    % FBA
    dist_data{row_idx, 1} = model_name;
    dist_data{row_idx, 2} = 'FBA';
    dist_data{row_idx, 3} = length(fba_nonzero);
    dist_data{row_idx, 4} = min(fba_nonzero);
    dist_data{row_idx, 5} = quantile(fba_nonzero, 0.25);
    dist_data{row_idx, 6} = median(fba_nonzero);
    dist_data{row_idx, 7} = mean(fba_nonzero);
    dist_data{row_idx, 8} = quantile(fba_nonzero, 0.75);
    dist_data{row_idx, 9} = max(fba_nonzero);
    dist_data{row_idx, 10} = std(fba_nonzero);
    row_idx = row_idx + 1;
    
    % pFBA
    dist_data{row_idx, 1} = model_name;
    dist_data{row_idx, 2} = 'pFBA';
    dist_data{row_idx, 3} = length(pfba_nonzero);
    dist_data{row_idx, 4} = min(pfba_nonzero);
    dist_data{row_idx, 5} = quantile(pfba_nonzero, 0.25);
    dist_data{row_idx, 6} = median(pfba_nonzero);
    dist_data{row_idx, 7} = mean(pfba_nonzero);
    dist_data{row_idx, 8} = quantile(pfba_nonzero, 0.75);
    dist_data{row_idx, 9} = max(pfba_nonzero);
    dist_data{row_idx, 10} = std(pfba_nonzero);
    row_idx = row_idx + 1;
    
    % sampling
    dist_data{row_idx, 1} = model_name;
    dist_data{row_idx, 2} = 'sampling';
    dist_data{row_idx, 3} = length(sampling_nonzero);
    dist_data{row_idx, 4} = min(sampling_nonzero);
    dist_data{row_idx, 5} = quantile(sampling_nonzero, 0.25);
    dist_data{row_idx, 6} = median(sampling_nonzero);
    dist_data{row_idx, 7} = mean(sampling_nonzero);
    dist_data{row_idx, 8} = quantile(sampling_nonzero, 0.75);
    dist_data{row_idx, 9} = max(sampling_nonzero);
    dist_data{row_idx, 10} = std(sampling_nonzero);
    row_idx = row_idx + 1;
end

dist_table = cell2table(dist_data(2:end, :), 'VariableNames', dist_data(1, :));
writetable(dist_table, excel_filename, 'Sheet', sheet_name);

%% Sheet 3: Total Flux Summary (Figures 152, 162)
sheet_name = '03_Total_Flux_Summary';
fprintf('Writing sheet: %s\n', sheet_name);

total_data = {};
row_idx = 1;
total_data{row_idx, 1} = 'Model';
total_data{row_idx, 2} = 'Method';
total_data{row_idx, 3} = 'Total_Flux';
total_data{row_idx, 4} = 'Average_Flux';
total_data{row_idx, 5} = 'N_Active_Reactions';
row_idx = row_idx + 1;

for model_idx = 1:2
    if model_idx == 1
        model_data = comparison_data.NSD;
        model_name = 'NSD';
    else
        model_data = comparison_data.HSD;
        model_name = 'HSD';
    end
    
    fba_fluxes = abs(model_data.fba_fluxes);
    pfba_fluxes = abs(model_data.pfba_fluxes);
    sampling_fluxes = abs(model_data.sampling_median);
    
    fba_nonzero = fba_fluxes(fba_fluxes > threshold);
    pfba_nonzero = pfba_fluxes(pfba_fluxes > threshold);
    sampling_nonzero = sampling_fluxes(sampling_fluxes > threshold);
    
    % FBA
    total_data{row_idx, 1} = model_name;
    total_data{row_idx, 2} = 'FBA';
    total_data{row_idx, 3} = sum(fba_fluxes);
    total_data{row_idx, 4} = mean(fba_nonzero);
    total_data{row_idx, 5} = length(fba_nonzero);
    row_idx = row_idx + 1;
    
    % pFBA
    total_data{row_idx, 1} = model_name;
    total_data{row_idx, 2} = 'pFBA';
    total_data{row_idx, 3} = sum(pfba_fluxes);
    total_data{row_idx, 4} = mean(pfba_nonzero);
    total_data{row_idx, 5} = length(pfba_nonzero);
    row_idx = row_idx + 1;
    
    % sampling
    total_data{row_idx, 1} = model_name;
    total_data{row_idx, 2} = 'sampling';
    total_data{row_idx, 3} = sum(sampling_fluxes);
    total_data{row_idx, 4} = mean(sampling_nonzero);
    total_data{row_idx, 5} = length(sampling_nonzero);
    row_idx = row_idx + 1;
end

total_table = cell2table(total_data(2:end, :), 'VariableNames', total_data(1, :));
writetable(total_table, excel_filename, 'Sheet', sheet_name);

%% Sheet 4: Flux Reduction Analysis
sheet_name = '04_Flux_Reduction';
fprintf('Writing sheet: %s\n', sheet_name);

reduction_data = {};
row_idx = 1;
reduction_data{row_idx, 1} = 'Model';
reduction_data{row_idx, 2} = 'Comparison';
reduction_data{row_idx, 3} = 'Total_Flux_Reduction';
reduction_data{row_idx, 4} = 'Percent_Reduction';
reduction_data{row_idx, 5} = 'Active_Reactions_Change';
row_idx = row_idx + 1;

for model_idx = 1:2
    if model_idx == 1
        model_data = comparison_data.NSD;
        model_name = 'NSD';
    else
        model_data = comparison_data.HSD;
        model_name = 'HSD';
    end
    
    fba_fluxes = abs(model_data.fba_fluxes);
    pfba_fluxes = abs(model_data.pfba_fluxes);
    sampling_fluxes = abs(model_data.sampling_median);
    
    fba_total = sum(fba_fluxes);
    pfba_total = sum(pfba_fluxes);
    sampling_total = sum(sampling_fluxes);
    
    fba_active = sum(fba_fluxes > threshold);
    pfba_active = sum(pfba_fluxes > threshold);
    sampling_active = sum(sampling_fluxes > threshold);
    
    % FBA vs pFBA
    reduction_data{row_idx, 1} = model_name;
    reduction_data{row_idx, 2} = 'FBA_vs_pFBA';
    reduction_data{row_idx, 3} = fba_total - pfba_total;
    reduction_data{row_idx, 4} = ((fba_total - pfba_total) / fba_total) * 100;
    reduction_data{row_idx, 5} = fba_active - pfba_active;
    row_idx = row_idx + 1;
    
    % FBA vs sampling
    reduction_data{row_idx, 1} = model_name;
    reduction_data{row_idx, 2} = 'FBA_vs_sampling';
    reduction_data{row_idx, 3} = fba_total - sampling_total;
    reduction_data{row_idx, 4} = ((fba_total - sampling_total) / fba_total) * 100;
    reduction_data{row_idx, 5} = fba_active - sampling_active;
    row_idx = row_idx + 1;
    
    % pFBA vs sampling
    reduction_data{row_idx, 1} = model_name;
    reduction_data{row_idx, 2} = 'pFBA_vs_sampling';
    reduction_data{row_idx, 3} = pfba_total - sampling_total;
    reduction_data{row_idx, 4} = ((pfba_total - sampling_total) / pfba_total) * 100;
    reduction_data{row_idx, 5} = pfba_active - sampling_active;
    row_idx = row_idx + 1;
end

reduction_table = cell2table(reduction_data(2:end, :), 'VariableNames', reduction_data(1, :));
writetable(reduction_table, excel_filename, 'Sheet', sheet_name);

%% Sheet 5: Summary for Manuscript
sheet_name = '05_Manuscript_Summary';
fprintf('Writing sheet: %s\n', sheet_name);

summary_data = {
    'Metric', 'NSD', 'HSD', 'Unit/Description';
    '', '', '', '';
    '=== INACTIVE REACTIONS (Figures 150, 160) ===', '', '', '';
    'Zero-flux % (FBA)', '', '', '%';
    'Zero-flux % (pFBA)', '', '', '%';
    'Zero-flux % (sampling)', '', '', '%';
    '', '', '', '';
    '=== TOTAL FLUX (Figures 152, 162) ===', '', '', '';
    'Total Flux (FBA)', '', '', 'a.u.';
    'Total Flux (pFBA)', '', '', 'a.u.';
    'Total Flux (sampling)', '', '', 'a.u.';
    '', '', '', '';
    '=== FLUX REDUCTION ===', '', '', '';
    'FBA to pFBA reduction', '', '', '%';
    'pFBA to sampling increase', '', '', '%';
    '', '', '', '';
    '=== ACTIVE REACTIONS ===', '', '', '';
    'Active reactions (FBA)', '', '', 'count';
    'Active reactions (pFBA)', '', '', 'count';
    'Active reactions (sampling)', '', '', 'count'
};

% Fill in NSD data
for model_idx = 1:2
    if model_idx == 1
        model_data = comparison_data.NSD;
        col = 2;
    else
        model_data = comparison_data.HSD;
        col = 3;
    end
    
    fba_fluxes = abs(model_data.fba_fluxes);
    pfba_fluxes = abs(model_data.pfba_fluxes);
    sampling_fluxes = abs(model_data.sampling_median);
    
    total_reactions = length(fba_fluxes);
    fba_zeros = sum(fba_fluxes <= threshold);
    pfba_zeros = sum(pfba_fluxes <= threshold);
    sampling_zeros = sum(sampling_fluxes <= threshold);
    
    summary_data{4, col} = (fba_zeros/total_reactions) * 100;
    summary_data{5, col} = (pfba_zeros/total_reactions) * 100;
    summary_data{6, col} = (sampling_zeros/total_reactions) * 100;
    
    summary_data{9, col} = sum(fba_fluxes);
    summary_data{10, col} = sum(pfba_fluxes);
    summary_data{11, col} = sum(sampling_fluxes);
    
    summary_data{14, col} = ((sum(fba_fluxes) - sum(pfba_fluxes)) / sum(fba_fluxes)) * 100;
    summary_data{15, col} = ((sum(sampling_fluxes) - sum(pfba_fluxes)) / sum(pfba_fluxes)) * 100;
    
    summary_data{18, col} = sum(fba_fluxes > threshold);
    summary_data{19, col} = sum(pfba_fluxes > threshold);
    summary_data{20, col} = sum(sampling_fluxes > threshold);
end

writecell(summary_data, excel_filename, 'Sheet', sheet_name);

fprintf('\n✓ Export complete! File saved as: %s\n', excel_filename);
fprintf('Total sheets created: 5\n\n');

