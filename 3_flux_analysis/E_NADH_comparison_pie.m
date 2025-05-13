%% Goal:
%   - Extract raw simulation data for corresponding reactions
%   - Generate boxplots comparing CTR vs. HSD for each reaction
%   - Compute Wilcoxon rank-sum test for CTR vs. HSD
%   - Compare total NAD/NADH flux with bar chart
%   - Create a pie chart comparing mean fluxes
%   - Generate boxplots comparing NAD-regeneration vs. NADH-generation

close all
%% Define Paths
pathway = pwd;
save_dirFlux = 'E';
subfolder = fullfile(pathway, save_dirFlux);
if ~exist(subfolder, 'dir'), mkdir(subfolder), end

load_folder = fullfile(pathway, 'D_flux_analysis', 'log', 'nadC');
excel_file = fullfile(load_folder, 'flux_data_raw.xlsx');

% Load data
CTR_all = readtable(excel_file, 'Sheet', 'NSD'); % CTR instead of NSD
HSD_all = readtable(excel_file, 'Sheet', 'HSD');

% Define reaction IDs and their labels
rxn_id = {'MAR04683', 'MAR08507', 'MAR04388', 'MAR00479'}; % Aldh3, T3dh, Ldh, Gpdh1
rxn_labels = {'Aldh3', 'T3dh', 'Ldh', 'Gpdh1'};

% Extract data for the selected reactions
CTR_data = CTR_all(:, ismember(CTR_all.Properties.VariableNames, rxn_id));
HSD_data = HSD_all(:, ismember(HSD_all.Properties.VariableNames, rxn_id));

% Compute mean values for pie chart
CTR_means = mean(table2array(CTR_data), 'omitnan');
HSD_means = mean(table2array(HSD_data), 'omitnan');

% Compute total NAD/NADH flux contribution
total_CTR = sum(CTR_means);
total_HSD = sum(HSD_means);

% Perform statistical tests
p_values = zeros(length(rxn_id), 1);
for i = 1:length(rxn_id)
    [p_values(i), ~] = ranksum(CTR_data{:, i}, HSD_data{:, i}); % Wilcoxon rank-sum test
end

%% Box Plots for Each Reaction (CTR vs HSD)
figure;
t = tiledlayout(2,2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

for i = 1:length(rxn_id)
    nexttile
    
    % Format p-value display
    if p_values(i) < 1e-6
        p_text = 'p < 1e-6';
    else
        p_text = sprintf('p = %.4f', p_values(i));
    end
    
    % Create categorical labels to ensure correct placement
    group_labels = categorical([repmat("CTR", size(CTR_data{:, i})); repmat("HSD", size(HSD_data{:, i}))]);
    
    boxchart(group_labels, [CTR_data{:, i}; HSD_data{:, i}])
    title([rxn_labels{i}, ' (', p_text, ')']) % Add formatted p-value to title
    ylabel('Flux')
    grid on
end

% Global labels
xlabel(t, 'Condition')
ylabel(t, 'Flux')
sgtitle('Flux Comparisons for NAD/NADH-Related Reactions')

% Save figures
saveas(gcf, fullfile(subfolder, 'boxplots_NAD_NADH.png'))
saveas(gcf, fullfile(subfolder, 'boxplots_NAD_NADH.svg'))


%% Combined Data for NAD-Regeneration and NADH-Generation

% Combine individual data points rather than summing
CTR_NAD_regeneration = [CTR_data{:,'MAR04388'}; CTR_data{:,'MAR00479'}]; % Ldh + Gpdh1
HSD_NAD_regeneration = [HSD_data{:,'MAR04388'}; HSD_data{:,'MAR00479'}];

CTR_NADH_generation = [CTR_data{:,'MAR04683'}; CTR_data{:,'MAR08507'}]; % Aldh3 + T3dh
HSD_NADH_generation = [HSD_data{:,'MAR04683'}; HSD_data{:,'MAR08507'}];

% Perform statistical tests
[p_NAD_reg, ~] = ranksum(CTR_NAD_regeneration, HSD_NAD_regeneration);
[p_NADH_gen, ~] = ranksum(CTR_NADH_generation, HSD_NADH_generation);
% Create a grouped boxplot comparing NAD-regeneration and NADH-generation

figure;
hold on

% Define categories for box grouping
group_labels = [repmat("NAD Regeneration", size(CTR_NAD_regeneration));
                repmat("NAD Regeneration", size(HSD_NAD_regeneration));
                repmat("NADH Generation", size(CTR_NADH_generation));
                repmat("NADH Generation", size(HSD_NADH_generation))];

condition_labels = [repmat("CTR", size(CTR_NAD_regeneration));
                    repmat("HSD", size(HSD_NAD_regeneration));
                    repmat("CTR", size(CTR_NADH_generation));
                    repmat("HSD", size(HSD_NADH_generation))];

data_values = [CTR_NAD_regeneration; HSD_NAD_regeneration; ...
               CTR_NADH_generation; HSD_NADH_generation];

% Convert categorical labels
group_labels = categorical(group_labels);
condition_labels = categorical(condition_labels);

% Create boxchart with groups
h = boxchart(group_labels, data_values, 'GroupByColor', condition_labels);

% Assign colors: black for CTR, red for HSD
hold on
for i = 1:length(h)
    if h(i).SeriesIndex == 1
        h(i).BoxFaceColor = [0 0 0]; % Black for CTR
        h(i).MarkerColor = [0 0 0]; % Keep outlier bubbles black
    elseif h(i).SeriesIndex == 2
        h(i).BoxFaceColor = [1 0 0]; % Red for HSD
        h(i).MarkerColor = [1 0 0]; % Keep outlier bubbles red
    end
end
hold off

% Formatting
ylabel('Flux')
legend({'CTR', 'HSD'}, 'Location', 'eastoutside')
set(gca, 'FontSize', 13)
legend boxoff
grid on

hold off

% Save figure
saveas(gcf, fullfile(subfolder, 'boxplots_combined_NAD_NADH.png'))
saveas(gcf, fullfile(subfolder, 'boxplots_combined_NAD_NADH.svg'))

%% Box Chart for Total NAD/NADH Flux
% Perform statistical test (Wilcoxon rank-sum test) for total NAD/NADH flux
[p_total, ~] = ranksum(sum(table2array(CTR_data), 2), sum(table2array(HSD_data), 2));

% Format p-value display
if p_total < 1e-6
    p_text = 'p < 1e-6';
else
    p_text = sprintf('p = %.4f', p_total);
end

% Prepare data for boxchart
total_CTR_values = sum(table2array(CTR_data), 2);
total_HSD_values = sum(table2array(HSD_data), 2);
group_labels = [repmat("CTR", size(total_CTR_values)); repmat("HSD", size(total_HSD_values))];
flux_values = [total_CTR_values; total_HSD_values];

% Create boxchart
figure;
boxchart(categorical(group_labels), flux_values)
title(['Total NAD/NADH Flux Comparison (', p_text, ')'])
ylabel('Total NAD/NADH Flux')
grid on

% Save figure
saveas(gcf, fullfile(subfolder, 'boxchart_total_NAD_NADH.png'))
saveas(gcf, fullfile(subfolder, 'boxchart_total_NAD_NADH.svg'))

%% Pie Chart Comparing Mean Fluxes
rxn_labels_new = {'GPDH1', 'Ldh', 'Aldh3', 'T3dh'};

figure;
t = tiledlayout(1,2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% CTR Pie Chart
nexttile
pie(CTR_means)
title('CTR (Mean Flux)')

% HSD Pie Chart
nexttile
pie(HSD_means)
title('HSD (Mean Flux)')

% Create a single legend
lgd = legend(rxn_labels_new, 'Location', 'southoutside', 'Orientation', 'horizontal');
legend boxoff
lgd.Layout.Tile = 'south'; % Position the legend below both pies

% Save figures
saveas(gcf, fullfile(subfolder, 'piechart_means.png'))
saveas(gcf, fullfile(subfolder, 'piechart_means.svg'))


%% Save Mean Values and Statistics to Excel
results_table = table(rxn_labels', CTR_means', HSD_means', p_values, ...
    'VariableNames', {'Reaction', 'CTR_Mean', 'HSD_Mean', 'p_value'});

% Summed Flux Results
summed_flux_table = table(["NAD Regeneration"; "NADH Generation"], ...
    [mean(CTR_NAD_regeneration); mean(CTR_NADH_generation)], ...
    [mean(HSD_NAD_regeneration); mean(HSD_NADH_generation)], ...
    [p_NAD_reg; p_NADH_gen], ...
    'VariableNames', {'Category', 'CTR_Summed_Mean', 'HSD_Summed_Mean', 'p_value'});

% Save to Excel
writetable(results_table, fullfile(subfolder, 'NAD_NADH_flux_means.xlsx'), 'Sheet', 'Means');
writetable(summed_flux_table, fullfile(subfolder, 'NAD_NADH_flux_means.xlsx'), 'Sheet', 'Summed_Flux');

disp('Figures and statistical results saved successfully.');
