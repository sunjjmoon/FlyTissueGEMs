%% script_02_NADHdemand_pFBA_vis: Visualize NADH maximum Capacity Results

clc;
clear;
close all;

%% 1. Load Results from Script 01
results_dir = fullfile(pwd, '01_NADH_demand_pFBA');
load(fullfile(results_dir, 'enhanced_FBA_results.mat'));
load(fullfile(results_dir, 'comparative_results.mat'));

output_dir = fullfile(pwd, '02_NADH_demand_pFBA_vis');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% 2. Extract Data for Visualization
n_models = length(enhanced_results);
model_names = cell(n_models, 1);
nadh_oxidation_capacity = zeros(n_models, 1);
biomass_rates = zeros(n_models, 1);
nadh_to_biomass_ratio = zeros(n_models, 1);

for i = 1:n_models
    model_names{i} = enhanced_results{i}.modelID;
    nadh_oxidation_capacity(i) = enhanced_results{i}.max_nadh_oxidation_capacity;
    biomass_rates(i) = enhanced_results{i}.biomass_growth_rate;
    
    % Calculate NADH oxidation capacity per unit biomass
    if biomass_rates(i) > 0
        nadh_to_biomass_ratio(i) = nadh_oxidation_capacity(i) / biomass_rates(i);
    else
        nadh_to_biomass_ratio(i) = 0;
    end
end

%% 3. Figure 1: NADH Oxidation Capacity Comparison
% figure('Position', [100, 100, 800, 600]);
nadh_oxidation_capacity_rel = nadh_oxidation_capacity(2)/nadh_oxidation_capacity(1);
figure(1)
bar_handle = bar(nadh_oxidation_capacity, 'EdgeColor', 'black', 'LineWidth', 1);
set(gca, 'XTickLabel', model_names, 'FontSize', 14, 'FontWeight', 'bold');
% ylabel(['Flux (a.u.)'], 'FontSize', 14, 'FontWeight', 'bold');
ylabel(['Objective function flux value (a.u.)'], 'FontSize', 16, 'FontWeight', 'bold');
ylabel(['Maximized NADH demand flux (a.u.)'], 'FontSize', 16, 'FontWeight', 'bold');

% title(['Maximized NADH demand flux' newline 'simulated by pFBA'], 'FontSize', 14, 'FontWeight', 'bold');
% title(['Maximized NADH demand flux'], 'FontSize', 14, 'FontWeight', 'bold');
% title(['NAD^{+} -> NADH'], 'FontSize', 16, 'FontWeight', 'bold');

ylim([0 100])
% title('Network NADH Oxidation Capacity (pFBA)', 'FontSize', 16, 'FontWeight', 'bold');
grid on;
box on;

bar_handle.FaceColor = 'flat';
bar_handle.CData(1,:) = [0.5, 0.5, 0.5]; % NSD 
bar_handle.CData(2,:) = [0.9, 0.4, 0.3]; % HSD - Red
    
% Add value labels on bars
% for i = 1:n_models
%     if i == 2
%     text(i, nadh_oxidation_capacity(i) + max(nadh_oxidation_capacity)*0.04, ...
%          sprintf('x%.2f', nadh_oxidation_capacity_rel), ...
%          'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
%     end
% end

% Save figure
saveas(gcf, fullfile(output_dir, 'NADH_Oxidation_Capacity_Bar.fig'));
saveas(gcf, fullfile(output_dir, 'NADH_Oxidation_Capacity_Bar.png'));
print(gcf, fullfile(output_dir, 'NADH_Oxidation_Capacity_Bar.svg'), '-dsvg');
set(gcf, 'Units', 'inches');
set(gcf, 'PaperUnits', 'inches');

% Define figure size
width = 4;
height = 4;

% Apply size to paper
set(gcf, 'PaperPosition', [0 0 width height]);
set(gcf, 'PaperSize', [width height]);

% Export to PDF
print(gcf, fullfile(output_dir, 'NADH_Oxidation_Capacity_Bar.pdf'), '-dpdf', '-painters');


%% 4. Figure 2: Fold Change Analysis (if multiple models)
if n_models > 1
    % Calculate fold changes relative to first model (control)
    control_capacity = nadh_oxidation_capacity(1);
    fold_changes = nadh_oxidation_capacity ./ control_capacity;
    
    figure(2)
%     figure('Position', [150, 150, 800, 600]);
    bar(fold_changes, 'FaceColor', [0.8, 0.2, 0.4], 'EdgeColor', 'black', 'LineWidth', 1.5);
    hold on;
    plot([0, n_models+1], [1, 1], 'k--', 'LineWidth', 2); % Reference line at 1
    hold off;
    
    set(gca, 'XTickLabel', model_names, 'FontSize', 12, 'FontWeight', 'bold');
    ylabel('Fold Change in NADH Oxidation Capacity', 'FontSize', 14, 'FontWeight', 'bold');
    title('NADH Oxidation Capacity Fold Changes (pFBA)', 'FontSize', 16, 'FontWeight', 'bold');
    ylim([0, max(fold_changes)*1.2]);
    grid on;
    box on;
    
    % Add value labels
    for i = 1:n_models
        text(i, fold_changes(i) + max(fold_changes)*0.02, ...
             sprintf('%.2fx', fold_changes(i)), ...
             'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
    end
    
    % Save figure
    saveas(gcf, fullfile(output_dir, 'NADH_Fold_Changes.fig'));
    saveas(gcf, fullfile(output_dir, 'NADH_Fold_Changes.png'));
    print(gcf, fullfile(output_dir, 'NADH_Fold_Changes.svg'), '-dsvg');
end

%% 7. Generate Summary Table
summary_table = table(model_names, nadh_oxidation_capacity, biomass_rates, nadh_to_biomass_ratio, ...
    'VariableNames', {'ModelID', 'NADH_Oxidation_Capacity', 'Biomass_Rate', 'NADH_Efficiency'});

if n_models > 1
    fold_changes_table = nadh_oxidation_capacity ./ nadh_oxidation_capacity(1);
    summary_table.Fold_Change = fold_changes_table;
end

% Save summary table
writetable(summary_table, fullfile(output_dir, 'NADH_Analysis_Summary.csv'));

%% 8. Display Results Summary
fprintf('\n=== NADH OXIDATION CAPACITY ANALYSIS SUMMARY ===\n');
fprintf('Results saved to: %s\n\n', output_dir);

for i = 1:n_models
    fprintf('Model: %s\n', model_names{i});
    fprintf('  NADH Max Capacity: %.3f a.u.\n', nadh_oxidation_capacity(i));
    fprintf('  Biomass Rate: %.3f 1/h\n', biomass_rates(i));
    if n_models > 1
        fprintf('  Fold Change: %.2fx\n', nadh_oxidation_capacity(i)/nadh_oxidation_capacity(1));
    end
    fprintf('\n');
end

fprintf('\nVisualization complete!\n');