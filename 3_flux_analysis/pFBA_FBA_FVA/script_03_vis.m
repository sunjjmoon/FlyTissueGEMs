%% script_03_vis

clc;
close all
clear;
% clear all;

%% Setup directories
current_path = pwd;
results_dir = fullfile(current_path, '03_vis');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
end

%% Load analysis results
load(fullfile(current_path,'01_results','model_constrained_out.mat'));
load(fullfile(current_path,'02_comp_fba_pFBA','FBA_pFBA_comparison.mat'));
load(fullfile(current_path,'02_1_comp_nad', 'NAD_flux_analysis.mat'));

%% Fig 1. Flux Comparison - FBA vs pFBA
figure(1);
clf;

% Prepare data
models = {'NSD', 'HSD'};
fba_total = [comparison_results{1}.FBA.total_flux, comparison_results{2}.FBA.total_flux];
pfba_total = [comparison_results{1}.pFBA.total_flux, comparison_results{2}.pFBA.total_flux];

fba_total_avg = [comparison_results{1}.FBA.total_flux_avg, comparison_results{2}.FBA.total_flux_avg];
pfba_total_avg = [comparison_results{1}.pFBA.total_flux_avg, comparison_results{2}.pFBA.total_flux_avg];

% Total flux comparison
fig = figure(1);
subplot(1,2,1);
bar_data = [fba_total; pfba_total]';
b = bar(bar_data, 'grouped');
set(gca, 'XTickLabel', models);
ylabel('Total flux (a.u.)','fontsize',14,'fontweight','bold');
% title('Total Flux: FBA vs pFBA');
legend({'FBA', 'pFBA'}, 'Location', 'northeast');
grid on;


% Flux reduction percentages
subplot(1,2,2);
bar_data_2 = [fba_total_avg; pfba_total_avg]';

total_reduction = [comparison_results{1}.flux_reduction_percent, comparison_results{2}.flux_reduction_percent];
flux_reduction_percent_avg = [comparison_results{1}.flux_reduction_percent_avg, comparison_results{2}.flux_reduction_percent_avg];

b2 = bar(bar_data_2);
set(gca, 'XTickLabel', models);
ylabel('Average flux (a.u.)','fontsize',14,'fontweight','bold');
legend({'FBA', 'pFBA'}, 'Location', 'northeast');
% title('Total Flux Reduction');
% b2.FaceColor = [0.6 0.8 0.4]; % Green color
grid on;
ylim([0 80]);


disp(['FLux reduction_percentage NSD: '  num2str(flux_reduction_percent_avg(1))])
disp(['FLux reduction_percentage HSD: '  num2str(flux_reduction_percent_avg(2))])

savefig(fig,fullfile(results_dir, ['01_flux_comp.fig']));
saveas(fig,fullfile(results_dir, ['01_flux_comp.png']), 'png');
% --- Save high-quality PDF version ---
set(fig, 'Units', 'inches');
screenposition = get(fig, 'Position');
set(fig, 'PaperPosition', [0 0 screenposition(3:4)], ...
         'PaperSize', [screenposition(3:4)], ...
         'PaperUnits', 'inches', ...
         'PaperOrientation', 'portrait');

% Export using vector graphics
print(fig, fullfile(results_dir, '01_flux_comp.pdf'), '-dpdf', '-painters');


%% Fig 2. NAD Flux Analysis
fig2 = figure(2);
clf;

fba_nad = [nad_analysis{1}.fba_nad_total, nad_analysis{2}.fba_nad_total];
pfba_nad = [nad_analysis{1}.pfba_nad_total, nad_analysis{2}.pfba_nad_total];

% NAD flux comparison
subplot(1,2,1);
bar_data_nad = [fba_nad; pfba_nad]';
b_nad = bar(bar_data_nad, 'grouped');
set(gca, 'XTickLabel', models);
ylabel('total NAD(P)-dependent flux (a.u.)','fontsize',14,'fontweight','bold');
legend({'FBA', 'pFBA'}, 'Location', 'northeast');
ylim([0 2.5e4]);

grid on;

% Colors
% b_nad(1).FaceColor = [0.2 0.8 0.6]; % Green for FBA NAD
% b_nad(2).FaceColor = [0.8 0.2 0.6]; % Pink for pFBA NAD

fba_nad_avg = [nad_analysis{1}.fba_nad_total_avg, nad_analysis{2}.fba_nad_total_avg];
pfba_nad_avg = [nad_analysis{1}.pfba_nad_total_avg, nad_analysis{2}.pfba_nad_total_avg];

% NAD flux comparison
subplot(1,2,2);
bar_data_nad_avg = [fba_nad_avg; pfba_nad_avg]';
b_nad_avg = bar(bar_data_nad_avg, 'grouped');
set(gca, 'XTickLabel', models);
ylabel('Average NAD(P)-dependent flux (a.u.)','fontsize',14,'fontweight','bold');
legend({'FBA', 'pFBA'}, 'Location', 'northeast');
ylim([0 80]);

grid on;

nad_reduction = [nad_analysis{1}.nad_flux_reduction, nad_analysis{2}.nad_flux_reduction];
nad_reduction_avg = [nad_analysis{1}.nad_flux_reduction_avg, nad_analysis{2}.nad_flux_reduction_avg];

disp([newline 'NAD FLux reduction_percentage avg NSD: '  num2str(nad_reduction_avg(1))])
disp(['NAD FLux reduction_percentage avg HSD: '  num2str(nad_reduction_avg(2))])

savefig(fig2,fullfile(results_dir, ['02_flux_comp_nad.fig']));
saveas(fig2,fullfile(results_dir, ['02_flux_comp_nad.png']), 'png');
% --- Save high-quality PDF version ---
set(fig2, 'Units', 'inches');
screenposition = get(fig2, 'Position');
set(fig2, 'PaperPosition', [0 0 screenposition(3:4)], ...
         'PaperSize', [screenposition(3:4)], ...
         'PaperUnits', 'inches', ...
         'PaperOrientation', 'portrait');

% Export using vector graphics
print(fig2, fullfile(results_dir, '02_flux_comp_nad.pdf'), '-dpdf', '-painters');


%% Fig 3. NAD reactions 
% NAD flux reduction
fig3= figure(3);
% NAD as percentage of total flux
subplot(1,2,1);
nad_percent_fba = (fba_nad ./ fba_total) * 100;
nad_percent_pfba = (pfba_nad ./ pfba_total) * 100;
bar_data3 = [nad_percent_fba; nad_percent_pfba]';
b3 = bar(bar_data3, 'grouped');
set(gca, 'XTickLabel', models);
ylabel('Proportion of NAD(P)-depdent flux (%)', 'FontSize', 14, 'FontWeight', 'bold');
% title('NAD Flux Contribution');
legend({'FBA', 'pFBA'}, 'Location', 'northeast');
ylim([0 30])
grid on;

% Add percentage labels
for i = 1:length(models)
    text(i-0.15, nad_percent_fba(i)+1, sprintf('%.1f%%', nad_percent_fba(i)), ...
        'HorizontalAlignment', 'center');
    text(i+0.15, nad_percent_pfba(i)+1, sprintf('%.1f%%', nad_percent_pfba(i)), ...
        'HorizontalAlignment', 'center');
end

% Cycling reactions count
subplot(1,2,2);
cycling_counts = [length(nad_analysis{1}.cycling_candidates), length(nad_analysis{2}.cycling_candidates)];
b4 = bar(cycling_counts);
set(gca, 'XTickLabel', models);
ylabel('Number of reactions', 'FontSize', 14, 'FontWeight', 'bold');
ylim([0 20])
b4.FaceColor = [0.7, 0.7, 0.7]; %gray
grid on;

savefig(fig3,fullfile(results_dir, ['03_flux_comp_nad_proportion.fig']));
saveas(fig3,fullfile(results_dir, ['03_flux_comp_nad_numb_rxns.png']), 'png');
% --- Save high-quality PDF version ---
set(fig3, 'Units', 'inches');
screenposition = get(fig3, 'Position');
set(fig3, 'PaperPosition', [0 0 screenposition(3:4)], ...
         'PaperSize', [screenposition(3:4)], ...
         'PaperUnits', 'inches', ...
         'PaperOrientation', 'portrait');

% Export using vector graphics
print(fig3, fullfile(results_dir, '03_flux_comp_nad_numb_rxns.pdf'), '-dpdf', '-painters');

%% Fig 4: Individual Model Analysis - Top NAD Reactions
fig4= figure(4);
clf;

for model_idx = 1:2
    model_id = nad_analysis{model_idx}.modelID;
    
    % Get top 15 NAD reactions with highest flux changes
    flux_changes = abs(nad_analysis{model_idx}.fba_nad_fluxes) - abs(nad_analysis{model_idx}.pfba_nad_fluxes);
    [sorted_changes, sorted_idx] = sort(flux_changes, 'descend');
    
    % Select top reactions for visualization
    n_top = min(15, sum(sorted_changes > 5)); % Show reactions with >5 flux change
    if n_top == 0
        n_top = min(10, length(sorted_idx)); % Show top 10 if no significant changes
    end
    
    subplot(1, 2, model_idx);
    
    if n_top > 0
        top_idx = sorted_idx(1:n_top);
        rxn_names_short = cell(n_top, 1);
        
        % Shorten reaction names for display
        for i = 1:n_top
            full_name = nad_analysis{model_idx}.nad_rxn_names{top_idx(i)};
            if length(full_name) > 10
                rxn_names_short{i} = [full_name(1:7) '...'];
            else
                rxn_names_short{i} = full_name;
            end
        end
        
        fba_values = abs(nad_analysis{model_idx}.fba_nad_fluxes(top_idx));
        pfba_values = abs(nad_analysis{model_idx}.pfba_nad_fluxes(top_idx));
        
        % Create horizontal bar chart
        y_pos = 1:n_top;
        barh(y_pos - 0.2, fba_values, 0.4, 'DisplayName', 'FBA');
        hold on;
        barh(y_pos + 0.2, pfba_values, 0.4, 'DisplayName', 'pFBA');
        xlim([0 1100])
        
        set(gca, 'YTick', y_pos, 'YTickLabel', rxn_names_short);
        xlabel('Flux (a.u.)');
        title(sprintf('%s: Top NAD Reactions', model_id));
%         title(model_id)
        legend('Location', 'southeast');
        grid on;
        
        % Flip y-axis to show highest changes at top
        set(gca, 'YDir', 'reverse');
    else
        text(0.5, 0.5, 'No significant NAD flux changes', ...
            'HorizontalAlignment', 'center', 'FontSize', 14);
        title(sprintf('%s: NAD Reactions', model_id));
%         title(model_id)
    end
end

savefig(fig4,fullfile(results_dir, ['04_flux_comp_nad_top.fig']));
saveas(fig4,fullfile(results_dir, ['04_flux_comp_nad_top.png']), 'png');
% --- Save high-quality PDF version ---
set(fig4, 'Units', 'inches');
screenposition = get(fig4, 'Position');
set(fig4, 'PaperPosition', [0 0 screenposition(3:4)], ...
         'PaperSize', [screenposition(3:4)], ...
         'PaperUnits', 'inches', ...
         'PaperOrientation', 'portrait');

% Export using vector graphics
print(fig4, fullfile(results_dir, '04_flux_comp_nad_top.pdf'), '-dpdf', '-painters');


%% 4: NAD Flux Distribution (Histograms)
fig5 = figure(5); 
clf;

for model_idx = 1:2
    model_id = nad_analysis{model_idx}.modelID;

    subplot(1,2,model_idx);

    fba_fluxes  = abs(nad_analysis{model_idx}.fba_nad_fluxes);
    pfba_fluxes = abs(nad_analysis{model_idx}.pfba_nad_fluxes);

    % Remove near-zero
    fba_nonzero  = fba_fluxes(fba_fluxes > 1e-6);
    pfba_nonzero = pfba_fluxes(pfba_fluxes > 1e-6);

    % Data + categorical groups
    flux_values = [fba_nonzero(:); pfba_nonzero(:)];
    groups = [repmat("FBA",  numel(fba_nonzero), 1);
              repmat("pFBA", numel(pfba_nonzero), 1)];
    groups = categorical(groups);
    groups = reordercats(groups, ["FBA","pFBA"]);   % ensure order

    boxchart(groups, flux_values);                  % <- groups must be categorical/numeric
    set(gca, 'YScale', 'log');
    ylabel('Flux (a.u.)', 'FontSize', 14, 'FontWeight', 'bold');
    grid on;
    set(gca, 'XTickLabel', {'FBA','pFBA'}, 'FontSize', 12);
    title(sprintf('%s: NAD Flux Distribution', model_id), 'FontSize', 12, 'FontWeight', 'bold');


end

savefig(fig5,fullfile(results_dir, ['05_flux_comp_nad_box.fig']));
saveas(fig5,fullfile(results_dir, ['05_flux_comp_nad_box.png']), 'png');
% --- Save high-quality PDF version ---
set(fig5, 'Units', 'inches');
screenposition = get(fig5, 'Position');
set(fig5, 'PaperPosition', [0 0 screenposition(3:4)], ...
         'PaperSize', [screenposition(3:4)], ...
         'PaperUnits', 'inches', ...
         'PaperOrientation', 'portrait');

% Export using vector graphics
print(fig5, fullfile(results_dir, '05_flux_comp_nad_box.pdf'), '-dpdf', '-painters');


%% Check reactions with lower fluxes
% Initialize cycling data for export
cycling_export_data = {};
export_row = 1;
cycling_export_data{export_row, 1} = 'Model';
cycling_export_data{export_row, 2} = 'Reaction';
cycling_export_data{export_row, 3} = 'FBA_Flux';
cycling_export_data{export_row, 4} = 'pFBA_Flux';
cycling_export_data{export_row, 5} = 'Flux_Reduction';
cycling_export_data{export_row, 6} = 'Percent_Reduction';
export_row = export_row + 1;

for model_idx = 1:2
    model_id = nad_analysis{model_idx}.modelID;
    
    subplot(1, 2, model_idx);
    
    fba_fluxes = abs(nad_analysis{model_idx}.fba_nad_fluxes);
    pfba_fluxes = abs(nad_analysis{model_idx}.pfba_nad_fluxes);
    
    
    % Highlight those with lower fluxes in pFBA compared to FBA
    if ~isempty(nad_analysis{model_idx}.cycling_candidates)
        cycling_idx = nad_analysis{model_idx}.cycling_candidates;
%         scatter(fba_fluxes(cycling_idx), pfba_fluxes(cycling_idx), 80, 'r', ...
%             'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
        
        % Export cycling data
        for i = 1:length(cycling_idx)
            rxn_idx = cycling_idx(i);
            cycling_export_data{export_row, 1} = model_id;
            cycling_export_data{export_row, 2} = nad_analysis{model_idx}.nad_rxn_names{rxn_idx};
            cycling_export_data{export_row, 3} = nad_analysis{model_idx}.fba_nad_fluxes(rxn_idx);
            cycling_export_data{export_row, 4} = nad_analysis{model_idx}.pfba_nad_fluxes(rxn_idx);
            cycling_export_data{export_row, 5} = fba_fluxes(rxn_idx) - pfba_fluxes(rxn_idx);
            cycling_export_data{export_row, 6} = ((fba_fluxes(rxn_idx) - pfba_fluxes(rxn_idx)) / fba_fluxes(rxn_idx)) * 100;
            export_row = export_row + 1;
        end
    end
    

end

%% Fig 6: Number of reactions with reduced flux (pFBA vs FBA)
fig6 = figure(6);
clf;

% Count reactions with reduced flux for each model
reduced_counts = zeros(1, 2);
for model_idx = 1:2
    fba_fluxes = abs(comparison_results{model_idx}.FBA.solution.x);
    pfba_fluxes = abs(comparison_results{model_idx}.pFBA.solution.x);
    reduced_counts(model_idx) = sum(pfba_fluxes < fba_fluxes);
end

b6 = bar(reduced_counts);
set(gca, 'XTickLabel', models);
ylabel('Number of reactions with reduced flux', 'FontSize', 14, 'FontWeight', 'bold');
b6.FaceColor = [0.7 0.7 0.7]; % Pink color
grid on;

% % Add count labels on bars
% for i = 1:length(models)
%     text(i, reduced_counts(i)+20, sprintf('%d', reduced_counts(i)), ...
%         'HorizontalAlignment', 'center', 'FontWeight', 'bold');
% end

savefig(fig6, fullfile(results_dir, '06_reduced_flux_count.fig'));
saveas(fig6, fullfile(results_dir, '06_reduced_flux_count.png'), 'png');
% --- Save high-quality PDF version ---
set(fig6, 'Units', 'inches');
screenposition = get(fig6, 'Position');
set(fig6, 'PaperPosition', [0 0 screenposition(3:4)], ...
         'PaperSize', [screenposition(3:4)], ...
         'PaperUnits', 'inches', ...
         'PaperOrientation', 'portrait');

% Export using vector graphics
print(fig6, fullfile(results_dir, '06_reduced_flux_count.pdf'), '-dpdf', '-painters');

%% Fig 7: Percentage of reactions with reduced flux (pFBA vs FBA)
fig7 = figure(7);
clf;

% Calculate percentage of reactions with reduced flux
reduced_percent = zeros(1, 2);
for model_idx = 1:2
    fba_fluxes = abs(comparison_results{model_idx}.FBA.solution.x);
    pfba_fluxes = abs(comparison_results{model_idx}.pFBA.solution.x);
    total_reactions = length(fba_fluxes);
    reduced_count = sum(pfba_fluxes < fba_fluxes);
    reduced_percent(model_idx) = (reduced_count / total_reactions) * 100;
end

b7 = bar(reduced_percent);
set(gca, 'XTickLabel', models);
ylabel('Percentage of reactions with reduced flux (%)', 'FontSize', 14, 'FontWeight', 'bold');
b6.FaceColor = [0.7 0.7 0.7]; % Pink color
grid on;
ylim([0 100]);

savefig(fig7, fullfile(results_dir, '07_reduced_flux_percent.fig'));
saveas(fig7, fullfile(results_dir, '07_reduced_flux_percent.png'), 'png');
% --- Save high-quality PDF version ---
set(fig7, 'Units', 'inches');
screenposition = get(fig7, 'Position');
set(fig7, 'PaperPosition', [0 0 screenposition(3:4)], ...
         'PaperSize', [screenposition(3:4)], ...
         'PaperUnits', 'inches', ...
         'PaperOrientation', 'portrait');

% Export using vector graphics
print(fig7, fullfile(results_dir, '07_reduced_flux_percent.pdf'), '-dpdf', '-painters');


%% Fig 8: Categorize all reactions (reduced, unchanged, zero)
fig8 = figure(8);
clf;

% Calculate categories for each model
categories = {'Reduced', 'Unchanged', 'Zero in both'};
category_counts = zeros(2, 3); % 2 models x 3 categories

for model_idx = 1:2
    fba_fluxes = abs(comparison_results{model_idx}.FBA.solution.x);
    pfba_fluxes = abs(comparison_results{model_idx}.pFBA.solution.x);
    
    % Define thresholds
    flux_threshold = 1e-6; % Consider fluxes below this as zero
    
    % Categorize reactions
    reduced_idx = pfba_fluxes < fba_fluxes & fba_fluxes > flux_threshold;
    zero_both_idx = fba_fluxes <= flux_threshold & pfba_fluxes <= flux_threshold;
    unchanged_idx = ~reduced_idx & ~zero_both_idx;
    
    category_counts(model_idx, 1) = sum(reduced_idx);      % Reduced
    category_counts(model_idx, 2) = sum(unchanged_idx);    % Unchanged  
    category_counts(model_idx, 3) = sum(zero_both_idx);    % Zero in both
end

% Create stacked bar chart
b8 = bar(category_counts, 'stacked');
set(gca, 'XTickLabel', models);
ylabel('Number of reactions', 'FontSize', 14, 'FontWeight', 'bold');
legend(categories, 'Location', 'northeast');
grid on;

% Add percentage labels
for model_idx = 1:2
    total_rxns = sum(category_counts(model_idx, :));
    cumsum_counts = cumsum(category_counts(model_idx, :));
    
    for cat_idx = 1:3
        if category_counts(model_idx, cat_idx) > 0
            y_pos = cumsum_counts(cat_idx) - category_counts(model_idx, cat_idx)/2;
            percent = (category_counts(model_idx, cat_idx) / total_rxns) * 100;
            text(model_idx, y_pos, sprintf('%.1f%%', percent), ...
                'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'white');
        end
    end
end

savefig(fig8, fullfile(results_dir, '08_flux_categories.fig'));
saveas(fig8, fullfile(results_dir, '08_flux_categories.png'), 'png');
% --- Save high-quality PDF version ---
set(fig8, 'Units', 'inches');
screenposition = get(fig8, 'Position');
set(fig8, 'PaperPosition', [0 0 screenposition(3:4)], ...
         'PaperSize', [screenposition(3:4)], ...
         'PaperUnits', 'inches', ...
         'PaperOrientation', 'portrait');

% Export using vector graphics
print(fig8, fullfile(results_dir, '08_flux_categories.pdf'), '-dpdf', '-painters');


%% Fig 9: Percentage of reactions
fig9 = figure(9);
clf;

% Calculate percentages for active reactions only (excluding zero in both)
active_categories = {'reduced', 'unchanged'};
active_percent = zeros(2, 2); % 2 models x 2 categories

for model_idx = 1:2
    fba_fluxes = abs(comparison_results{model_idx}.FBA.solution.x);
    pfba_fluxes = abs(comparison_results{model_idx}.pFBA.solution.x);
    
    flux_threshold = 1e-6;
    
    % Only consider active reactions (non-zero in at least one method)
    active_idx = fba_fluxes > flux_threshold | pfba_fluxes > flux_threshold;
    
    fba_active = fba_fluxes(active_idx);
    pfba_active = pfba_fluxes(active_idx);
    
    % Categorize active reactions
    reduced_count = sum(pfba_active < fba_active);
    unchanged_count = sum(pfba_active >= fba_active);
    total_active = reduced_count + unchanged_count;
    
    active_percent(model_idx, 1) = (reduced_count / total_active) * 100;   % Reduced
    active_percent(model_idx, 2) = (unchanged_count / total_active) * 100; % Unchanged
end

% Create bar chart
b9 = bar(active_percent, 'stacked');
set(gca, 'XTickLabel', models);
ylabel('Percentage of reactions affected by pFBA (%)', 'FontSize', 14, 'FontWeight', 'bold');
legend(active_categories, 'Location', 'eastoutside');
grid on;
ylim([0 100]);

% Add percentage labels
for model_idx = 1:2
    cumsum_percent = cumsum(active_percent(model_idx, :));
    
    for cat_idx = 1:2
        if active_percent(model_idx, cat_idx) > 0
            y_pos = cumsum_percent(cat_idx) - active_percent(model_idx, cat_idx)/2;
            text(model_idx, y_pos, sprintf('%.1f%%', active_percent(model_idx, cat_idx)), ...
                'HorizontalAlignment', 'center', 'FontWeight', 'bold', 'Color', 'white');
        end
    end
end

savefig(fig9, fullfile(results_dir, '09_reactions_affected.fig'));
saveas(fig9, fullfile(results_dir, '09_reactions_affected.png'), 'png');
% --- Save high-quality PDF version ---
set(fig9, 'Units', 'inches');
screenposition = get(fig9, 'Position');
set(fig9, 'PaperPosition', [0 0 screenposition(3:4)], ...
         'PaperSize', [screenposition(3:4)], ...
         'PaperUnits', 'inches', ...
         'PaperOrientation', 'portrait');

% Export using vector graphics
print(fig9, fullfile(results_dir, '09_reactions_affected.pdf'), '-dpdf', '-painters');

%% Export reduced flux reactions to Excel
excel_filename = fullfile(results_dir, 'reduced_flux_reactions.xlsx');

for model_idx = 1:2
    model_id = models{model_idx};
    fba_fluxes = abs(comparison_results{model_idx}.FBA.solution.x);
    pfba_fluxes = abs(comparison_results{model_idx}.pFBA.solution.x);
    
    % Find reactions with reduced flux
    reduced_idx = find(pfba_fluxes < fba_fluxes);
    
    if ~isempty(reduced_idx)
        % Create data for export
        rxn_names = model_constrained_out{model_idx}.rxns(reduced_idx);
        fba_vals = fba_fluxes(reduced_idx);
        pfba_vals = pfba_fluxes(reduced_idx);
        relative_vals = pfba_vals ./ fba_vals; % pFBA/FBA ratio
        
        % Create table
        export_table = table(rxn_names, fba_vals, pfba_vals, relative_vals, ...
            'VariableNames', {'Reaction', 'FBA_Flux', 'pFBA_Flux', 'Relative_Flux'});
        
        % Write to Excel sheet
        writetable(export_table, excel_filename, 'Sheet', model_id);
    end
end

%% Fig 11. Flux Comparison - FBA vs pFBA with Flux Distribution
figure(11);
clf;

% Prepare data
models = {'NSD', 'HSD'};
fba_total = [comparison_results{1}.FBA.total_flux, comparison_results{2}.FBA.total_flux];
pfba_total = [comparison_results{1}.pFBA.total_flux, comparison_results{2}.pFBA.total_flux];
fba_total_avg = [comparison_results{1}.FBA.total_flux_avg, comparison_results{2}.FBA.total_flux_avg];
pfba_total_avg = [comparison_results{1}.pFBA.total_flux_avg, comparison_results{2}.pFBA.total_flux_avg];

fig11 = figure(11);

% Subplot 1: Total flux comparison (bar chart)
subplot(1,2,1);
bar_data = [fba_total; pfba_total]';
b = bar(bar_data, 'grouped');
set(gca, 'XTickLabel', models);
ylabel('Total flux (a.u.)','fontsize',14,'fontweight','bold');
legend({'FBA', 'pFBA'}, 'Location', 'northeast');
grid on;

% Subplot 2: Flux distribution (boxplot)
subplot(1,2,2);
hold on;

for model_idx = 1:2
    % Get flux data for each model (use .x field from solution)
    fba_fluxes  = abs(comparison_results{model_idx}.FBA.solution.x);
    pfba_fluxes = abs(comparison_results{model_idx}.pFBA.solution.x);
    
    % Remove near-zero fluxes
    fba_nonzero  = fba_fluxes(fba_fluxes > 1e-6);
    pfba_nonzero = pfba_fluxes(pfba_fluxes > 1e-6);
    
    % Combine data with categorical groups
    flux_values = [fba_nonzero(:); pfba_nonzero(:)];
    
    % Create model-specific labels
    if model_idx == 1
        groups = [repmat("NSD-FBA",  numel(fba_nonzero), 1);
                  repmat("NSD-pFBA", numel(pfba_nonzero), 1)];
    else
        groups = [repmat("HSD-FBA",  numel(fba_nonzero), 1);
                  repmat("HSD-pFBA", numel(pfba_nonzero), 1)];
    end
    groups = categorical(groups);
    
    % Create boxplot
    boxchart(groups, flux_values);
end

set(gca, 'YScale', 'log');
ylabel('Flux (a.u.)', 'FontSize', 14, 'FontWeight', 'bold');
grid on;
% title('Flux Distribution', 'FontSize', 12, 'FontWeight', 'bold');
xtickangle(45);

% Display flux reduction percentages
flux_reduction_percent_avg = [comparison_results{1}.flux_reduction_percent_avg, comparison_results{2}.flux_reduction_percent_avg];
disp(['Flux reduction_percentage NSD: '  num2str(flux_reduction_percent_avg(1))])
disp(['Flux reduction_percentage HSD: '  num2str(flux_reduction_percent_avg(2))])

% Save figures
savefig(fig11, fullfile(results_dir, '11_flux_comp_distribution.fig'));
saveas(fig11, fullfile(results_dir, '11_flux_comp_distribution.png'), 'png');

% Save high-quality PDF version
set(fig11, 'Units', 'inches');
screenposition = get(fig11, 'Position');
set(fig11, 'PaperPosition', [0 0 screenposition(3:4)], ...
         'PaperSize', [screenposition(3:4)], ...
         'PaperUnits', 'inches', ...
         'PaperOrientation', 'portrait');
print(fig11, fullfile(results_dir, '11_flux_comp_distribution.pdf'), '-dpdf', '-painters');

%% Export cycling reactions to CSV
cycling_filename = fullfile(results_dir, 'cycling_reactions.csv');
if size(cycling_export_data, 1) > 1 % Check if we have data beyond headers
    % Convert cell array to table for better CSV export
    cycling_table = cell2table(cycling_export_data(2:end, :), 'VariableNames', cycling_export_data(1, :));
    writetable(cycling_table, cycling_filename);
%     fprintf('Cycling reactions exported to: %s\n', cycling_filename);
else
    fprintf('No cycling reactions found to export.\n');
end


%% Export cycling reactions to CSV with detailed information


%% Export cycling reactions to CSV with detailed information (WITH EQUATIONS)
cycling_filename = fullfile(results_dir, 'cycling_reactions_detailed.csv');

if size(cycling_export_data, 1) > 1 % Check if we have data beyond headers
    
    % Enhanced export with gene and reaction information
    enhanced_cycling_data = {};
    export_row = 1;
    
    % Headers
    enhanced_cycling_data{export_row, 1} = 'Model';
    enhanced_cycling_data{export_row, 2} = 'Reaction_ID';
    enhanced_cycling_data{export_row, 3} = 'Reaction_Name';
    enhanced_cycling_data{export_row, 4} = 'Reaction_Equation';
    enhanced_cycling_data{export_row, 5} = 'Subsystem';
    enhanced_cycling_data{export_row, 6} = 'Genes';
    enhanced_cycling_data{export_row, 7} = 'Gene_Rules';
    enhanced_cycling_data{export_row, 8} = 'Reversible';
    enhanced_cycling_data{export_row, 9} = 'FBA_Flux';
    enhanced_cycling_data{export_row, 10} = 'pFBA_Flux';
    enhanced_cycling_data{export_row, 11} = 'Flux_Reduction';
    enhanced_cycling_data{export_row, 12} = 'Percent_Reduction';
    export_row = export_row + 1;
    
    % Process each model
    for model_idx = 1:2
        model = model_constrained_out{model_idx};
        model_id = nad_analysis{model_idx}.modelID;
        
        if ~isempty(nad_analysis{model_idx}.cycling_candidates)
            cycling_idx = nad_analysis{model_idx}.cycling_candidates;
            
            for i = 1:length(cycling_idx)
                rxn_idx = cycling_idx(i);
                
                % Get reaction index in NAD analysis
                nad_rxn_idx = nad_analysis{model_idx}.nad_rxn_indices(rxn_idx);
                
                % Get reaction ID
                rxn_id = model.rxns{nad_rxn_idx};
                
                % Basic info
                enhanced_cycling_data{export_row, 1} = model_id;
                enhanced_cycling_data{export_row, 2} = rxn_id;
                enhanced_cycling_data{export_row, 3} = model.rxnNames{nad_rxn_idx};
                
                % Reaction equation
                enhanced_cycling_data{export_row, 4} = makeReactionEquation(model, nad_rxn_idx);
                
                % Subsystem
                if ~isempty(model.subSystems{nad_rxn_idx})
                    if iscell(model.subSystems{nad_rxn_idx})
                        enhanced_cycling_data{export_row, 5} = strjoin(model.subSystems{nad_rxn_idx}, '; ');
                    else
                        enhanced_cycling_data{export_row, 5} = model.subSystems{nad_rxn_idx};
                    end
                else
                    enhanced_cycling_data{export_row, 5} = '';
                end
                
                % Genes associated with reaction
                gene_idx = find(model.rxnGeneMat(nad_rxn_idx, :));
                if ~isempty(gene_idx)
                    gene_list = model.genes(gene_idx);
                    enhanced_cycling_data{export_row, 6} = strjoin(gene_list, '; ');
                else
                    enhanced_cycling_data{export_row, 6} = '';
                end
                
                % Gene rules (grRules)
                if ~isempty(model.grRules{nad_rxn_idx})
                    enhanced_cycling_data{export_row, 7} = model.grRules{nad_rxn_idx};
                else
                    enhanced_cycling_data{export_row, 7} = '';
                end
                
                % Reversibility
                enhanced_cycling_data{export_row, 8} = model.rev(nad_rxn_idx);
                
                % Flux values
                fba_flux = nad_analysis{model_idx}.fba_nad_fluxes(rxn_idx);
                pfba_flux = nad_analysis{model_idx}.pfba_nad_fluxes(rxn_idx);
                
                enhanced_cycling_data{export_row, 9} = fba_flux;
                enhanced_cycling_data{export_row, 10} = pfba_flux;
                enhanced_cycling_data{export_row, 11} = abs(fba_flux) - abs(pfba_flux);
                enhanced_cycling_data{export_row, 12} = ((abs(fba_flux) - abs(pfba_flux)) / abs(fba_flux)) * 100;
                
                export_row = export_row + 1;
            end
        end
    end
    
    % Convert to table and export
    if export_row > 2
        cycling_table = cell2table(enhanced_cycling_data(2:end, :), ...
            'VariableNames', enhanced_cycling_data(1, :));
        writetable(cycling_table, cycling_filename);
        fprintf('Enhanced cycling reactions exported to: %s\n', cycling_filename);
        fprintf('Total cycling reactions exported: %d\n', export_row - 2);
    else
        fprintf('No cycling reactions found to export.\n');
    end
    
else
    fprintf('No cycling reactions found to export.\n');
end

%% Also export ALL NAD reactions (not just cycling) for reference
all_nad_filename = fullfile(results_dir, 'all_nad_reactions_detailed.csv');

all_nad_data = {};
export_row = 1;

% Headers
all_nad_data{export_row, 1} = 'Model';
all_nad_data{export_row, 2} = 'Reaction_ID';
all_nad_data{export_row, 3} = 'Reaction_Name';
all_nad_data{export_row, 4} = 'Reaction_Equation';
all_nad_data{export_row, 5} = 'Subsystem';
all_nad_data{export_row, 6} = 'Genes';
all_nad_data{export_row, 7} = 'Gene_Rules';
all_nad_data{export_row, 8} = 'Reversible';
all_nad_data{export_row, 9} = 'FBA_Flux';
all_nad_data{export_row, 10} = 'pFBA_Flux';
all_nad_data{export_row, 11} = 'Flux_Change';
all_nad_data{export_row, 12} = 'Abs_Flux_Change';
all_nad_data{export_row, 13} = 'Percent_Change';
all_nad_data{export_row, 14} = 'Is_Cycling_Candidate';
export_row = export_row + 1;

for model_idx = 1:2
    model = model_constrained_out{model_idx};
    model_id = nad_analysis{model_idx}.modelID;
    
    for i = 1:length(nad_analysis{model_idx}.nad_rxn_indices)
        nad_rxn_idx = nad_analysis{model_idx}.nad_rxn_indices(i);
        rxn_id = model.rxns{nad_rxn_idx};
        
        % Check if cycling candidate
        is_cycling = ismember(i, nad_analysis{model_idx}.cycling_candidates);
        
        % Basic info
        all_nad_data{export_row, 1} = model_id;
        all_nad_data{export_row, 2} = rxn_id;
        all_nad_data{export_row, 3} = model.rxnNames{nad_rxn_idx};
        
        % Reaction equation
        all_nad_data{export_row, 4} = makeReactionEquation(model, nad_rxn_idx);
        
        % Subsystem
        if ~isempty(model.subSystems{nad_rxn_idx})
            if iscell(model.subSystems{nad_rxn_idx})
                all_nad_data{export_row, 5} = strjoin(model.subSystems{nad_rxn_idx}, '; ');
            else
                all_nad_data{export_row, 5} = model.subSystems{nad_rxn_idx};
            end
        else
            all_nad_data{export_row, 5} = '';
        end
        
        % Genes
        gene_idx = find(model.rxnGeneMat(nad_rxn_idx, :));
        if ~isempty(gene_idx)
            gene_list = model.genes(gene_idx);
            all_nad_data{export_row, 6} = strjoin(gene_list, '; ');
        else
            all_nad_data{export_row, 6} = '';
        end
        
        % Gene rules
        if ~isempty(model.grRules{nad_rxn_idx})
            all_nad_data{export_row, 7} = model.grRules{nad_rxn_idx};
        else
            all_nad_data{export_row, 7} = '';
        end
        
        % Reversibility
        all_nad_data{export_row, 8} = model.rev(nad_rxn_idx);
        
        % Flux values
        fba_flux = nad_analysis{model_idx}.fba_nad_fluxes(i);
        pfba_flux = nad_analysis{model_idx}.pfba_nad_fluxes(i);
        flux_change = fba_flux - pfba_flux;
        
        all_nad_data{export_row, 9} = fba_flux;
        all_nad_data{export_row, 10} = pfba_flux;
        all_nad_data{export_row, 11} = flux_change;
        all_nad_data{export_row, 12} = abs(flux_change);
        
        if abs(fba_flux) > 1e-6
            all_nad_data{export_row, 13} = (flux_change / abs(fba_flux)) * 100;
        else
            all_nad_data{export_row, 13} = 0;
        end
        
        all_nad_data{export_row, 14} = is_cycling;
        
        export_row = export_row + 1;
    end
end

% Export all NAD reactions
all_nad_table = cell2table(all_nad_data(2:end, :), 'VariableNames', all_nad_data(1, :));
writetable(all_nad_table, all_nad_filename);
fprintf('All NAD reactions exported to: %s\n', all_nad_filename);
fprintf('Total NAD reactions: %d\n', export_row - 2);






%% Export All Analysis Results to Excel
% Add this section at the end of your existing script (Script 04 v2)

%% Comprehensive Excel Export
excel_filename = fullfile(results_dir, 'Complete_Analysis_Results.xlsx');
fprintf('\nExporting comprehensive results to: %s\n', excel_filename);

%% Sheet 1: Total Flux Comparison (Figure 1)
sheet_name = '01_Total_Flux_Comparison';
flux_table = table(models', fba_total', pfba_total', fba_total_avg', pfba_total_avg', flux_reduction_percent_avg', ...
    'VariableNames', {'Model', 'FBA_Total_Flux', 'pFBA_Total_Flux', 'FBA_Avg_Flux', 'pFBA_Avg_Flux', 'Flux_Reduction_Percent'});
writetable(flux_table, excel_filename, 'Sheet', sheet_name);
fprintf('  - %s written\n', sheet_name);

%% Sheet 2: NAD Flux Comparison (Figure 2)
sheet_name = '02_NAD_Flux_Comparison';
nad_table = table(models', fba_nad', pfba_nad', fba_nad_avg', pfba_nad_avg', nad_reduction_avg', ...
    'VariableNames', {'Model', 'FBA_NAD_Total', 'pFBA_NAD_Total', 'FBA_NAD_Avg', 'pFBA_NAD_Avg', 'NAD_Reduction_Percent'});
writetable(nad_table, excel_filename, 'Sheet', sheet_name);
fprintf('  - %s written\n', sheet_name);

%% Sheet 3: NAD Proportion Analysis (Figure 3)
sheet_name = '03_NAD_Proportion';
nad_percent_fba = (fba_nad ./ fba_total) * 100;
nad_percent_pfba = (pfba_nad ./ pfba_total) * 100;
cycling_counts = [length(nad_analysis{1}.cycling_candidates), length(nad_analysis{2}.cycling_candidates)];

proportion_table = table(models', nad_percent_fba', nad_percent_pfba', cycling_counts', ...
    'VariableNames', {'Model', 'NAD_Percent_FBA', 'NAD_Percent_pFBA', 'Cycling_Reactions_Count'});
writetable(proportion_table, excel_filename, 'Sheet', sheet_name);
fprintf('  - %s written\n', sheet_name);

%% Sheet 4: Top NAD Reactions (Figure 4)
sheet_name = '04_Top_NAD_Reactions';
top_nad_data = {};
row_idx = 1;
top_nad_data{row_idx, 1} = 'Model';
top_nad_data{row_idx, 2} = 'Rank';
top_nad_data{row_idx, 3} = 'Reaction_ID';
top_nad_data{row_idx, 4} = 'Reaction_Name';
top_nad_data{row_idx, 5} = 'FBA_Flux';
top_nad_data{row_idx, 6} = 'pFBA_Flux';
top_nad_data{row_idx, 7} = 'Flux_Change';
row_idx = row_idx + 1;

for model_idx = 1:2
    flux_changes = abs(nad_analysis{model_idx}.fba_nad_fluxes) - abs(nad_analysis{model_idx}.pfba_nad_fluxes);
    [sorted_changes, sorted_idx] = sort(flux_changes, 'descend');
    n_top = min(15, sum(sorted_changes > 5));
    if n_top == 0
        n_top = min(10, length(sorted_idx));
    end
    
    for i = 1:n_top
        idx = sorted_idx(i);
        % Get reaction ID from model using nad_rxn_indices
        rxn_model_idx = nad_analysis{model_idx}.nad_rxn_indices(idx);
        rxn_id = model_constrained_out{model_idx}.rxns{rxn_model_idx};
        
        top_nad_data{row_idx, 1} = nad_analysis{model_idx}.modelID;
        top_nad_data{row_idx, 2} = i;
        top_nad_data{row_idx, 3} = rxn_id;
        top_nad_data{row_idx, 4} = nad_analysis{model_idx}.nad_rxn_names{idx};
        top_nad_data{row_idx, 5} = abs(nad_analysis{model_idx}.fba_nad_fluxes(idx));
        top_nad_data{row_idx, 6} = abs(nad_analysis{model_idx}.pfba_nad_fluxes(idx));
        top_nad_data{row_idx, 7} = sorted_changes(i);
        row_idx = row_idx + 1;
    end
end

top_nad_table = cell2table(top_nad_data(2:end, :), 'VariableNames', top_nad_data(1, :));
writetable(top_nad_table, excel_filename, 'Sheet', sheet_name);
fprintf('  - %s written\n', sheet_name);

%% Sheet 5: NAD Flux Distribution Statistics (Figure 5)
sheet_name = '05_NAD_Flux_Distribution';
dist_data = {};
row_idx = 1;
dist_data{row_idx, 1} = 'Model';
dist_data{row_idx, 2} = 'Method';
dist_data{row_idx, 3} = 'Median_Flux';
dist_data{row_idx, 4} = 'Mean_Flux';
dist_data{row_idx, 5} = 'Std_Flux';
dist_data{row_idx, 6} = 'Q25';
dist_data{row_idx, 7} = 'Q75';
dist_data{row_idx, 8} = 'N_Reactions';
row_idx = row_idx + 1;

for model_idx = 1:2
    fba_fluxes = abs(nad_analysis{model_idx}.fba_nad_fluxes);
    pfba_fluxes = abs(nad_analysis{model_idx}.pfba_nad_fluxes);
    
    fba_nonzero = fba_fluxes(fba_fluxes > 1e-6);
    pfba_nonzero = pfba_fluxes(pfba_fluxes > 1e-6);
    
    % FBA stats
    dist_data{row_idx, 1} = nad_analysis{model_idx}.modelID;
    dist_data{row_idx, 2} = 'FBA';
    dist_data{row_idx, 3} = median(fba_nonzero);
    dist_data{row_idx, 4} = mean(fba_nonzero);
    dist_data{row_idx, 5} = std(fba_nonzero);
    dist_data{row_idx, 6} = quantile(fba_nonzero, 0.25);
    dist_data{row_idx, 7} = quantile(fba_nonzero, 0.75);
    dist_data{row_idx, 8} = length(fba_nonzero);
    row_idx = row_idx + 1;
    
    % pFBA stats
    dist_data{row_idx, 1} = nad_analysis{model_idx}.modelID;
    dist_data{row_idx, 2} = 'pFBA';
    dist_data{row_idx, 3} = median(pfba_nonzero);
    dist_data{row_idx, 4} = mean(pfba_nonzero);
    dist_data{row_idx, 5} = std(pfba_nonzero);
    dist_data{row_idx, 6} = quantile(pfba_nonzero, 0.25);
    dist_data{row_idx, 7} = quantile(pfba_nonzero, 0.75);
    dist_data{row_idx, 8} = length(pfba_nonzero);
    row_idx = row_idx + 1;
end

dist_table = cell2table(dist_data(2:end, :), 'VariableNames', dist_data(1, :));
writetable(dist_table, excel_filename, 'Sheet', sheet_name);
fprintf('  - %s written\n', sheet_name);

%% Sheet 6: Reduced Flux Count (Figure 6)
sheet_name = '06_Reduced_Flux_Count';
count_table = table(models', reduced_counts', ...
    'VariableNames', {'Model', 'Reactions_With_Reduced_Flux'});
writetable(count_table, excel_filename, 'Sheet', sheet_name);
fprintf('  - %s written\n', sheet_name);

%% Sheet 7: Reduced Flux Percentage (Figure 7)
sheet_name = '07_Reduced_Flux_Percent';
percent_table = table(models', reduced_percent', ...
    'VariableNames', {'Model', 'Percent_Reactions_Reduced'});
writetable(percent_table, excel_filename, 'Sheet', sheet_name);
fprintf('  - %s written\n', sheet_name);

%% Sheet 8: Flux Categories (Figure 8)
sheet_name = '08_Flux_Categories';
cat_table = table(models', category_counts(:,1), category_counts(:,2), category_counts(:,3), ...
    'VariableNames', {'Model', 'Reduced', 'Unchanged', 'Zero_In_Both'});
writetable(cat_table, excel_filename, 'Sheet', sheet_name);
fprintf('  - %s written\n', sheet_name);

%% Sheet 9: Active Reactions Affected (Figure 9)
sheet_name = '09_Active_Reactions';
active_table = table(models', active_percent(:,1), active_percent(:,2), ...
    'VariableNames', {'Model', 'Percent_Reduced', 'Percent_Unchanged'});
writetable(active_table, excel_filename, 'Sheet', sheet_name);
fprintf('  - %s written\n', sheet_name);

%% Sheet 10: Summary Statistics
sheet_name = '10_Summary';
summary_data = {
    'Metric', 'NSD', 'HSD';
    '', '', '';
    'Total Flux (FBA)', fba_total(1), fba_total(2);
    'Total Flux (pFBA)', pfba_total(1), pfba_total(2);
    'Total Flux Reduction (%)', flux_reduction_percent_avg(1), flux_reduction_percent_avg(2);
    '', '', '';
    'NAD Flux (FBA)', fba_nad(1), fba_nad(2);
    'NAD Flux (pFBA)', pfba_nad(1), pfba_nad(2);
    'NAD Flux Reduction (%)', nad_reduction_avg(1), nad_reduction_avg(2);
    '', '', '';
    'NAD as % of Total (FBA)', nad_percent_fba(1), nad_percent_fba(2);
    'NAD as % of Total (pFBA)', nad_percent_pfba(1), nad_percent_pfba(2);
    '', '', '';
    'Cycling Reactions Count', cycling_counts(1), cycling_counts(2);
    'Reactions with Reduced Flux', reduced_counts(1), reduced_counts(2);
    'Percent Reactions Reduced', reduced_percent(1), reduced_percent(2)
};

writecell(summary_data, excel_filename, 'Sheet', sheet_name);
fprintf('  - %s written\n', sheet_name);

%% Sheet 11: All Reactions Detailed
sheet_name = '11_All_Reactions_Detail';
all_rxn_data = {};
row_idx = 1;
all_rxn_data{row_idx, 1} = 'Model';
all_rxn_data{row_idx, 2} = 'Reaction_ID';
all_rxn_data{row_idx, 3} = 'FBA_Flux';
all_rxn_data{row_idx, 4} = 'pFBA_Flux';
all_rxn_data{row_idx, 5} = 'Flux_Change';
all_rxn_data{row_idx, 6} = 'Percent_Change';
all_rxn_data{row_idx, 7} = 'Category';
row_idx = row_idx + 1;

for model_idx = 1:2
    fba_fluxes = abs(comparison_results{model_idx}.FBA.solution.x);
    pfba_fluxes = abs(comparison_results{model_idx}.pFBA.solution.x);
    rxn_ids = model_constrained_out{model_idx}.rxns;
    
    flux_threshold = 1e-6;
    
    for rxn_idx = 1:length(rxn_ids)
        fba_val = fba_fluxes(rxn_idx);
        pfba_val = pfba_fluxes(rxn_idx);
        
        % Categorize
        if fba_val <= flux_threshold && pfba_val <= flux_threshold
            category = 'Zero_Both';
        elseif pfba_val < fba_val
            category = 'Reduced';
        else
            category = 'Unchanged';
        end
        
        flux_change = fba_val - pfba_val;
        if fba_val > flux_threshold
            percent_change = (flux_change / fba_val) * 100;
        else
            percent_change = 0;
        end
        
        all_rxn_data{row_idx, 1} = models{model_idx};
        all_rxn_data{row_idx, 2} = rxn_ids{rxn_idx};
        all_rxn_data{row_idx, 3} = fba_val;
        all_rxn_data{row_idx, 4} = pfba_val;
        all_rxn_data{row_idx, 5} = flux_change;
        all_rxn_data{row_idx, 6} = percent_change;
        all_rxn_data{row_idx, 7} = category;
        row_idx = row_idx + 1;
    end
end

all_rxn_table = cell2table(all_rxn_data(2:end, :), 'VariableNames', all_rxn_data(1, :));
writetable(all_rxn_table, excel_filename, 'Sheet', sheet_name);
fprintf('  - %s written (this may take a moment)\n', sheet_name);

fprintf('\n✓ Export complete! File saved as: %s\n', excel_filename);
fprintf('Total sheets created: 11\n\n');

%% Helper function to construct reaction equation
function equation = makeReactionEquation(model, rxn_idx)
    % Get stoichiometry for this reaction
    stoich = full(model.S(:, rxn_idx)); % Convert sparse to full
    
    % Find reactants (negative stoich) and products (positive stoich)
    reactant_idx = find(stoich < 0);
    product_idx = find(stoich > 0);
    
    % Build reactants string
    reactants = {};
    for i = 1:length(reactant_idx)
        met_idx = reactant_idx(i);
        coeff = abs(stoich(met_idx));
        
        % Use metabolite NAME instead of ID
        met_name = model.metNames{met_idx};
        % Also keep compartment info from met ID
        met_id = model.mets{met_idx};
        % Extract compartment (everything after last '[')
        comp_start = find(met_id == '[', 1, 'last');
        if ~isempty(comp_start)
            compartment = met_id(comp_start:end);
        else
            compartment = '';
        end
        
        full_met_name = [met_name compartment];
        
        if abs(coeff - 1) < 1e-10  % Check if coefficient is 1
            reactants{end+1} = full_met_name;
        else
            reactants{end+1} = sprintf('%.2g %s', coeff, full_met_name);
        end
    end
    
    % Build products string
    products = {};
    for i = 1:length(product_idx)
        met_idx = product_idx(i);
        coeff = stoich(met_idx);
        
        % Use metabolite NAME instead of ID
        met_name = model.metNames{met_idx};
        % Also keep compartment info from met ID
        met_id = model.mets{met_idx};
        % Extract compartment (everything after last '[')
        comp_start = find(met_id == '[', 1, 'last');
        if ~isempty(comp_start)
            compartment = met_id(comp_start:end);
        else
            compartment = '';
        end
        
        full_met_name = [met_name compartment];
        
        if abs(coeff - 1) < 1e-10  % Check if coefficient is 1
            products{end+1} = full_met_name;
        else
            products{end+1} = sprintf('%.2g %s', coeff, full_met_name);
        end
    end
    
    % Combine into equation
    if isempty(reactants)
        reactants_str = '';
    else
        reactants_str = strjoin(reactants, ' + ');
    end
    
    if isempty(products)
        products_str = '';
    else
        products_str = strjoin(products, ' + ');
    end
    
    % Determine arrow based on reversibility
    if model.rev(rxn_idx) == 1
        arrow = ' <=> ';
    else
        arrow = ' => ';
    end
    
    equation = [reactants_str arrow products_str];
end