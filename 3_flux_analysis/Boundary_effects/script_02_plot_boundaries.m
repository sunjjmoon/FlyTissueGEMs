%% Script_02_plot boundaries
clc; 
clear; 
close all;

%% Set up parameters
bound_conditions = [1000, 10000, 50000];
pathway = pwd;

% Output
analysis_dir = fullfile(pathway,'results', '04_eval_simplified');
if ~exist(analysis_dir, 'dir')
    mkdir(analysis_dir);
end

saturation_val = bound_conditions*0.01;  % 1% of the boundary

%% 1. Load sampling data
comparison_data = struct();
model_names = {'NSD', 'HSD'};
saturation_percentage = 0.99;

for bound_idx = 1:length(bound_conditions)
    current_bound = bound_conditions(bound_idx);
    dir_sampling = fullfile('results', sprintf('03_FSA_bounds_%d', current_bound));

    fprintf('Loading data for boundary ±%d...\n', current_bound);

    for model_idx = 1:length(model_names)
        model_name = model_names{model_idx};
        sample_file = fullfile(pathway, dir_sampling, [model_name, '.mat']);

        if exist(sample_file, 'file')
            data = load(sample_file);

            if isfield(data, 'samples')
                samples_struct = data.samples;
                
                % Extract sampling points
                if isstruct(samples_struct) && isfield(samples_struct, 'points')
                    samples = samples_struct.points;
                    if issparse(samples)
                        samples = full(samples);
                    end
                else
                    continue;
                end

                % Store statistics
                bound_field = sprintf('bound_%d', current_bound);
                max_vals    = max(samples, [], 2);
                min_vals    = min(samples, [], 2);
                median_vals = median(samples, 2);
                n_total     = size(samples, 1);

                saturation_val_in = saturation_val(bound_idx);
                
                % Absolute saturation
                abs_upper = (max_vals >= (current_bound - saturation_val_in));
                abs_lower = (min_vals <= (-current_bound + saturation_val_in));

                % Relative saturation
                rel_tol   = current_bound * (1 - saturation_percentage);
                rel_upper = (median_vals >= ( current_bound - rel_tol));
                rel_lower = (median_vals <= (-current_bound + rel_tol));

                % Store percentages only
                comparison_data.(bound_field).(model_name).abs_upper_pct = 100 * sum(abs_upper) / n_total;
                comparison_data.(bound_field).(model_name).abs_lower_pct = 100 * sum(abs_lower) / n_total;
                comparison_data.(bound_field).(model_name).rel_upper_pct = 100 * sum(rel_upper) / n_total;
                comparison_data.(bound_field).(model_name).rel_lower_pct = 100 * sum(rel_lower) / n_total;
                comparison_data.(bound_field).(model_name).total_reactions = n_total;
            end
        end
    end
    fprintf('\n');
end

%% Extract data for plotting
nsd_abs_upper = [];
nsd_abs_lower = [];
nsd_rel_upper = [];
nsd_rel_lower = [];

hsd_abs_upper = [];
hsd_abs_lower = [];
hsd_rel_upper = [];
hsd_rel_lower = [];

for bound_idx = 1:length(bound_conditions)
    bound_field = sprintf('bound_%d', bound_conditions(bound_idx));
    
    if isfield(comparison_data, bound_field)
        if isfield(comparison_data.(bound_field), 'NSD')
            nsd_abs_upper = [nsd_abs_upper, comparison_data.(bound_field).NSD.abs_upper_pct];
            nsd_abs_lower = [nsd_abs_lower, comparison_data.(bound_field).NSD.abs_lower_pct];
            nsd_rel_upper = [nsd_rel_upper, comparison_data.(bound_field).NSD.rel_upper_pct];
            nsd_rel_lower = [nsd_rel_lower, comparison_data.(bound_field).NSD.rel_lower_pct];
        end
        
        if isfield(comparison_data.(bound_field), 'HSD')
            hsd_abs_upper = [hsd_abs_upper, comparison_data.(bound_field).HSD.abs_upper_pct];
            hsd_abs_lower = [hsd_abs_lower, comparison_data.(bound_field).HSD.abs_lower_pct];
            hsd_rel_upper = [hsd_rel_upper, comparison_data.(bound_field).HSD.rel_upper_pct];
            hsd_rel_lower = [hsd_rel_lower, comparison_data.(bound_field).HSD.rel_lower_pct];
        end
    end
end

%% KEY FIGURE 1: Combined Absolute Saturation (NSD + HSD)
fig1 = figure(1);
clf;
% set(fig1, 'Position', [100, 100, 800, 600]);

% Prepare data: [NSD_upper, NSD_lower, HSD_upper, HSD_lower] x 3 boundaries
data_matrix = [nsd_abs_upper', nsd_abs_lower', hsd_abs_upper', hsd_abs_lower'];

b = bar(data_matrix, 'grouped');
set(gca, 'XTickLabel', {'±1,000', '±10,000', '±50,000'}, 'FontSize', 14);
xlabel('Boundary condition (a.u.)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Percentage of saturated reactions (%)', 'FontSize', 14, 'FontWeight', 'bold');
% title(['Reactions approaching boundaries' newline '(|max/min flux| ≥99% of boundary)'], 'FontSize', 14, 'FontWeight', 'bold');
title(['(|max/min flux| ≥99% of boundary)'], 'FontSize', 18, 'FontWeight', 'bold');
legend({'NSD Upper', 'NSD Lower', 'HSD Upper', 'HSD Lower'}, ...
    'Location', 'northeast', 'FontSize', 8);
grid on;
ylim([0 3]);

% Color scheme
b(1).FaceColor = [0.2 0.4 0.8];  % NSD upper - blue
b(2).FaceColor = [0.4 0.6 1.0];  % NSD lower - light blue
b(3).FaceColor = [0.8 0.4 0.2];  % HSD upper - orange
b(4).FaceColor = [1.0 0.6 0.4];  % HSD lower - light orange

% Save
saveas(fig1, fullfile(analysis_dir, 'Fig_Boundary_Absolute_Saturation.png'));
% saveas(fig1, fullfile(analysis_dir, 'Fig_Boundary_Absolute_Saturation.pdf'));
savefig(fig1, fullfile(analysis_dir, 'Fig_Boundary_Absolute_Saturation.fig'));
% Set up high-resolution vector PDF export
set(fig1, 'Units', 'inches');
screenposition = get(fig1, 'Position');
set(fig1, ...
    'PaperPosition', [0 0 screenposition(3:4)], ...
    'PaperSize', [screenposition(3:4)]);

% Export to high-quality vector PDF
print(fig1, [fullfile(analysis_dir,'Fig_Boundary_Absolute_Saturation.pdf')], '-dpdf', '-painters');

%% KEY FIGURE 2: Combined Relative Saturation (NSD + HSD)
fig2 = figure(2);
clf;
% set(fig2, 'Position', [100, 100, 800, 600]);

% Prepare data
data_matrix_rel = [nsd_rel_upper', nsd_rel_lower', hsd_rel_upper', hsd_rel_lower'];

b = bar(data_matrix_rel, 'grouped');
set(gca, 'XTickLabel', {'±1,000', '±10,000', '±50,000'},'fontsize',14);
xlabel('Boundary condition (a.u.)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Percentage of saturated reactions (%)', 'FontSize', 14, 'FontWeight', 'bold');
% title(['FVA-sampling reactions with consistently high flux' newline '(|median flux| ≥99% of boundary)'], 'FontSize', 14, 'FontWeight', 'bold');
title(['(|median flux| ≥99% of boundary)'], 'FontSize', 18, 'FontWeight', 'bold');
legend({'NSD Upper', 'NSD Lower', 'HSD Upper', 'HSD Lower'}, ...
    'Location', 'northeast', 'FontSize', 8);
grid on;
ylim([0 1]);

% Color scheme
b(1).FaceColor = [0.2 0.4 0.8];
b(2).FaceColor = [0.4 0.6 1.0];
b(3).FaceColor = [0.8 0.4 0.2];
b(4).FaceColor = [1.0 0.6 0.4];

% Save
saveas(fig2, fullfile(analysis_dir, 'Fig_Boundary_Relative_Saturation.png'));
% saveas(fig2, fullfile(analysis_dir, 'Fig_Boundary_Relative_Saturation.pdf'));
savefig(fig2, fullfile(analysis_dir, 'Fig_Boundary_Relative_Saturation.fig'));
% Set up high-resolution vector PDF export
set(fig2, 'Units', 'inches');
screenposition = get(fig2, 'Position');
set(fig2, ...
    'PaperPosition', [0 0 screenposition(3:4)], ...
    'PaperSize', [screenposition(3:4)]);

% Export to high-quality vector PDF
print(fig2, [fullfile(analysis_dir,'Fig_Boundary_Relative_Saturation.pdf')], '-dpdf', '-painters');


%% Export summary CSV
fprintf('Exporting summary CSV...\n');

csv_data = [];
row_labels = {};

for bound_idx = 1:length(bound_conditions)
    current_bound = bound_conditions(bound_idx);
    bound_field = sprintf('bound_%d', current_bound);

    % NSD row
    csv_data = [csv_data; current_bound, ...
        comparison_data.(bound_field).NSD.abs_upper_pct, ...
        comparison_data.(bound_field).NSD.abs_lower_pct, ...
        comparison_data.(bound_field).NSD.rel_upper_pct, ...
        comparison_data.(bound_field).NSD.rel_lower_pct];
    row_labels{end+1} = sprintf('NSD_±%d', current_bound);

    % HSD row
    csv_data = [csv_data; current_bound, ...
        comparison_data.(bound_field).HSD.abs_upper_pct, ...
        comparison_data.(bound_field).HSD.abs_lower_pct, ...
        comparison_data.(bound_field).HSD.rel_upper_pct, ...
        comparison_data.(bound_field).HSD.rel_lower_pct];
    row_labels{end+1} = sprintf('HSD_±%d', current_bound);
end

headers = {'Boundary_Condition', 'Abs_Upper_Pct', 'Abs_Lower_Pct', 'Rel_Upper_Pct', 'Rel_Lower_Pct'};
csv_table = array2table(csv_data, 'VariableNames', headers, 'RowNames', row_labels);
writetable(csv_table, fullfile(analysis_dir, 'Boundary_Summary.csv'), 'WriteRowNames', true);

%% Print summary statistics
fprintf('\n=== BOUNDARY ANALYSIS SUMMARY ===\n\n');

fprintf('ABSOLUTE SATURATION (within 1%% of boundary):\n');
fprintf('  NSD Upper: %.2f%% → %.2f%% → %.2f%%\n', nsd_abs_upper);
fprintf('  NSD Lower: %.2f%% → %.2f%% → %.2f%%\n', nsd_abs_lower);
fprintf('  HSD Upper: %.2f%% → %.2f%% → %.2f%%\n', hsd_abs_upper);
fprintf('  HSD Lower: %.2f%% → %.2f%% → %.2f%%\n\n', hsd_abs_lower);

fprintf('RELATIVE SATURATION (median ≥99%% of boundary):\n');
fprintf('  NSD Upper: %.2f%% → %.2f%% → %.2f%%\n', nsd_rel_upper);
fprintf('  NSD Lower: %.2f%% → %.2f%% → %.2f%%\n', nsd_rel_lower);
fprintf('  HSD Upper: %.2f%% → %.2f%% → %.2f%%\n', hsd_rel_upper);
fprintf('  HSD Lower: %.2f%% → %.2f%% → %.2f%%\n\n', hsd_rel_lower);


fprintf('Files saved to: %s\n', analysis_dir);
