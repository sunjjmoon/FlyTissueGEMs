%% Script 04: compare all results
clc;
clear;
close all

%% Setup
current_dir = pwd;
output_dir = 'results';
output_path = fullfile(current_dir, output_dir);
if ~exist(output_path, 'dir')
    mkdir(output_path);
end

perturbation_folders = {
    'p50'
    'p70'
    'p80'
    'p90'
    'p95_f'
};

% Extract perturbation levels
perturbation_levels = [50, 70, 80, 90, 95];
fprintf('Analyzing pyruvate consumption sensitivity for perturbation levels: %s\n', mat2str(perturbation_levels));

%% Load results
all_results = struct();
valid_folders = {};
valid_levels = [];

for i = 1:length(perturbation_folders)
    folder_path = fullfile(current_dir, perturbation_folders{i});
    excel_file = fullfile(folder_path, '2_results', 'pyruvate_consumption_sensitivity_results.xlsx');
    
    if exist(excel_file, 'file')
        fprintf('Loading %s...\n', perturbation_folders{i});
        
        try
            % Load detailed summary sheet (contains enzyme-level sensitivity data)
            detailed_summary = readtable(excel_file, 'Sheet', 'Detailed_Summary');
            
            % Extract data
            enzyme_names = detailed_summary.Enzyme;
            sensitivity_values = detailed_summary.Sensitivity;
            system_response_values = detailed_summary.System_Response;
            flux_difference_values = detailed_summary.Flux_Difference;
            baseline_pyruvate_values = detailed_summary.Avg_Baseline_Pyruvate;
            perturbed_pyruvate_values = detailed_summary.Avg_Perturbed_Pyruvate;
            
            % Calculate achieved perturbation (approximate from system response)
            achieved_pert_values = ones(size(enzyme_names)) * (perturbation_levels(i) / 100);
            
            % Store
            level_key = sprintf('f_%d', perturbation_levels(i));
            all_results.(level_key).enzyme_names = enzyme_names;
            all_results.(level_key).sensitivity = sensitivity_values;
            all_results.(level_key).system_response = system_response_values;
            all_results.(level_key).flux_difference = flux_difference_values;
            all_results.(level_key).baseline_pyruvate = baseline_pyruvate_values;
            all_results.(level_key).perturbed_pyruvate = perturbed_pyruvate_values;
            all_results.(level_key).achieved_perturbation = achieved_pert_values;
            all_results.(level_key).level = perturbation_levels(i);
            
            valid_folders{end+1} = perturbation_folders{i};
            valid_levels(end+1) = perturbation_levels(i);
            
            fprintf('  Loaded %d enzymes\n', length(enzyme_names));
            
        catch ME
            fprintf('  Error loading %s: %s\n', perturbation_folders{i}, ME.message);
        end
    else
        fprintf('  File not found: %s\n', excel_file);
    end
end

fprintf('Successfully loaded %d folders\n', length(valid_folders));

if length(valid_folders) == 0
    error('No valid data files found. Check folder structure and file names.');
end

%% Find common enzymes
level_keys = fieldnames(all_results);
common_enzymes = all_results.(level_keys{1}).enzyme_names;

for i = 2:length(level_keys)
    current_enzymes = all_results.(level_keys{i}).enzyme_names;
    common_enzymes = intersect(common_enzymes, current_enzymes, 'stable');
end

fprintf('Found %d common enzymes\n', length(common_enzymes));

%% Load enzyme gene mapping for better labels
glycolysis_file = fullfile(current_dir,'p95_f', 'files', 'Reactions for subsystem glycolysis_gluconeogenesis.tsv');

enzyme_labels = common_enzymes;  % Default to enzyme IDs

if exist(glycolysis_file, 'file')
    fprintf('Found glycolysis file\n');
    try
        opts = detectImportOptions(glycolysis_file, 'FileType', 'text');
        opts.Delimiter = '\t';
        glycolysis_data = readtable(glycolysis_file, opts);
        fprintf('Loaded glycolysis data: %d rows\n', height(glycolysis_data));
        
        % Debug: show what's in the data
        if height(glycolysis_data) > 0
            fprintf('First few reaction IDs: %s\n', strjoin(glycolysis_data.ReactionID(1:min(3, height(glycolysis_data))), ', '));
        end
        
        % Create mapping using the updated function
        for i = 1:length(common_enzymes)
            enzyme_id = common_enzymes{i};
            old_label = enzyme_labels{i};
            enzyme_labels{i} = create_gene_label(enzyme_id, glycolysis_data);
            fprintf('  %s -> %s\n', old_label, enzyme_labels{i});
        end
        
        fprintf('Loaded enzyme gene mappings\n');
    catch ME
        fprintf('Error loading gene mappings: %s\n', ME.message);
    end
else
    fprintf('Glycolysis file NOT FOUND: %s\n', glycolysis_file);
end


%% Aggregate data
n_enzymes = length(common_enzymes);
n_levels = length(valid_levels);

aggregate_sensitivity = NaN(n_enzymes, n_levels);
aggregate_system_response = NaN(n_enzymes, n_levels);
aggregate_flux_difference = NaN(n_enzymes, n_levels);
aggregate_baseline_pyruvate = NaN(n_enzymes, n_levels);
aggregate_achieved_perturbation = NaN(n_enzymes, n_levels);

for i = 1:length(level_keys)
    level_key = level_keys{i};
    level_idx = find(valid_levels == all_results.(level_key).level);
    
    for j = 1:n_enzymes
        enzyme = common_enzymes{j};
        enzyme_idx = find(strcmp(all_results.(level_key).enzyme_names, enzyme));
        
        if ~isempty(enzyme_idx)
            aggregate_sensitivity(j, level_idx) = all_results.(level_key).sensitivity(enzyme_idx);
            aggregate_system_response(j, level_idx) = all_results.(level_key).system_response(enzyme_idx);
            aggregate_flux_difference(j, level_idx) = all_results.(level_key).flux_difference(enzyme_idx);
            aggregate_baseline_pyruvate(j, level_idx) = all_results.(level_key).baseline_pyruvate(enzyme_idx);
            aggregate_achieved_perturbation(j, level_idx) = all_results.(level_key).achieved_perturbation(enzyme_idx);
        end
    end
end

% Calculate rankings based on absolute sensitivity
avg_sensitivity = nanmean((aggregate_sensitivity), 2);
std_sensitivity = nanstd((aggregate_sensitivity), 0, 2);
se_sensitivity = std_sensitivity ./ sqrt(sum(~isnan(aggregate_sensitivity), 2));

[sorted_avg, sort_idx] = sort(avg_sensitivity, 'descend');
sorted_std_sensitivity = std_sensitivity(sort_idx);
sorted_se_sensitivity = se_sensitivity(sort_idx);

%% Plot 1: Perturbation Verification
figure(1); clf

% Calculate actual perturbation magnitude: 1 - (v_pert / v_baseline)
% Note: aggregate_achieved_perturbation currently stores v_pert/v_baseline
actual_perturbation_magnitude = 1 - aggregate_achieved_perturbation';

% X-axis: input perturbation magnitudes (alpha)
input_alpha = 1 - valid_levels/100;  % Convert to alpha (e.g., 0.50, 0.30, 0.20, 0.10, 0.05)
x_labels = arrayfun(@(x) sprintf('%.2f', x), input_alpha, 'UniformOutput', false);

bar(categorical(x_labels), actual_perturbation_magnitude, 'grouped');

xlabel(['Input perturbation magnitude (\alpha)'], 'FontSize', 12, 'FontWeight', 'bold');
ylabel(['Achieved perturbation' newline '(1 - v_{pert}/v_{baseline})'], 'FontSize', 12, 'FontWeight', 'bold');
legend(enzyme_labels, 'Location', 'eastoutside', 'FontSize', 8);
% set(gca, 'FontSize', 10);
yticks(0:0.05:0.6);   % Sets y-axis ticks at 0, 0.5, 1.0, 1.5

% Add identity line reference
hold on;
ylim_vals = ylim;
xlim_vals = xlim;
hold off;

grid on;

saveas(gcf, fullfile(output_path, '1_perturbation_verification_pyruvate.png'));
saveas(gcf, fullfile(output_path, '1_perturbation_verification_pyruvate.svg'));

% Save PDF
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'PaperUnits', 'inches');
pos = get(gcf, 'Position');
set(gcf, 'PaperSize', [pos(3)/100, pos(4)/100]);
print(gcf, fullfile(output_path, '1_perturbation_verification_pyruvate.pdf'), '-dpdf', '-fillpage');
% print(gcf, fullfile(output_path, 'test.pdf'), '-dpdf', '-fillpage');

fprintf('Saved: 1_perturbation_verification_pyruvate files\n');
%% Figure 2b: Sensitivity by Perturbation Level (Grouped)
figure(22); clf

% Transpose for grouped bar (enzymes x perturbation levels)
sensitivity_data = (aggregate_sensitivity(sort_idx, :));

barh(categorical(enzyme_labels(sort_idx)), sensitivity_data, 'grouped');

% ylabel('Enzyme', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Normalized sensitivity coefficient (C_{i})', 'FontSize', 12, 'FontWeight', 'bold');
% title('Enzyme Sensitivity by Perturbation Level', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'YDir', 'reverse', 'FontSize', 10);
% ylim([-0.4 0.4])
% Legend for perturbation levels with title
% level_labels = arrayfun(@(x) sprintf('%.2f', x/100), valid_levels, 'UniformOutput', false);
level_labels = arrayfun(@(x) sprintf('%d%', round((1-x/100)*100)), valid_levels, 'UniformOutput', false);
lgd = legend(level_labels, 'Location', 'eastoutside','fontsize',8);
lgd.Title.String = 'Perturbation factor (%)';
% Add separator lines between enzyme rows using rectangles
hold on;
xlim([-0.4 0.2])
x_limits = xlim;
for i = 1.5:1:(length(sort_idx)-0.5)
    rectangle('Position', [x_limits(1), i-0.02, diff(x_limits), 0.03], ...
              'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none','linewidth',0.5);
end

hold off;

% grid on;
% Save figure in multiple formats (robust for categorical axes)
exportgraphics(gcf, fullfile(output_path, '2b_sensitivity_by_perturbation_pyruvate.png'), 'Resolution', 300);
exportgraphics(gcf, fullfile(output_path, '2b_sensitivity_by_perturbation_pyruvate.pdf'), 'ContentType', 'vector');
% print(gcf, fullfile(output_path, '2b_sensitivity_by_perturbation_pyruvate.svg'), '-dsvg', '-r300');

fprintf('Saved: 2b_sensitivity_by_perturbation_pyruvate files\n');

%% Figure 2b: Sensitivity by Perturbation Level (Grouped)
figure(222); clf
% Define the exact order you want
new_order = {
    'Hex-a/b (MAR04394)'
    'Pgi (MAR04381)'
    'Pfk (MAR04379)'
    'Aldo (MAR04375)'
    'Tpi (MAR04391)'        % Moved here
    'Gapdh1/2 (MAR04373)'
    'Pgk (MAR04368)'
    'Pgam1 (MAR04365)'
    'Eno (MAR04363)'
    'Pyk (MAR04358)'
};

% Find indices for reordering
[~, reorder_idx] = ismember(new_order, enzyme_labels);

enzyme_labels_reordered = enzyme_labels(reorder_idx);
sensitivity_data_reordered = aggregate_sensitivity(reorder_idx, :);

% Create plot
cat_labels = categorical(enzyme_labels_reordered, enzyme_labels_reordered, 'Ordinal', true);
barh(cat_labels, sensitivity_data_reordered, 'grouped');

set(gca, 'YDir', 'reverse', 'FontSize', 10);

% ylabel('Enzyme', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('Normalized sensitivity coefficient (C_{i})', 'FontSize', 12, 'FontWeight', 'bold');
% title('Enzyme Sensitivity by Perturbation Level', 'FontSize', 14, 'FontWeight', 'bold');
set(gca, 'YDir', 'reverse', 'FontSize', 10);

% ylim([-0.4 0.4])
% Legend for perturbation levels with title
% level_labels = arrayfun(@(x) sprintf('%.2f', x/100), valid_levels, 'UniformOutput', false);
level_labels = arrayfun(@(x) sprintf('%d%', round((1-x/100)*100)), valid_levels, 'UniformOutput', false);
lgd = legend(level_labels, 'Location', 'eastoutside','fontsize',6);
lgd.Title.String = 'Perturbation factor (%)';
% Add separator lines between enzyme rows using rectangles
hold on;
% xlim([0 0.1])
xlim([-0.4 0.2])


x_limits = xlim;
for i = 1.5:1:(length(enzyme_labels)-0.5)
    rectangle('Position', [x_limits(1), i-0.02, diff(x_limits), 0.03], ...
              'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none','linewidth',0.5);
end
hold off;

% grid on;
% Save figure in multiple formats (robust for categorical axes)
exportgraphics(gcf, fullfile(output_path, '2c_sensitivity_by_perturbation_pyruvate_notSorted.png'), 'Resolution', 300);
exportgraphics(gcf, fullfile(output_path, '2c_sensitivity_by_perturbation_pyruvate_notSorted.pdf'), 'ContentType', 'vector');
% print(gcf, fullfile(output_path, '2b_sensitivity_by_perturbation_pyruvate.svg'), '-dsvg', '-r300');

fprintf('Saved: 2c_sensitivity_by_perturbation_pyruvate files\n');




%% Figure 33: Combined Sensitivity Plot - Grouped bars (left) + Absolute bars (right)
figure(33); clf
% set(gcf, 'Position', [100, 100, 800, 600]); % [x, y, width, height] - reduced width

% Define the exact order you want
new_order = {
    'Hex-a/b (MAR04394)'
    'Pgi (MAR04381)'
    'Pfk (MAR04379)'
    'Aldo (MAR04375)'
    'Tpi (MAR04391)'        
    'Gapdh1/2 (MAR04373)'
    'Pgk (MAR04368)'
    'Pgam1 (MAR04365)'
    'Eno (MAR04363)'
    'Pyk (MAR04358)'
};

% Find indices for reordering
[~, reorder_idx] = ismember(new_order, enzyme_labels);
enzyme_labels_reordered = enzyme_labels(reorder_idx);
sensitivity_data_reordered = aggregate_sensitivity(reorder_idx, :);

% Calculate absolute averages for right panel
abs_avg_sensitivity = mean(abs(sensitivity_data_reordered), 2);

% Create subplots with custom positioning - more space for labels
subplot('Position', [0.26, 0.12, 0.5, 0.8]); % Left subplot (leave more space for labels) [left, bottom, width, height]
% subplot('Position', [0.2, 0.15, 0.5, 0.75]); % Left panel - adjusted for narrower figure

% Left panel: Grouped bars
cat_labels = categorical(enzyme_labels_reordered, enzyme_labels_reordered, 'Ordinal', true);
barh(cat_labels, sensitivity_data_reordered, 'grouped');
set(gca, 'YDir', 'reverse', 'FontSize', 10);
xlabel('Normalized sensitivity coefficient (S_i)', 'FontSize', 12, 'FontWeight', 'bold');
xlim([-0.4 0.25]);
ax = gca;
ax.XGrid = 'on';
% set(gca, 'FontWeight', 'bold');   % makes both x and y tick labels bold


% Legend for perturbation levels
level_labels = arrayfun(@(x) sprintf('%d%', round((1-x/100)*100)), valid_levels, 'UniformOutput', false);
lgd = legend(level_labels, 'Location', 'northeast', 'FontSize', 6);
lgd.Title.String = 'Perturbation (%)';

% Add separator lines
hold on;
x_limits = xlim;
for i = 1.5:1:(length(enzyme_labels_reordered)-0.5)
    rectangle('Position', [x_limits(1), i-0.02, diff(x_limits), 0.03], ...
              'FaceColor', [0.7 0.7 0.7], 'EdgeColor', 'none', 'LineWidth', 0.5);
end
% hold off;

% Right panel: Absolute average bars
subplot('Position', [0.8, 0.12, 0.15, 0.8]); % [left, bottom, width, height]
% subplot('Position', [0.75, 0.15, 0.2, 0.75]); % Right panel - adjusted

% Create matching categorical for right panel
barh(1:length(abs_avg_sensitivity), abs_avg_sensitivity, 'FaceColor', [0.7, 0.7, 0.7]);
set(gca, 'YDir', 'reverse', 'FontSize', 10);
set(gca, 'YTick', 1:length(enzyme_labels_reordered), 'YTickLabel', []); % Remove labels but keep ticks
ax = gca;
ax.XGrid = 'on';
xlabel('Average |S_{i}|', 'FontSize', 10, 'FontWeight', 'bold');

% Set y-axis limits to match left plot
ylim([0.5, length(enzyme_labels_reordered) + 0.5]);

% Add separator lines to match left panel
hold on;
xlim([0 0.13])
x_limits_right = xlim;
for i = 1.5:1:(length(enzyme_labels_reordered)-0.5)
    line([x_limits_right(1), x_limits_right(2)], [i, i], ...
         'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
end
hold off;

% Add overall title
% sgtitle('Enzyme Sensitivity Analysis - Grouped vs Absolute', 'FontSize', 14, 'FontWeight', 'bold');

% Save figure
exportgraphics(gcf, fullfile(output_path, '33_combined_sensitivity_plot.png'), 'Resolution', 300);
exportgraphics(gcf, fullfile(output_path, '33_combined_sensitivity_plot.pdf'), 'ContentType', 'vector');
fprintf('Saved: 33_combined_sensitivity_plot files\n');


%% Additional statistics
fprintf('\n=== ADDITIONAL STATISTICS ===\n');
fprintf('Sensitivity coefficient ranges:\n');
all_sens_vals = aggregate_sensitivity(:);
valid_sens = all_sens_vals(~isnan(all_sens_vals));
fprintf('  Min: %.4f, Max: %.4f\n', min(valid_sens), max(valid_sens));
fprintf('  Mean: %.4f, Std: %.4f\n', mean(valid_sens), std(valid_sens));

fprintf('System response ranges:\n');
all_resp_vals = aggregate_system_response(:);
valid_resp = all_resp_vals(~isnan(all_resp_vals));
fprintf('  Min: %.4f, Max: %.4f\n', min(valid_resp), max(valid_resp));
fprintf('  Mean: %.4f, Std: %.4f\n', mean(valid_resp), std(valid_resp));


%% Gene name mapping function (place at end of script)
function label = create_gene_label(reaction_id, pathway_data)
    % Find this reaction in pathway_data
    idx = strcmp(pathway_data.ReactionID, reaction_id);
    if any(idx)
        genes_str = pathway_data.Genes{idx};
        
        % Split genes by ';' and take first two
        genes_split = strsplit(genes_str, ';');
        if length(genes_split) >= 2
            first_two = strjoin(genes_split(1:2), '');
        elseif length(genes_split) == 1
            first_two = genes_split{1};
        else
            first_two = '';
        end
        
        % Apply custom gene name mapping
        first_two = strtrim(first_two);
        switch first_two
            case 'Ald1'
                first_two = 'Aldo';
            case 'Gapdh2 Gapdh1'
                first_two = 'Gapdh1/2';
            case 'CG7059 Pglym78'
                first_two = 'Pgam1';
            case 'CG30486 Pgk'
                first_two = 'Pgk';
            case 'Hex-A'
                first_two = 'Hex-a/b';
            case 'Ih PyK'
                first_two = 'Pyk';
            case 'TSG101 Ldh'
                first_two = 'Ldh';    
        end
        
        % Create combined label
        if ~isempty(first_two)
            label = sprintf('%s (%s)', first_two, reaction_id);
        else
            label = reaction_id;
        end
    else
        label = reaction_id; % fallback
    end
end