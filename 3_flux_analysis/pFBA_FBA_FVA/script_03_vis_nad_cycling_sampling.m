%% script_03_vis_nad_cycling_sampling

clc;
clear;
close all;

%% Load data
current_path = pwd;
vis_dir = fullfile(current_path, '03_vis');
output_dir = fullfile(current_path, '03_vis_cycling_sampling');
sampling_dir = fullfile(current_path, 'models');

if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Load manually curated cycling reactions
cycling_data = readtable(fullfile(vis_dir, 'cycling_reactions_detailed_analz.csv'));

% Load sampling results
disp('Loading sampling results...');
try
    nsd_sampling = load(fullfile(sampling_dir, 'nsd.mat'));
    hsd_sampling = load(fullfile(sampling_dir, 'hsd.mat'));
    disp('Sampling data loaded successfully.');
catch ME
    disp(['Error loading sampling data: ' ME.message]);
    disp(['Please check if files exist at: ' sampling_dir]);
    return;
end

sampling_data = {nsd_sampling, hsd_sampling};

fprintf('Loaded %d manually curated cycling reactions\n', height(cycling_data));

%% Separate by model
nsd_data = cycling_data(strcmp(cycling_data.Model, 'NSD'), :);
hsd_data = cycling_data(strcmp(cycling_data.Model, 'HSD'), :);

fprintf('NSD: %d reactions\n', height(nsd_data));
fprintf('HSD: %d reactions\n', height(hsd_data));

%% Add FVA-sampling flux to data
for model_idx = 1:2
    if model_idx == 1
        data = nsd_data;
        model_id = 'NSD';
    else
        data = hsd_data;
        model_id = 'HSD';
    end
    
    % Get sampling data for this model
    samp = sampling_data{model_idx}.samples;
    model_rxns = samp.rxns;
    
    % Use median of sampled fluxes (points field contains the sampled fluxes)
    sampling_flux = median(abs(samp.points), 2);
    
    % Match cycling reactions to sampling data
    fva_flux_values = zeros(height(data), 1);
    for i = 1:height(data)
        rxn_id = data.Reaction_ID{i};
        rxn_idx = find(strcmp(model_rxns, rxn_id));
        if ~isempty(rxn_idx)
            fva_flux_values(i) = sampling_flux(rxn_idx);
        end
    end
    
    % Add to data table
    if model_idx == 1
        nsd_data.FVA_Sampling_Flux = fva_flux_values;
    else
        hsd_data.FVA_Sampling_Flux = fva_flux_values;
    end
end

%% Figure: NSD and HSD Cycling Reactions - Now with FVA
fig1 = figure(1);
clf;
set(fig1, 'Position', [100, 100, 1400, 500]);

for model_idx = 1:2
    if model_idx == 1
        data = nsd_data;
        model_id = 'NSD';
    else
        data = hsd_data;
        model_id = 'HSD';
    end
    
    subplot(1, 2, model_idx);
    
    if height(data) > 0
        % Get top reactions (sort by flux reduction)
        [~, sorted_idx] = sort(data.Flux_Reduction, 'descend');
        n_rxns = min(15, height(data)); % Show top 15 or all if less
        top_idx = sorted_idx(1:n_rxns);
        
        % Create labels: gene (reaction_id)
        rxn_labels = cell(n_rxns, 1);
        for i = 1:n_rxns
            gene = data.Genes{top_idx(i)};
            rxn_id = data.Reaction_ID{top_idx(i)};
            
            % Handle empty or multiple genes
            if isempty(gene) || strcmp(gene, '')
                gene_str = 'No gene';
            else
                % If multiple genes separated by ';', take first one
                genes_split = strsplit(gene, ';');
                gene_str = strtrim(genes_split{1});
            end
            
            % Format: gene (reaction_id)
            rxn_labels{i} = sprintf('%s (%s)', gene_str, rxn_id);
        end
        
        % Apply abs() to all flux values
        fba_values = abs(data.FBA_Flux(top_idx));
        pfba_values = abs(data.pFBA_Flux(top_idx));
        fva_values = abs(data.FVA_Sampling_Flux(top_idx));
        
        % Create horizontal bar chart with three methods
        y_pos = 1:n_rxns;
        barh(y_pos - 0.25, fba_values, 0.25, 'DisplayName', 'FBA');
        hold on;
        barh(y_pos, pfba_values, 0.25, 'DisplayName', 'pFBA');
        barh(y_pos + 0.25, fva_values, 0.25, 'DisplayName', 'FVA-Sampling');
        
        set(gca, 'YTick', y_pos, 'YTickLabel', rxn_labels,'fontsize',14);
        xlabel('Flux (a.u.)', 'FontSize', 16, 'FontWeight', 'bold');
        title(sprintf('%s: NAD(P)-dependent cycling Reactions', model_id), 'FontSize', 18, 'FontWeight', 'bold');
        legend('Location', 'southeast', 'FontSize', 11);
        grid on;
        
        % Flip y-axis to show highest at top
        set(gca, 'YDir', 'reverse');
    else
        text(0.5, 0.5, 'No cycling reactions', ...
            'HorizontalAlignment', 'center', 'FontSize', 14);
        title(sprintf('%s: Cycling Reactions', model_id), 'FontSize', 14, 'FontWeight', 'bold');
    end
end

% Save figure
savefig(fig1, fullfile(output_dir, '01_Cycling_Reactions_Comparison.fig'));
saveas(fig1, fullfile(output_dir, '01_Cycling_Reactions_Comparison.png'), 'png');
set(fig1, 'Units', 'inches');
screenposition = get(fig1, 'Position');
set(fig1, 'PaperPosition', [0 0 screenposition(3:4)], ...
         'PaperSize', [screenposition(3:4)], ...
         'PaperUnits', 'inches', ...
         'PaperOrientation', 'portrait');
print(fig1, fullfile(output_dir, '01_Cycling_Reactions_Comparison.pdf'), '-dpdf', '-painters');

%% Figure 2: Flux Reduction Summary 


% Calculate reductions for pFBA
nsd_total_reduction = sum(nsd_data.Flux_Reduction);
hsd_total_reduction = sum(hsd_data.Flux_Reduction);
nsd_avg_reduction = mean(nsd_data.Percent_Reduction);
hsd_avg_reduction = mean(hsd_data.Percent_Reduction);


%% Export results to Excel
excel_file = fullfile(output_dir, 'cycling_reactions_with_sampling.xlsx');

% Combine NSD and HSD data
combined_data = [nsd_data; hsd_data];

% Sort by model and flux reduction
combined_data = sortrows(combined_data, {'Model', 'Flux_Reduction'}, {'ascend', 'descend'});

% Write to Excel
writetable(combined_data, excel_file, 'Sheet', 'All_Reactions');

% Write NSD and HSD separately
writetable(nsd_data, excel_file, 'Sheet', 'NSD');
writetable(hsd_data, excel_file, 'Sheet', 'HSD');

% Create summary table
summary_table = table();
summary_table.Model = {'NSD'; 'HSD'};
summary_table.Num_Reactions = [height(nsd_data); height(hsd_data)];
summary_table.pFBA_Total_Reduction = [nsd_total_reduction; hsd_total_reduction];
summary_table.pFBA_Avg_Percent_Reduction = [nsd_avg_reduction; hsd_avg_reduction];

writetable(summary_table, excel_file, 'Sheet', 'Summary');

fprintf('Results exported to: %s\n', excel_file);