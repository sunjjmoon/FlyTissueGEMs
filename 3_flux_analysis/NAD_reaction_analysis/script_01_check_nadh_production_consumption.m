%% script_1_check_nadh_production_consumption
% NADH Production/Consumption Analysis with Directionality 
% Goal: Analyze NADH production vs consumption

clc
close all

pathway = pwd;
[parentFolder, ~, ~] = fileparts(pathway);

save_dirFlux = '01_nad_balance';

subfolder = [pathway '\' save_dirFlux];
if ~exist(subfolder, 'dir')
    mkdir(subfolder)
end

model_sampling_in_name = {'NSD','HSD'};

%% Load data (if not already loaded)
if ~exist('model_sampling_in', 'var') || isempty(model_sampling_in)
    model_sampling_in = [];
    for i = 1:length(model_sampling_in_name)
       fileName = fullfile(pathway,'C_sampling', model_sampling_in_name{i});
       fileName = strcat(fileName,'.mat');
       model_sampling_in{i} = load(fileName);

    end
end

fprintf('=== CORRECTED NADH DIRECTIONALITY ANALYSIS ===\n');

%% Step 1: Find all NADH-containing reactions
fprintf('=== Step 1: Finding NADH-containing reactions ===\n');

model = model_sampling_in{1}.samples;
compartments = model.comps;

% Find all NADH metabolites (MAM02553) across compartments
nadh_metabolites = {};
nadh_compartments = {};

for i = 1:length(model.mets)
    met_id = model.mets{i};
    if contains(met_id, 'MAM02553')
        nadh_metabolites{end+1} = met_id;
        % Extract compartment
        comp_match = regexp(met_id, '\[([a-z])\]', 'tokens');
        if ~isempty(comp_match)
            nadh_compartments{end+1} = comp_match{1}{1};
        end
    end
end

fprintf('NADH metabolites found: %d\n', length(nadh_metabolites));
disp('NADH metabolites:'); 
disp(nadh_metabolites');

% Find all reactions involving NADH
nadh_reactions = [];
for i = 1:length(nadh_metabolites)
    met_idx = find(strcmp(model.mets, nadh_metabolites{i}));
    if ~isempty(met_idx)
        rxn_indices = find(model.S(met_idx, :) ~= 0);
        nadh_reactions = [nadh_reactions, rxn_indices];
    end
end

nadh_reactions = unique(nadh_reactions);


fprintf('Total NADH reactions (excluding demand): %d\n', length(nadh_reactions));

%% Step 2: Analyze NADH directionality for each reaction (using median fluxes)
fprintf('\n=== Step 2: Analyzing NADH directionality using median fluxes ===\n');

reaction_analysis = struct();
comp_names = {'c', 'cytosol'; 'm', 'mitochondria'; 'p', 'peroxisome'; 'r', 'ER'; 'e', 'extracellular'; 'l', 'lysosome'; 'g', 'golgi'; 'n', 'nucleus'; 'i', 'intermembrane'};

for i = 1:length(nadh_reactions)
    rxn_idx = nadh_reactions(i);
    rxn_id = model.rxns{rxn_idx};
    
    % Find NADH metabolites in this reaction and their stoichiometry
    nadh_info = struct();
    nadh_stoich = [];
    nadh_mets_in_rxn = {};
    nadh_comps_in_rxn = {};
    
    for j = 1:length(nadh_metabolites)
        met_idx = find(strcmp(model.mets, nadh_metabolites{j}));
        if ~isempty(met_idx)
            stoich_coeff = model.S(met_idx, rxn_idx);
            if stoich_coeff ~= 0
                nadh_stoich = [nadh_stoich, stoich_coeff];
                nadh_mets_in_rxn{end+1} = nadh_metabolites{j};
                
                % Extract compartment
                comp_match = regexp(nadh_metabolites{j}, '\[([a-z])\]', 'tokens');
                if ~isempty(comp_match)
                    nadh_comps_in_rxn{end+1} = comp_match{1}{1};
                end
            end
        end
    end
    
    % Store basic reaction information
    reaction_analysis(i).rxn_idx = rxn_idx;
    reaction_analysis(i).rxn_id = rxn_id;
    reaction_analysis(i).rxn_name = model.rxnNames{rxn_idx};
    reaction_analysis(i).nadh_metabolites = nadh_mets_in_rxn;
    reaction_analysis(i).nadh_stoichiometry = nadh_stoich;
    reaction_analysis(i).compartments = unique(nadh_comps_in_rxn);
    
    % Handle subsystem safely
    subsystem_val = model.subSystems{rxn_idx};
    if isempty(subsystem_val) || (iscell(subsystem_val) && isempty(subsystem_val{1}))
        reaction_analysis(i).subsystem = 'Unknown';
    elseif iscell(subsystem_val)
        reaction_analysis(i).subsystem = subsystem_val{1};
    else
        reaction_analysis(i).subsystem = char(subsystem_val);
    end
end

%% Step 4: Calculate condition-specific effects
fprintf('\n=== Step 4: Calculating condition-specific effects ===\n');
model_names = {'NSD','HSD'};
conditions = model_names;
for cond = 1:length(conditions)
    condition_name = conditions{cond};
    points = model_sampling_in{cond}.samples.points;
    
    fprintf('Processing %s condition...\n', condition_name);
    
    % Calculate median flux for each NADH reaction
    for i = 1:length(nadh_reactions)
        rxn_idx = nadh_reactions(i);
        flux_values = points(rxn_idx, :);
        median_flux = mean(flux_values); % median
        std_flux = std(flux_values);
        
        % Store flux information
        reaction_analysis(i).([condition_name '_median_flux']) = median_flux;
        reaction_analysis(i).([condition_name '_std_flux']) = std_flux;
        reaction_analysis(i).([condition_name '_flux_range']) = [min(flux_values), max(flux_values)];
        
        % Calculate net NADH change for this reaction
        net_nadh_change = 0;
        nadh_effects = {};
        
        for j = 1:length(reaction_analysis(i).nadh_stoichiometry)
            stoich = reaction_analysis(i).nadh_stoichiometry(j);
            nadh_met = reaction_analysis(i).nadh_metabolites{j};
            
            % Calculate net NADH change: stoichiometry * flux (ensure full matrix)
            nadh_change = full(stoich * median_flux);
            net_nadh_change = net_nadh_change + nadh_change;
            
            flux_threshold = 0;
            
            if abs(median_flux) > flux_threshold
                if nadh_change > flux_threshold
                    nadh_effects{end+1} = sprintf('%s: PRODUCTION (%.6f)', nadh_met, nadh_change);
                elseif nadh_change < -flux_threshold
                    nadh_effects{end+1} = sprintf('%s: CONSUMPTION (%.6f)', nadh_met, nadh_change);
                else
                    nadh_effects{end+1} = sprintf('%s: NEUTRAL (%.6f)', nadh_met, nadh_change);
                end
            else
                nadh_effects{end+1} = sprintf('%s: INACTIVE', nadh_met);
            end
        end
        
        reaction_analysis(i).([condition_name '_nadh_effects']) = nadh_effects;
        reaction_analysis(i).([condition_name '_net_nadh_change']) = net_nadh_change;
        
        % Classify overall reaction effect
        if abs(median_flux) > flux_threshold
            if net_nadh_change > flux_threshold
                reaction_analysis(i).([condition_name '_classification']) = 'NADH_PRODUCING';
            elseif net_nadh_change < -flux_threshold
                reaction_analysis(i).([condition_name '_classification']) = 'NADH_CONSUMING';
            else
                reaction_analysis(i).([condition_name '_classification']) = 'NADH_NEUTRAL';
            end
        else
            reaction_analysis(i).([condition_name '_classification']) = 'INACTIVE';
        end
    end
end

%% Step 5: Compare conditions
fprintf('\n=== Step 5: Comparing conditions ===\n');

for i = 1:length(nadh_reactions)
    nsd_production = reaction_analysis(i).NSD_net_nadh_change;
    hsd_production = reaction_analysis(i).HSD_net_nadh_change;
    
    reaction_analysis(i).production_difference = hsd_production - nsd_production;
    
    % Handle division by zero for fold change
    if abs(nsd_production) < flux_threshold
        if abs(hsd_production) < flux_threshold
            reaction_analysis(i).production_fold_change = 1;
        else
            reaction_analysis(i).production_fold_change = Inf;
        end
    else
        reaction_analysis(i).production_fold_change = hsd_production / nsd_production;
    end
end

%% Step 6: Create comprehensive tables
fprintf('\n=== Step 6: Creating results tables ===\n');

% Pre-allocate for table creation
n_reactions = length(nadh_reactions);
results_table = table();

% Basic reaction info
results_table.Reaction_ID = cell(n_reactions, 1);
results_table.Reaction_Name = cell(n_reactions, 1);
results_table.Subsystem = cell(n_reactions, 1);
results_table.Compartments = cell(n_reactions, 1);
results_table.NADH_Metabolites = cell(n_reactions, 1);
results_table.NADH_Stoichiometry = cell(n_reactions, 1);
results_table.Gene_Rules = cell(n_reactions, 1);

% Flux data
results_table.NSD_Median_Flux = zeros(n_reactions, 1);
results_table.HSD_Median_Flux = zeros(n_reactions, 1);
results_table.NSD_Net_NADH_Change = zeros(n_reactions, 1);
results_table.HSD_Net_NADH_Change = zeros(n_reactions, 1);
results_table.Production_Difference = zeros(n_reactions, 1);
results_table.Production_Fold_Change = zeros(n_reactions, 1);

% Classifications
results_table.NSD_Classification = cell(n_reactions, 1);
results_table.HSD_Classification = cell(n_reactions, 1);
results_table.NSD_Effects = cell(n_reactions, 1);
results_table.HSD_Effects = cell(n_reactions, 1);

% Fill table
for i = 1:n_reactions
    results_table.Reaction_ID{i} = reaction_analysis(i).rxn_id;
    results_table.Reaction_Name{i} = reaction_analysis(i).rxn_name;
    results_table.Subsystem{i} = reaction_analysis(i).subsystem;
    results_table.Compartments{i} = strjoin(reaction_analysis(i).compartments, '; ');
    results_table.NADH_Metabolites{i} = strjoin(reaction_analysis(i).nadh_metabolites, '; ');
    results_table.NADH_Stoichiometry{i} = num2str(reaction_analysis(i).nadh_stoichiometry);
    results_table.Gene_Rules{i} = model.grRules{reaction_analysis(i).rxn_idx};

    results_table.NSD_Median_Flux(i) = reaction_analysis(i).NSD_median_flux;
    results_table.HSD_Median_Flux(i) = reaction_analysis(i).HSD_median_flux;
    results_table.NSD_Net_NADH_Change(i) = reaction_analysis(i).NSD_net_nadh_change;
    results_table.HSD_Net_NADH_Change(i) = reaction_analysis(i).HSD_net_nadh_change;
    results_table.Production_Difference(i) = reaction_analysis(i).production_difference;
    results_table.Production_Fold_Change(i) = reaction_analysis(i).production_fold_change;
    
    results_table.NSD_Classification{i} = reaction_analysis(i).NSD_classification;
    results_table.HSD_Classification{i} = reaction_analysis(i).HSD_classification;
    results_table.NSD_Effects{i} = strjoin(reaction_analysis(i).NSD_nadh_effects, '; ');
    results_table.HSD_Effects{i} = strjoin(reaction_analysis(i).HSD_nadh_effects, '; ');
end

%% Step 7: Create compartment summary
fprintf('\n=== Step 7: Creating compartment summary ===\n');

% Get unique compartments
all_compartments = {};
for i = 1:n_reactions
    all_compartments = [all_compartments, reaction_analysis(i).compartments];
end
unique_compartments = unique(all_compartments);

% IMPORTANT: Track which reactions we've already counted to avoid double-counting
counted_reactions = [];

% Create compartment summary
comp_summary = table();
for i = 1:length(unique_compartments)
    comp = unique_compartments{i};
    comp_name_idx = find(strcmp(comp_names(:,1), comp));
    if ~isempty(comp_name_idx)
        comp_full_name = comp_names{comp_name_idx, 2};
    else
        comp_full_name = comp;
    end
    
    % Find reactions that PRIMARILY occur in this compartment
    % (to avoid double-counting multi-compartment reactions)
    comp_reactions = [];
    for j = 1:n_reactions
        % Check if this reaction involves NADH in this compartment
        involves_this_comp = false;
        for k = 1:length(reaction_analysis(j).nadh_metabolites)
            nadh_met = reaction_analysis(j).nadh_metabolites{k};
            if contains(nadh_met, sprintf('[%s]', comp))
                involves_this_comp = true;
                break;
            end
        end
        
        if involves_this_comp
            comp_reactions = [comp_reactions, j];
        end
    end
    
    % Count by classification for reactions in this compartment
    nsd_producing = 0;
    nsd_consuming = 0;
    nsd_neutral = 0;
    nsd_inactive = 0;
    hsd_producing = 0;
    hsd_consuming = 0;
    hsd_neutral = 0;
    hsd_inactive = 0;
    
    % Calculate NADH changes SPECIFIC to this compartment
    nsd_total_production = 0;
    nsd_total_consumption = 0;
    hsd_total_production = 0;
    hsd_total_consumption = 0;
    
    for j = comp_reactions
        % Get the NADH change ONLY for metabolites in THIS compartment
        nsd_nadh_change_this_comp = 0;
        hsd_nadh_change_this_comp = 0;
        
        for k = 1:length(reaction_analysis(j).nadh_metabolites)
            nadh_met = reaction_analysis(j).nadh_metabolites{k};
            
            % Check if this NADH metabolite is in the current compartment
            if contains(nadh_met, sprintf('[%s]', comp))
                stoich = reaction_analysis(j).nadh_stoichiometry(k);
                
                % Calculate compartment-specific NADH change
                nsd_flux = reaction_analysis(j).NSD_median_flux;
                hsd_flux = reaction_analysis(j).HSD_median_flux;
                
                nsd_change = stoich * nsd_flux;
                hsd_change = stoich * hsd_flux;
                
                nsd_nadh_change_this_comp = nsd_nadh_change_this_comp + nsd_change;
                hsd_nadh_change_this_comp = hsd_nadh_change_this_comp + hsd_change;
            end
        end
        
        % Update production/consumption for this compartment
        if nsd_nadh_change_this_comp > 0
            nsd_total_production = nsd_total_production + nsd_nadh_change_this_comp;
            nsd_producing = nsd_producing + 1;
        elseif nsd_nadh_change_this_comp < 0
            nsd_total_consumption = nsd_total_consumption + nsd_nadh_change_this_comp;
            nsd_consuming = nsd_consuming + 1;
        elseif abs(reaction_analysis(j).NSD_median_flux) > flux_threshold
            nsd_neutral = nsd_neutral + 1;
        else
            nsd_inactive = nsd_inactive + 1;
        end
        
        if hsd_nadh_change_this_comp > 0
            hsd_total_production = hsd_total_production + hsd_nadh_change_this_comp;
            hsd_producing = hsd_producing + 1;
        elseif hsd_nadh_change_this_comp < 0
            hsd_total_consumption = hsd_total_consumption + hsd_nadh_change_this_comp;
            hsd_consuming = hsd_consuming + 1;
        elseif abs(reaction_analysis(j).HSD_median_flux) > flux_threshold
            hsd_neutral = hsd_neutral + 1;
        else
            hsd_inactive = hsd_inactive + 1;
        end
    end
    
    % Store compartment data
    comp_summary_row = table({comp}, {comp_full_name}, length(comp_reactions), ...
        nsd_producing, nsd_consuming, nsd_neutral, nsd_inactive, ...
        hsd_producing, hsd_consuming, hsd_neutral, hsd_inactive, ...
        nsd_total_production, nsd_total_consumption, ...
        hsd_total_production, hsd_total_consumption, ...
        nsd_total_production + nsd_total_consumption, ...
        hsd_total_production + hsd_total_consumption);
    
    if i == 1
        comp_summary = comp_summary_row;
    else
        comp_summary = [comp_summary; comp_summary_row];
    end
end

comp_summary.Properties.VariableNames = {'Compartment', 'Full_Name', 'Total_Reactions', ...
    'NSD_Producing', 'NSD_Consuming', 'NSD_Neutral', 'NSD_Inactive', ...
    'HSD_Producing', 'HSD_Consuming', 'HSD_Neutral', 'HSD_Inactive', ...
    'NSD_Total_Production', 'NSD_Total_Consumption', ...
    'HSD_Total_Production', 'HSD_Total_Consumption', ...
    'NSD_Net_NADH_Change', 'HSD_Net_NADH_Change'};

% Add verification step
fprintf('\n=== Compartment Balance Verification ===\n');
total_nsd_net = sum(comp_summary.NSD_Net_NADH_Change);
total_hsd_net = sum(comp_summary.HSD_Net_NADH_Change);
fprintf('Sum of compartment NSD net changes: %.8f\n', total_nsd_net);
fprintf('Sum of compartment HSD net changes: %.8f\n', total_hsd_net);

% Compare with overall summary
if exist('overall_summary', 'var')
    overall_nsd_net = overall_summary.Net_Change(strcmp(overall_summary.Condition, 'NSD'));
    overall_hsd_net = overall_summary.Net_Change(strcmp(overall_summary.Condition, 'HSD'));
    fprintf('Overall summary NSD net change: %.8f\n', overall_nsd_net);
    fprintf('Overall summary HSD net change: %.8f\n', overall_hsd_net);
    
    fprintf('\nDiscrepancy check:\n');
    fprintf('NSD difference: %.8f\n', abs(total_nsd_net - overall_nsd_net));
    fprintf('HSD difference: %.8f\n', abs(total_hsd_net - overall_hsd_net));
    
    if abs(total_nsd_net - overall_nsd_net) > 1e-6 || abs(total_hsd_net - overall_hsd_net) > 1e-6
        fprintf('WARNING: Compartment totals do not match overall totals!\n');
        fprintf('This may indicate transport reactions between compartments.\n');
    end
end

%% Step 8: Create overall summary based on reaction-level analysis
fprintf('\n=== Step 8: Creating overall summary ===\n');

overall_summary = table();

% Calculate totals from reaction-level analysis
for cond = 1:length(conditions)
    condition_name = conditions{cond};
    
    % Sum up production and consumption from all reactions
    total_production = 0;
    total_consumption = 0;
    
    for i = 1:length(nadh_reactions)
        nadh_change = reaction_analysis(i).([condition_name '_net_nadh_change']);
        if nadh_change > 0
            total_production = total_production + nadh_change;
        else
            total_consumption = total_consumption + nadh_change;
        end
    end
    
    net_change = total_production + total_consumption;
    
    overall_summary_row = table({condition_name}, total_production, total_consumption, net_change);
    
    if cond == 1
        overall_summary = overall_summary_row;
    else
        overall_summary = [overall_summary; overall_summary_row];
    end
end

overall_summary.Properties.VariableNames = {'Condition', 'Total_Production', 'Total_Consumption', 'Net_Change'};

%% Step 9: Save all results
fprintf('\n=== Step 9: Saving results ===\n');

excel_file = fullfile(subfolder, 'NADH_Production_Consumption_Analysis.xlsx');
writetable(results_table, excel_file, 'Sheet', 'All_NADH_Reactions');
writetable(comp_summary, excel_file, 'Sheet', 'Compartment_Summary');
writetable(overall_summary, excel_file, 'Sheet', 'Overall_Summary');

fprintf('Results saved to: %s\n', excel_file);

%% Step 10: Print summary
fprintf('\n=== NADH DIRECTIONALITY ANALYSIS SUMMARY ===\n');
fprintf('Total NADH reactions analyzed: %d\n', n_reactions);

% Extract totals from overall summary
nsd_data = overall_summary(strcmp(overall_summary.Condition, 'NSD'), :);
hsd_data = overall_summary(strcmp(overall_summary.Condition, 'HSD'), :);

if ~isempty(nsd_data) && ~isempty(hsd_data)
    fprintf('\nNSD Condition:\n');
    fprintf('  Total NADH production: %.8f\n', nsd_data.Total_Production);
    fprintf('  Total NADH consumption: %.8f\n', nsd_data.Total_Consumption);
    fprintf('  Net NADH change: %.8f\n', nsd_data.Net_Change);
    
    fprintf('\nHSD Condition:\n');
    fprintf('  Total NADH production: %.8f\n', hsd_data.Total_Production);
    fprintf('  Total NADH consumption: %.8f\n', hsd_data.Total_Consumption);
    fprintf('  Net NADH change: %.8f\n', hsd_data.Net_Change);
    
    fprintf('\nComparison (HSD vs NSD):\n');
    if nsd_data.Total_Production ~= 0
        fprintf('  Production fold change: %.6f\n', hsd_data.Total_Production / nsd_data.Total_Production);
    end
    if nsd_data.Total_Consumption ~= 0
        fprintf('  Consumption fold change: %.6f\n', hsd_data.Total_Consumption / nsd_data.Total_Consumption);
    end
end

fprintf('\nDirectionality analysis complete!\n');



%% Step 11: Create grouped bar graphs for NADH production and consumption
fprintf('\n=== Step 11: Creating grouped bar graphs ===\n');

% Extract data from overall summary
nsd_production = abs(overall_summary.Total_Production(strcmp(overall_summary.Condition, 'NSD')));
nsd_consumption = abs(overall_summary.Total_Consumption(strcmp(overall_summary.Condition, 'NSD')));
hsd_production = abs(overall_summary.Total_Production(strcmp(overall_summary.Condition, 'HSD')));
hsd_consumption = abs(overall_summary.Total_Consumption(strcmp(overall_summary.Condition, 'HSD')));

% Prepare data for grouped bar plot
% Each row represents a condition (NSD, HSD)
% Each column represents Production and Consumption
bar_data = [nsd_production, nsd_consumption;
            hsd_production, hsd_consumption];

% Create figure
figure('Position', [100, 100, 800, 600]);

% Create grouped bar graph
b = bar(bar_data, 'grouped');

% Customize bar colors
b(1).FaceColor = [0.2, 0.6, 0.8];  % Blue for Production
b(2).FaceColor = [0.8, 0.3, 0.3];  % Red for Consumption

% Set x-axis labels
set(gca, 'XTickLabel', {'NSD', 'HSD'});
% xlabel('Condition', 'FontSize', 12, 'FontWeight', 'bold');

% Set y-axis label
ylabel('Flux (a.u.)', 'FontSize', 12, 'FontWeight', 'bold');

% Set title
title('NADH-dependent reactions', 'FontSize', 14, 'FontWeight', 'bold');

% Add legend
legend({'Production', 'Consumption'}, 'Location', 'northeast', 'FontSize', 11);

% Add grid for better readability
grid on;
set(gca, 'GridLineStyle', '--', 'GridAlpha', 0.3);

% Improve aesthetics
set(gca, 'FontSize', 15);
box on;

% Add value labels on top of each bar
% for i = 1:size(bar_data, 1)
%     for j = 1:size(bar_data, 2)
%         % Get the x position for each bar
%         x_pos = b(j).XEndPoints(i);
%         y_pos = bar_data(i, j);
%         
%         % Add text label
%         text(x_pos, y_pos, sprintf('%.4f', y_pos), ...
%             'HorizontalAlignment', 'center', ...
%             'VerticalAlignment', 'bottom', ...
%             'FontSize', 9);
%     end
% end

% Adjust y-axis limits to accommodate text labels
ylim([0, max(bar_data(:))*1.1]);

% Save the figure
fig_filename = fullfile(subfolder, 'NADH_Production_Consumption_BarGraph.png');
saveas(gcf, fig_filename);
fprintf('Bar graph saved to: %s\n', fig_filename);

% Also save as high-quality figure for publication
% fig_filename_hq = fullfile(subfolder, 'NADH_Production_Consumption_BarGraph_HQ.eps');
% saveas(gcf, fig_filename_hq, 'epsc');
% fprintf('High-quality figure saved to: %s\n', fig_filename_hq);



% Define figure size
width = 4; %4
height = 4; % 4

% Apply size to paper
set(gcf, 'PaperPosition', [0 0 width height]);
set(gcf, 'PaperSize', [width height]);

% Export to PDF
print(gcf, fullfile(subfolder, 'NADH_Production_Consumption_BarGraph.pdf'), '-dpdf', '-painters');
