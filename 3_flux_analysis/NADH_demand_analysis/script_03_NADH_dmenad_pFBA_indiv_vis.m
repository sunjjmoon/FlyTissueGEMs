%% script_03_NADH_dmenad_pFBA_indiv_vis: Comprehensive NADH-dependent reactions analysis
% Analyzes all NADH-producing and NADH-consuming reactions from pFBA
% results under this NADH demand maximization simulation
% Verifies mass balance with NADH production capacity from Script 02

clc;
close all;
clear;

%% 1. Load Results from Script 01
results_dir = fullfile(pwd, '01_NADH_demand_pFBA');
load(fullfile(results_dir, 'enhanced_FBA_results.mat'));
load(fullfile(results_dir, 'model_constrained_out.mat'));

output_dir = fullfile(pwd, '03_NADH_demand_pFBA_vis_indiv');
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

%% 2. Define NADH/NAD+ metabolites
nadh_mets_cyto = {'MAM02553c[c]'};  % NADH cytosolic
nadh_mets_mito = {'MAM02553m[m]'};  % NADH mitochondrial

nad_mets_cyto = {'MAM02552c[c]'};   % NAD+ cytosolic  
nad_mets_mito = {'MAM02552m[m]'};   % NAD+ mitochondrial

nadh_mets_p = {'MAM02553p[p]'};   % NAD+ cytosolic  
nadh_mets_r = {'MAM02553r[r]'};   % NAD+ mitochondrial

nad_mets_p = {'MAM02552p[p]'};   % NAD+ cytosolic  
nad_mets_r = {'MAM02552r[r]'};   % NAD+ mitochondrial

all_nadh_mets = [nadh_mets_cyto, nadh_mets_mito,nadh_mets_p];
all_nad_mets = [nad_mets_cyto, nad_mets_mito,nad_mets_p];

%% 3. Extract NADH-dependent reactions for each model
n_models = length(enhanced_results);
model_names = cell(n_models, 1);
nadh_analysis = cell(n_models, 1);

for i = 1:n_models
    model_names{i} = enhanced_results{i}.modelID;
    model = model_constrained_out{i,1};
    biomass_solution = enhanced_results{i}.biomass_solution;
    nad_demand_solution = enhanced_results{i}.nad_demand_solution;
    
    fprintf('Processing model %s (%d/%d)\n', model_names{i}, i, n_models);
    
    if biomass_solution.stat == 1
        flux_values = biomass_solution.x;
        flux_values = nad_demand_solution.x;

    else
        flux_values = zeros(length(model.rxns), 1);
        fprintf('Warning: No valid solution for model %s\n', model_names{i});
        continue;
    end
    
    %% Find all reactions involving NADH/NAD+
    nadh_rxn_indices = [];
    
    % Find reactions involving NADH (cytosolic and mitochondrial)
    for met = 1:length(all_nadh_mets)
        met_idx = find(strcmp(model.mets, all_nadh_mets{met}));
        if ~isempty(met_idx)
            rxn_indices = find(model.S(met_idx, :) ~= 0);
            nadh_rxn_indices = [nadh_rxn_indices, rxn_indices];
        end
    end
    
    % Find reactions involving NAD+
    for met = 1:length(all_nad_mets)
        met_idx = find(strcmp(model.mets, all_nad_mets{met}));
        if ~isempty(met_idx)
            rxn_indices = find(model.S(met_idx, :) ~= 0);
            nadh_rxn_indices = [nadh_rxn_indices, rxn_indices];
        end
    end
    
    nadh_rxn_indices = unique(nadh_rxn_indices);
    
    fprintf('  Found %d NADH/NAD+-related reactions\n', length(nadh_rxn_indices));
    
    %% Categorize reactions by NADH production/consumption
    nadh_producing_rxns = {};
    nadh_consuming_rxns = {};
    nadh_producing_fluxes = [];
    nadh_consuming_fluxes = [];
    nadh_producing_stoich = [];
    nadh_consuming_stoich = [];
    
    total_nadh_production = 0;
    total_nadh_consumption = 0;
    
    for j = 1:length(nadh_rxn_indices)
        rxn_idx = nadh_rxn_indices(j);
        flux = flux_values(rxn_idx);
        
        % Calculate net NADH change for this reaction
        net_nadh_stoich = 0;
        
        % Check NADH stoichiometry (both compartments)
        for met = 1:length(all_nadh_mets)
            met_idx = find(strcmp(model.mets, all_nadh_mets{met}));
            if ~isempty(met_idx)
                net_nadh_stoich = net_nadh_stoich + model.S(met_idx, rxn_idx);
            end
        end
        
        % Check NAD+ stoichiometry (negative contribution to NADH)
%         for met = 1:length(all_nad_mets)
%             met_idx = find(strcmp(model.mets, all_nad_mets{met}));
%             if ~isempty(met_idx)
%                 net_nadh_stoich = net_nadh_stoich - model.S(met_idx, rxn_idx);
%             end
%         end
        
        net_nadh_flux = net_nadh_stoich * flux;
        
        % Categorize based on net NADH production/consumption
        if abs(net_nadh_flux) > 1e-9  % Only consider significant fluxes
            if net_nadh_flux > 0  % Net NADH production
                nadh_producing_rxns{end+1} = model.rxns{rxn_idx};
                nadh_producing_fluxes(end+1) = abs(net_nadh_flux);
                nadh_producing_stoich(end+1) = net_nadh_stoich;
                total_nadh_production = total_nadh_production + abs(net_nadh_flux);
            else  % Net NADH consumption
                nadh_consuming_rxns{end+1} = model.rxns{rxn_idx};
                nadh_consuming_fluxes(end+1) = abs(net_nadh_flux);
                nadh_consuming_stoich(end+1) = net_nadh_stoich;
                total_nadh_consumption = total_nadh_consumption + abs(net_nadh_flux);
            end
        end
    end
    
    %% Store analysis results
    nadh_analysis{i}.model_id = model_names{i};
    nadh_analysis{i}.nadh_producing_rxns = nadh_producing_rxns;
    nadh_analysis{i}.nadh_producing_fluxes = nadh_producing_fluxes;
    nadh_analysis{i}.nadh_producing_stoich = nadh_producing_stoich;
    nadh_analysis{i}.nadh_consuming_rxns = nadh_consuming_rxns;
    nadh_analysis{i}.nadh_consuming_fluxes = nadh_consuming_fluxes;
    nadh_analysis{i}.nadh_consuming_stoich = nadh_consuming_stoich;
    nadh_analysis{i}.total_nadh_production = total_nadh_production;
    nadh_analysis{i}.total_nadh_consumption = total_nadh_consumption;
    nadh_analysis{i}.mass_balance_error = abs(total_nadh_production - total_nadh_consumption);
    nadh_analysis{i}.nadh_oxidation_capacity = enhanced_results{i}.max_nadh_oxidation_capacity;
    
    fprintf('  NADH production: %.3f mmol/gDW/h\n', total_nadh_production);
    fprintf('  NADH consumption: %.3f mmol/gDW/h\n', total_nadh_consumption);
    fprintf('  Mass balance error: %.6f\n', nadh_analysis{i}.mass_balance_error);
    fprintf('  NADH oxidation capacity: %.3f mmol/gDW/h\n', nadh_analysis{i}.nadh_oxidation_capacity);
end

%% 4. Create comprehensive summary tables
% NADH-producing reactions table
if n_models >= 1
    all_producing_rxns = unique([nadh_analysis{1}.nadh_producing_rxns, ...
                                nadh_analysis{2}.nadh_producing_rxns]);
    
    producing_table = table();
    producing_table.ReactionID = all_producing_rxns';
    
    for i = 1:n_models
        model_fluxes = zeros(length(all_producing_rxns), 1);
        for j = 1:length(all_producing_rxns)
            rxn_idx = find(strcmp(nadh_analysis{i}.nadh_producing_rxns, all_producing_rxns{j}));
            if ~isempty(rxn_idx)
                model_fluxes(j) = nadh_analysis{i}.nadh_producing_fluxes(rxn_idx);
            end
        end
        producing_table.([model_names{i} '_NADH_Production']) = model_fluxes;
    end
    
    if n_models == 2
        fold_changes = producing_table.([model_names{2} '_NADH_Production']) ./ ...
                      producing_table.([model_names{1} '_NADH_Production']);
        fold_changes(isinf(fold_changes) | isnan(fold_changes)) = 0;
        producing_table.Fold_Change_HSD_vs_NSD = fold_changes;
    end
    
end

% NADH-consuming reactions table  
if n_models >= 1
    all_consuming_rxns = unique([nadh_analysis{1}.nadh_consuming_rxns, ...
                                nadh_analysis{2}.nadh_consuming_rxns]);
    
    consuming_table = table();
    consuming_table.ReactionID = all_consuming_rxns';
    
    for i = 1:n_models
        model_fluxes = zeros(length(all_consuming_rxns), 1);
        for j = 1:length(all_consuming_rxns)
            rxn_idx = find(strcmp(nadh_analysis{i}.nadh_consuming_rxns, all_consuming_rxns{j}));
            if ~isempty(rxn_idx)
                model_fluxes(j) = nadh_analysis{i}.nadh_consuming_fluxes(rxn_idx);
            end
        end
        consuming_table.([model_names{i} '_NADH_Consumption']) = model_fluxes;
    end
    
    if n_models == 2
        fold_changes = consuming_table.([model_names{2} '_NADH_Consumption']) ./ ...
                      consuming_table.([model_names{1} '_NADH_Consumption']);
        fold_changes(isinf(fold_changes) | isnan(fold_changes)) = 0;
        consuming_table.Fold_Change_HSD_vs_NSD = fold_changes;
    end
    
end

%% Export Raw Reaction Flux Values (just flux, no stoichiometry)
if n_models == 2
    fprintf('\n=== Exporting raw reaction flux values ===\n');
    
    % Get ALL reactions (both producing and consuming combined)
    all_reactions = unique([nadh_analysis{1}.nadh_producing_rxns, ...
                           nadh_analysis{1}.nadh_consuming_rxns, ...
                           nadh_analysis{2}.nadh_producing_rxns, ...
                           nadh_analysis{2}.nadh_consuming_rxns]);
    
    combined_table = table();
    combined_table.ReactionID = all_reactions';
    
    for i = 1:n_models
        model_fluxes_raw = zeros(length(all_reactions), 1);
        model = model_constrained_out{i,1};
        flux_values = enhanced_results{i}.nad_demand_solution.x;
        
        for j = 1:length(all_reactions)
            rxn_idx = find(strcmp(model.rxns, all_reactions{j}));
            if ~isempty(rxn_idx)
                % Just the raw flux value, no stoichiometry multiplication
                model_fluxes_raw(j) = flux_values(rxn_idx);
            end
        end
        combined_table.([model_names{i} '_Raw_Flux']) = model_fluxes_raw;
    end
    
end

%% NADH Production vs Consumption grouped by condition
if n_models == 2
    figure(13);
    categories = {'NSD', 'HSD'}; % x-axis = condition
    
    % Values for each condition
    nsd_values = [nadh_analysis{1}.total_nadh_production, ...
                  nadh_analysis{1}.total_nadh_consumption];
    
    hsd_values = [nadh_analysis{2}.total_nadh_production, ...
                  nadh_analysis{2}.total_nadh_consumption];
    
    % Arrange as rows = condition, cols = [Production, Consumption]
    data_matrix = [nsd_values; hsd_values];
    
    % Plot grouped bars (x = NSD/HSD)
    bar_handle = bar(data_matrix, 'grouped');
    bar_handle(1).FaceColor = [0.3, 0.6, 0.9]; % Production - Blue
    bar_handle(2).FaceColor = [0.9, 0.4, 0.3]; % Consumption - Red
    
    set(gca, 'XTickLabel', categories, 'FontSize', 16);
    ylabel('Flux (a.u.)', 'FontSize', 18, 'FontWeight', 'bold');
    ylim([0 1600])
%     title('NADH production vs oxidation', ...
%           'FontSize', 14, 'FontWeight', 'bold');
    legend({'Production', 'Consumption'}, 'Location', 'eastoutside', 'FontSize', 11);
    grid on;
    box on;
    
    saveas(gcf, fullfile(output_dir, 'NADH_Prod_Cons_byCondition.fig'));
    saveas(gcf, fullfile(output_dir, 'NADH_Prod_Cons_byCondition.png'));
    print(gcf, fullfile(output_dir, 'NADH_Prod_Cons_byCondition.svg'), '-dsvg');
    
    % ---- Export PDF with true figure dimensions ----
    set(gcf, 'Units', 'inches');
    fig_pos = get(gcf, 'Position');               % [x y width height] in inches
    set(gcf, 'PaperUnits', 'inches');
    set(gcf, 'PaperPosition', [0 0 fig_pos(3) fig_pos(4)]);  % match figure size
    set(gcf, 'PaperSize', [fig_pos(3) fig_pos(4)]);           % remove margins
    
    print(gcf, fullfile(output_dir, 'NADH_Prod_Cons_byCondition.pdf'), ...
          '-dpdf', '-r300');  % high-resolution PDF export
end




%% Change TSG101 to LDH for multiple reactions
model = model_constrained_out{1,1};

% List of reactions to change from TSG101 to LDH
ldh_reactions = {'MAR04388', 'MAR04280', 'MAR04281'};

for i = 1:length(ldh_reactions)
    rxn_idx = find(strcmp(model.rxns, ldh_reactions{i}));
    if ~isempty(rxn_idx)
        model.grRules{rxn_idx} = 'LDH';
        fprintf('Changed %s gene association from TSG101 to LDH\n', ldh_reactions{i});
    else
        fprintf('Reaction %s not found\n', ldh_reactions{i});
    end
end

% Update the model in both model slots
model_constrained_out{1,1} = model;
if length(model_constrained_out) > 1
    model2 = model_constrained_out{2,1};
    for i = 1:length(ldh_reactions)
        rxn_idx = find(strcmp(model2.rxns, ldh_reactions{i}));
        if ~isempty(rxn_idx)
            model2.grRules{rxn_idx} = 'LDH';
        end
    end
    model_constrained_out{2,1} = model2;
end

%% Figure 24: Cumulative NADH Production by Gene
if n_models == 2
    figure(24);
    
    % Get all producing reactions and their fluxes
    all_producing_rxns = unique([nadh_analysis{1}.nadh_producing_rxns, ...
                                nadh_analysis{2}.nadh_producing_rxns]);
    
    nsd_prod_fluxes = zeros(length(all_producing_rxns), 1);
    hsd_prod_fluxes = zeros(length(all_producing_rxns), 1);
    
    for j = 1:length(all_producing_rxns)
        % NSD fluxes
        nsd_idx = find(strcmp(nadh_analysis{1}.nadh_producing_rxns, all_producing_rxns{j}));
        if ~isempty(nsd_idx)
            nsd_prod_fluxes(j) = nadh_analysis{1}.nadh_producing_fluxes(nsd_idx);
        end
        
        % HSD fluxes
        hsd_idx = find(strcmp(nadh_analysis{2}.nadh_producing_rxns, all_producing_rxns{j}));
        if ~isempty(hsd_idx)
            hsd_prod_fluxes(j) = nadh_analysis{2}.nadh_producing_fluxes(hsd_idx);
        end
    end
    
    % Sort by maximum flux across both conditions (descending)
    max_fluxes = max([nsd_prod_fluxes, hsd_prod_fluxes], [], 2);
    [~, sort_idx] = sort(max_fluxes, 'descend');
    
    nsd_sorted = nsd_prod_fluxes(sort_idx);
    hsd_sorted = hsd_prod_fluxes(sort_idx);
    rxns_sorted = all_producing_rxns(sort_idx);
    
    % Extract gene names for sorted reactions
    model = model_constrained_out{1,1};
    gene_labels = cell(length(rxns_sorted), 1);
    
    for j = 1:length(rxns_sorted)
        rxn_idx = find(strcmp(model.rxns, rxns_sorted{j}));
        if ~isempty(rxn_idx)
            gene_rule = model.grRules{rxn_idx};
            
            % Extract compartment info from NADH metabolites in this reaction
            comp_str = '';
            for met = 1:length(all_nadh_mets)
                met_idx = find(strcmp(model.mets, all_nadh_mets{met}));
                if ~isempty(met_idx)
                    stoich = full(model.S(met_idx, rxn_idx));
                    if stoich ~= 0
                        % Extract compartment from metabolite name (e.g., [c], [m])
                        comp_match = regexp(all_nadh_mets{met}, '\[([a-z])\]', 'tokens');
                        if ~isempty(comp_match)
                            comp_str = sprintf('[%s]', comp_match{1}{1});
                            break;
                        end
                    end
                end
            end
            
            if ~isempty(gene_rule)
                % Extract first gene name (handle "or" cases)
                gene_rule = strrep(gene_rule, ' or ', ' ');
                genes = regexp(gene_rule, '\w+', 'match');
                if ~isempty(genes)
                    gene_labels{j} = sprintf('%s (%s) %s', genes{1}, rxns_sorted{j}, comp_str);
                else
                    gene_labels{j} = sprintf('%s %s', rxns_sorted{j}, comp_str);
                end
            else
                gene_labels{j} = sprintf('%s %s', rxns_sorted{j}, comp_str);
            end
        else
            gene_labels{j} = rxns_sorted{j};
        end
    end
    
    
    % Calculate cumulative percentages
    nsd_total = nadh_analysis{1}.total_nadh_production;
    hsd_total = nadh_analysis{2}.total_nadh_production;
    
    nsd_cumulative = cumsum(nsd_sorted) / nsd_total * 100;
    hsd_cumulative = cumsum(hsd_sorted) / hsd_total * 100;
    
    % Plot cumulative curves
    plot(1:length(gene_labels), nsd_cumulative, 'o-', 'Color', [0.5, 0.5, 0.5], ...
         'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', [0.5, 0.5, 0.5]);
    hold on;
    plot(1:length(gene_labels), hsd_cumulative, 's-', 'Color', [0.9, 0.4, 0.3], ...
         'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', [0.9, 0.4, 0.3]);
    
    % Add horizontal line at 80%

    xlabel('gene (rxn id)', 'FontSize', 14, 'FontWeight', 'bold');
%     ylabel('Cumulative % of total NADH production', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(['Cumulative %' newline 'of total NADH production'], 'FontSize', 16, 'FontWeight', 'bold');

%     title('NADH production', 'FontSize', 14, 'FontWeight', 'bold');
    legend({'NSD', 'HSD'}, 'Location', 'southeast', 'FontSize', 11);
 
    
    grid on;
    box on;
    xlim([0, length(gene_labels)+1]);
    ylim([0, 105]);
    
    % Set x-tick labels to gene names (rotated)
    set(gca, 'XTick', 1:length(gene_labels), 'XTickLabel', gene_labels, 'XTickLabelRotation', 45);
    
    saveas(gcf, fullfile(output_dir, 'Cumulative_NADH_Production_Genes.fig'));
    saveas(gcf, fullfile(output_dir, 'Cumulative_NADH_Production_Genes.png'));
    print(gcf, fullfile(output_dir, 'Cumulative_NADH_Production_Genes.svg'), '-dsvg');
    set(gcf, 'Units', 'inches');                 % work in inches
    fig_pos = get(gcf, 'Position');              % current figure size
    set(gcf, 'PaperUnits', 'inches', ...
             'PaperSize', fig_pos(3:4), ...      % match paper size to figure
             'PaperPosition', [0 0 fig_pos(3:4)]); % no margins
    % Save as PDF (vector graphics, no margins)
    print(gcf, fullfile(output_dir, 'Cumulative_NADH_Production_Genes.pdf'), '-dpdf', '-painters');
    
    
        %% Save gene-based data to Excel
    gene_cumulative_data = table();
    gene_cumulative_data.Gene_Reaction = gene_labels;
    gene_cumulative_data.NSD_Flux = nsd_sorted;
    gene_cumulative_data.HSD_Flux = hsd_sorted;
    gene_cumulative_data.NSD_Cumulative_Percent = nsd_cumulative;
    gene_cumulative_data.HSD_Cumulative_Percent = hsd_cumulative;
    
    writetable(gene_cumulative_data, fullfile(output_dir, 'Cumulative_NADH_Analysis_Genes_prod.xlsx'));
    
end

%% Figure 25: Cumulative NADH Consumption by Gene
if n_models == 2
    figure(25);
    
    % Get all consuming reactions and their fluxes (excluding NADH_demand)
    all_consuming_rxns = unique([nadh_analysis{1}.nadh_consuming_rxns, ...
                                nadh_analysis{2}.nadh_consuming_rxns]);
    
    % Remove NADH_demand from the list
    nadh_demand_idx = find(strcmp(all_consuming_rxns, 'NADH_demand'));
    if ~isempty(nadh_demand_idx)
        all_consuming_rxns(nadh_demand_idx) = [];
    end
    
    nsd_cons_fluxes = zeros(length(all_consuming_rxns), 1);
    hsd_cons_fluxes = zeros(length(all_consuming_rxns), 1);
    
    for j = 1:length(all_consuming_rxns)
        % NSD fluxes
        nsd_idx = find(strcmp(nadh_analysis{1}.nadh_consuming_rxns, all_consuming_rxns{j}));
        if ~isempty(nsd_idx)
            nsd_cons_fluxes(j) = nadh_analysis{1}.nadh_consuming_fluxes(nsd_idx);
        end
        
        % HSD fluxes
        hsd_idx = find(strcmp(nadh_analysis{2}.nadh_consuming_rxns, all_consuming_rxns{j}));
        if ~isempty(hsd_idx)
            hsd_cons_fluxes(j) = nadh_analysis{2}.nadh_consuming_fluxes(hsd_idx);
        end
    end
    
    % Sort by maximum flux across both conditions (descending)
    max_fluxes = max([nsd_cons_fluxes, hsd_cons_fluxes], [], 2);
    [~, sort_idx] = sort(max_fluxes, 'descend');
    
    nsd_sorted = nsd_cons_fluxes(sort_idx);
    hsd_sorted = hsd_cons_fluxes(sort_idx);
    rxns_sorted = all_consuming_rxns(sort_idx);
    
    % Extract gene names for sorted reactions
    gene_labels = cell(length(rxns_sorted), 1);
    
    for j = 1:length(rxns_sorted)
        rxn_idx = find(strcmp(model.rxns, rxns_sorted{j}));
        if ~isempty(rxn_idx)
            gene_rule = model.grRules{rxn_idx};
            if ~isempty(gene_rule)
                % Extract first gene name (handle "or" cases)
                gene_rule = strrep(gene_rule, ' or ', ' '); % Replace "or" with space
                genes = regexp(gene_rule, '\w+', 'match');
                if ~isempty(genes)
                    gene_labels{j} = sprintf('%s (%s)', genes{1}, rxns_sorted{j});
                else
                    gene_labels{j} = rxns_sorted{j};
                end
            else
                gene_labels{j} = rxns_sorted{j}; % No gene association
            end
        else
            gene_labels{j} = rxns_sorted{j};
        end
    end
    
    % Calculate cumulative percentages (excluding NADH_demand from totals)
    nsd_demand_flux = 0;
    hsd_demand_flux = 0;
    nsd_demand_idx = find(strcmp(nadh_analysis{1}.nadh_consuming_rxns, 'NADH_demand'));
    if ~isempty(nsd_demand_idx)
        nsd_demand_flux = nadh_analysis{1}.nadh_consuming_fluxes(nsd_demand_idx);
    end
    hsd_demand_idx = find(strcmp(nadh_analysis{2}.nadh_consuming_rxns, 'NADH_demand'));
    if ~isempty(hsd_demand_idx)
        hsd_demand_flux = nadh_analysis{2}.nadh_consuming_fluxes(hsd_demand_idx);
    end
    
    nsd_total_excluding_demand = nadh_analysis{1}.total_nadh_consumption - nsd_demand_flux;
    hsd_total_excluding_demand = nadh_analysis{2}.total_nadh_consumption - hsd_demand_flux;
    
    nsd_cumulative = cumsum(nsd_sorted) / nsd_total_excluding_demand * 100;
    hsd_cumulative = cumsum(hsd_sorted) / hsd_total_excluding_demand * 100;
    
    % Plot cumulative curves
    plot(1:length(gene_labels), nsd_cumulative, 'o-', 'Color', [0.5, 0.5, 0.5], ...
         'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', [0.5, 0.5, 0.5]);
    hold on;
    plot(1:length(gene_labels), hsd_cumulative, 's-', 'Color', [0.9, 0.4, 0.3], ...
         'LineWidth', 2, 'MarkerSize', 6, 'MarkerFaceColor', [0.9, 0.4, 0.3]);
    

    grid on;
    box on;
    xlim([0, length(gene_labels)+1]);
    ylim([0, 105]);
    
    % Set x-tick labels to gene names (rotated)
    set(gca, 'XTick', 1:length(gene_labels), 'XTickLabel', gene_labels, 'XTickLabelRotation', 45,'FontSize', 8);

         % Add horizontal line at 80%
    xlabel('gene (rxn id)', 'FontSize', 16, 'FontWeight', 'bold');
%     ylabel('Cumulative % of total NADH oxidation', 'FontSize', 12, 'FontWeight', 'bold');
    ylabel(['Cumulative %' newline 'of total NADH oxidation'], 'FontSize', 16, 'FontWeight', 'bold');

%     title('NADH oxidation', 'FontSize', 14, 'FontWeight', 'bold');
    legend({'NSD', 'HSD'}, 'Location', 'southeast', 'FontSize', 11);
     
    saveas(gcf, fullfile(output_dir, 'Cumulative_NADH_Consumption_Genes.fig'));
    saveas(gcf, fullfile(output_dir, 'Cumulative_NADH_Consumption_Genes.png'));
    print(gcf, fullfile(output_dir, 'Cumulative_NADH_Consumption_Genes.svg'), '-dsvg');
    set(gcf, 'Units', 'inches');                 % work in inches
    fig_pos = get(gcf, 'Position');              % current figure size
    set(gcf, 'PaperUnits', 'inches', ...
             'PaperSize', fig_pos(3:4), ...      % match paper size to figure
             'PaperPosition', [0 0 fig_pos(3:4)]); % no margins
    % Save as PDF (vector graphics, no margins)
    print(gcf, fullfile(output_dir, 'Cumulative_NADH_Consumption_Genes.pdf'), '-dpdf', '-painters');
    
    %% Save gene-based data to Excel
    gene_cumulative_data = table();
    gene_cumulative_data.Gene_Reaction = gene_labels;
    gene_cumulative_data.NSD_Flux = nsd_sorted;
    gene_cumulative_data.HSD_Flux = hsd_sorted;
    gene_cumulative_data.NSD_Cumulative_Percent = nsd_cumulative;
    gene_cumulative_data.HSD_Cumulative_Percent = hsd_cumulative;
    
    writetable(gene_cumulative_data, fullfile(output_dir, 'Cumulative_NADH_Analysis_Genes.xlsx'));
end