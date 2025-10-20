%% script_05_compare_all_export.m

clc;
clear;
close all;

%% Load data
current_path = pwd;
results_dir = fullfile(current_path, '04_sampling_comp');

load(fullfile(results_dir, 'fba_pfba_sampling_comparison.mat'));
load(fullfile(current_path,'01_results','model_constrained_out.mat'));


%% Get reaction names and annotations
rxn_names = model_constrained_out{1}.rxns; % Using NSD as reference
subsystems = model_constrained_out{1}.subSystems;
rxn_formulas = model_constrained_out{1}.rxnNames; % Reaction formulas/names

% Extract subsystem names (handle nested cell arrays)
subsystem_names = cell(length(subsystems), 1);
for i = 1:length(subsystems)
    if iscell(subsystems{i}) && ~isempty(subsystems{i})
        subsystem_names{i} = subsystems{i}{1};
    elseif ischar(subsystems{i}) || isstring(subsystems{i})
        subsystem_names{i} = char(subsystems{i});
    else
        subsystem_names{i} = 'Unknown';
    end
end

% Get metabolite information for reactions
met_names = model_constrained_out{1}.metNames;
met_ids = model_constrained_out{1}.mets; % This should contain MAM-style IDs
S_matrix = model_constrained_out{1}.S;

% Create metabolite lists for each reaction (both names and MAM IDs)
rxn_metabolites = cell(length(rxn_names), 1);
rxn_metabolite_ids = cell(length(rxn_names), 1);

for i = 1:length(rxn_names)
    % Find metabolites involved in this reaction
    met_idx = find(abs(S_matrix(:, i)) > 0);
    if ~isempty(met_idx) && length(met_idx) <= 8 % Limit to avoid huge lists
        % Get metabolite names
        met_list = met_names(met_idx);
        met_list_clean = met_list(~cellfun(@isempty, met_list));
        
        % Get metabolite IDs (MAM style)
        met_id_list = met_ids(met_idx);
        met_id_list_clean = met_id_list(~cellfun(@isempty, met_id_list));
        
        if ~isempty(met_list_clean)
            rxn_metabolites{i} = strjoin(met_list_clean, '; ');
        else
            rxn_metabolites{i} = 'Metabolites not available';
        end
        
        if ~isempty(met_id_list_clean)
            rxn_metabolite_ids{i} = strjoin(met_id_list_clean, '; ');
        else
            rxn_metabolite_ids{i} = 'IDs not available';
        end
    else
        rxn_metabolites{i} = 'Multiple metabolites';
        rxn_metabolite_ids{i} = 'Multiple IDs';
    end
end

%% Get gene information for reactions
genes = model_constrained_out{1}.genes;
rxnGeneMat = model_constrained_out{1}.rxnGeneMat;
grRules = model_constrained_out{1}.grRules;

% Create gene lists for each reaction
rxn_genes = cell(length(rxn_names), 1);
rxn_gene_rules = cell(length(rxn_names), 1);

for i = 1:length(rxn_names)
    % Find genes associated with this reaction
    gene_idx = find(rxnGeneMat(i, :) > 0);
    
    if ~isempty(gene_idx)
        gene_list = genes(gene_idx);
        rxn_genes{i} = strjoin(gene_list, '; ');
        
        if ~isempty(grRules{i})
            rxn_gene_rules{i} = grRules{i};
        else
            rxn_gene_rules{i} = 'No rule available';
        end
    else
        rxn_genes{i} = 'No genes associated';
        rxn_gene_rules{i} = 'No rule available';
    end
end

%% Export Venn Diagram Data to Excel for R Analysis
% Prepare data for Venn diagram creation in R

%% Calculate disease differences (if not already done)
% fba_diff = abs(comparison_data.HSD.fba_fluxes) - abs(comparison_data.NSD.fba_fluxes);
% pfba_diff = abs(comparison_data.HSD.pfba_fluxes) - abs(comparison_data.NSD.pfba_fluxes);
% sampling_diff = abs(comparison_data.HSD.sampling_median) - abs(comparison_data.NSD.sampling_median);

% fba_diff = ((comparison_data.HSD.fba_fluxes) - (comparison_data.NSD.fba_fluxes));
% pfba_diff = ((comparison_data.HSD.pfba_fluxes) - (comparison_data.NSD.pfba_fluxes));
% sampling_diff = ((comparison_data.HSD.sampling_median) - (comparison_data.NSD.sampling_median));


% Initialize difference arrays
fba_diff = zeros(size(comparison_data.HSD.fba_fluxes));
pfba_diff = zeros(size(comparison_data.HSD.pfba_fluxes));
sampling_diff = zeros(size(comparison_data.HSD.sampling_median));

% Process each reaction for FBA
for i = 1:length(fba_diff)
    nsd_flux = comparison_data.NSD.fba_fluxes(i);
    hsd_flux = comparison_data.HSD.fba_fluxes(i);
    
    % Check if direction reversed (different signs)
    if (nsd_flux * hsd_flux) < 0
        % Reversal: use signed difference
        fba_diff(i) = hsd_flux - nsd_flux;
    else
        % No reversal: use magnitude difference
        fba_diff(i) = abs(hsd_flux) - abs(nsd_flux);
    end
end

% Process each reaction for pFBA
for i = 1:length(pfba_diff)
    nsd_flux = comparison_data.NSD.pfba_fluxes(i);
    hsd_flux = comparison_data.HSD.pfba_fluxes(i);
    
    if (nsd_flux * hsd_flux) < 0
        pfba_diff(i) = hsd_flux - nsd_flux;
    else
        pfba_diff(i) = abs(hsd_flux) - abs(nsd_flux);
    end
end

% Process each reaction for Sampling
for i = 1:length(sampling_diff)
    nsd_flux = comparison_data.NSD.sampling_median(i);
    hsd_flux = comparison_data.HSD.sampling_median(i);
    
    if (nsd_flux * hsd_flux) < 0
        sampling_diff(i) = hsd_flux - nsd_flux;
    else
        sampling_diff(i) = abs(hsd_flux) - abs(nsd_flux);
    end
end


% Set threshold for significant changes
threshold = 1;

% Get reaction names
rxn_names = model_constrained_out{1}.rxns; % Using NSD as reference

%% Create binary classification matrices
% For "increased" reactions (HSD > NSD)
fba_increased = fba_diff > threshold;
pfba_increased = pfba_diff > threshold;
sampling_increased = sampling_diff > threshold;

% For "decreased" reactions (HSD < NSD)
fba_decreased = fba_diff < -threshold;
pfba_decreased = pfba_diff < -threshold;
sampling_decreased = sampling_diff < -threshold;

%% Export 1: Increased Reactions Data (Enhanced)
% Create table for increased reactions with subsystem and metabolite info
increased_data = table();
increased_data.Reaction_ID = rxn_names;
increased_data.Subsystem = subsystem_names;
increased_data.Reaction_Formula = rxn_formulas;
increased_data.Associated_Metabolites = rxn_metabolites;
increased_data.Metabolite_IDs = rxn_metabolite_ids;  % New MAM-style column
increased_data.Associated_Genes = rxn_genes;         % NEW: Gene list
increased_data.Gene_Rules = rxn_gene_rules;          % NEW: Gene rules
increased_data.FBA_Difference = fba_diff;
increased_data.pFBA_Difference = pfba_diff;
increased_data.Sampling_Difference = sampling_diff;
increased_data.FBA_Increased = fba_increased;
increased_data.pFBA_Increased = pfba_increased;
increased_data.Sampling_Increased = sampling_increased;

% Export to Excel
venn_excel_file = fullfile(results_dir, 'venn_diagram_data.xlsx');
writetable(increased_data, venn_excel_file, 'Sheet', 'Increased_Reactions');

%% Export 2: Decreased Reactions Data (Enhanced)
decreased_data = table();
decreased_data.Reaction_ID = rxn_names;
decreased_data.Subsystem = subsystem_names;
decreased_data.Reaction_Formula = rxn_formulas;
decreased_data.Associated_Metabolites = rxn_metabolites;
decreased_data.Metabolite_IDs = rxn_metabolite_ids;  % New MAM-style column
decreased_data.Associated_Genes = rxn_genes;         % NEW: Gene list
decreased_data.Gene_Rules = rxn_gene_rules;          % NEW: Gene rules
decreased_data.FBA_Difference = fba_diff;
decreased_data.pFBA_Difference = pfba_diff;
decreased_data.Sampling_Difference = sampling_diff;
decreased_data.FBA_Decreased = fba_decreased;
decreased_data.pFBA_Decreased = pfba_decreased;
decreased_data.Sampling_Decreased = sampling_decreased;

writetable(decreased_data, venn_excel_file, 'Sheet', 'Decreased_Reactions');

%% Export 3: Summary Statistics for R
% Create summary table
summary_stats = table();

% For increased reactions
fba_inc_count = sum(fba_increased);
pfba_inc_count = sum(pfba_increased);
sampling_inc_count = sum(sampling_increased);

% For decreased reactions
fba_dec_count = sum(fba_decreased);
pfba_dec_count = sum(pfba_decreased);
sampling_dec_count = sum(sampling_decreased);

% Overlaps for increased reactions
fba_pfba_inc = sum(fba_increased & pfba_increased);
fba_sampling_inc = sum(fba_increased & sampling_increased);
pfba_sampling_inc = sum(pfba_increased & sampling_increased);
all_three_inc = sum(fba_increased & pfba_increased & sampling_increased);

% Overlaps for decreased reactions
fba_pfba_dec = sum(fba_decreased & pfba_decreased);
fba_sampling_dec = sum(fba_decreased & sampling_decreased);
pfba_sampling_dec = sum(pfba_decreased & sampling_decreased);
all_three_dec = sum(fba_decreased & pfba_decreased & sampling_decreased);

% Create summary table
summary_stats = table(...
    {'FBA_Increased'; 'pFBA_Increased'; 'Sampling_Increased'; 'FBA_Decreased'; 'pFBA_Decreased'; 'Sampling_Decreased'; ...
     'FBA_pFBA_Increased'; 'FBA_Sampling_Increased'; 'pFBA_Sampling_Increased'; 'All_Three_Increased'; ...
     'FBA_pFBA_Decreased'; 'FBA_Sampling_Decreased'; 'pFBA_Sampling_Decreased'; 'All_Three_Decreased'}, ...
    [fba_inc_count; pfba_inc_count; sampling_inc_count; fba_dec_count; pfba_dec_count; sampling_dec_count; ...
     fba_pfba_inc; fba_sampling_inc; pfba_sampling_inc; all_three_inc; ...
     fba_pfba_dec; fba_sampling_dec; pfba_sampling_dec; all_three_dec], ...
    'VariableNames', {'Category', 'Count'});

writetable(summary_stats, venn_excel_file, 'Sheet', 'Summary_Counts');

%% Export 4: Reaction Lists for Each Set (for easier R processing)
% Lists of reaction IDs for each category

% Increased reactions only
fba_inc_list = rxn_names(fba_increased);
pfba_inc_list = rxn_names(pfba_increased);
sampling_inc_list = rxn_names(sampling_increased);

% Decreased reactions only
fba_dec_list = rxn_names(fba_decreased);
pfba_dec_list = rxn_names(pfba_decreased);
sampling_dec_list = rxn_names(sampling_decreased);

% Pad shorter lists with empty cells for equal-length columns
max_inc_length = max([length(fba_inc_list), length(pfba_inc_list), length(sampling_inc_list)]);
max_dec_length = max([length(fba_dec_list), length(pfba_dec_list), length(sampling_dec_list)]);

% Pad increased lists
fba_inc_padded = [fba_inc_list; cell(max_inc_length - length(fba_inc_list), 1)];
pfba_inc_padded = [pfba_inc_list; cell(max_inc_length - length(pfba_inc_list), 1)];
sampling_inc_padded = [sampling_inc_list; cell(max_inc_length - length(sampling_inc_list), 1)];

% Pad decreased lists
fba_dec_padded = [fba_dec_list; cell(max_dec_length - length(fba_dec_list), 1)];
pfba_dec_padded = [pfba_dec_list; cell(max_dec_length - length(pfba_dec_list), 1)];
sampling_dec_padded = [sampling_dec_list; cell(max_dec_length - length(sampling_dec_list), 1)];

% Create tables
increased_lists = table(fba_inc_padded, pfba_inc_padded, sampling_inc_padded, ...
    'VariableNames', {'FBA_Increased', 'pFBA_Increased', 'Sampling_Increased'});
decreased_lists = table(fba_dec_padded, pfba_dec_padded, sampling_dec_padded, ...
    'VariableNames', {'FBA_Decreased', 'pFBA_Decreased', 'Sampling_Decreased'});

writetable(increased_lists, venn_excel_file, 'Sheet', 'Increased_Lists');
writetable(decreased_lists, venn_excel_file, 'Sheet', 'Decreased_Lists');

%% Export 5: R-Ready Format (Enhanced with annotations)
% Create enhanced binary matrix for easy R import
r_ready_data = table();
r_ready_data.Reaction_ID = rxn_names;
r_ready_data.Subsystem = subsystem_names;
r_ready_data.Reaction_Formula = rxn_formulas;
r_ready_data.Associated_Metabolites = rxn_metabolites;
r_ready_data.Metabolite_IDs = rxn_metabolite_ids;        % Metabolite IDs
r_ready_data.Associated_Genes = rxn_genes;               % NEW: Gene list
r_ready_data.Gene_Rules = rxn_gene_rules;                % NEW: Gene rules
r_ready_data.FBA_Difference = fba_diff;
r_ready_data.pFBA_Difference = pfba_diff;
r_ready_data.Sampling_Difference = sampling_diff;
r_ready_data.FBA_Inc = double(fba_increased);
r_ready_data.pFBA_Inc = double(pfba_increased);
r_ready_data.Sampling_Inc = double(sampling_increased);
r_ready_data.FBA_Dec = double(fba_decreased);
r_ready_data.pFBA_Dec = double(pfba_decreased);
r_ready_data.Sampling_Dec = double(sampling_decreased);

writetable(r_ready_data, venn_excel_file, 'Sheet', 'R_Binary_Matrix');

%% Print summary for user
disp('=== Venn Diagram Data Export Complete ===');
disp(['File saved to: ' venn_excel_file]);
disp(' ');
disp('Excel sheets created:');
disp('1. Increased_Reactions - Full data for reactions increased in HSD');
disp('2. Decreased_Reactions - Full data for reactions decreased in HSD');
disp('3. Summary_Counts - Overlap counts for Venn diagrams');
disp('4. Increased_Lists - Lists of reaction IDs for each method (increased)');
disp('5. Decreased_Lists - Lists of reaction IDs for each method (decreased)');
disp('6. R_Binary_Matrix - Simple binary format for R analysis');
disp(' ');
disp('For R Venn diagrams, you can use:');
disp('- VennDiagram package');
disp('- ggvenn package');
disp('- UpSetR package (for more complex intersections)');
disp(' ');
disp('Summary statistics:');
disp(['Increased in disease - FBA: ' num2str(fba_inc_count) ', pFBA: ' num2str(pfba_inc_count) ', Sampling: ' num2str(sampling_inc_count)]);
disp(['Decreased in disease - FBA: ' num2str(fba_dec_count) ', pFBA: ' num2str(pfba_dec_count) ', Sampling: ' num2str(sampling_dec_count)]);
disp(['All three agree (increased): ' num2str(all_three_inc)]);
disp(['All three agree (decreased): ' num2str(all_three_dec)]);
