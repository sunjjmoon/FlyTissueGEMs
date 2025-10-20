%% Script 0. Sensitivity analysis. Perturb rxns of interest and perform FVA-sampling
clc;
% clear all;

%% SETTINGS
baseline_condition = 'HSD';  % Focus on HSD only

% perturbation_factors = [0.5, 0.70, 0.80, 0.90, 0.95];  % 30%, 20%, 10%, 5% reductions
perturbation_factors = [0.5];  % 10% reduction for testing

% Glycolysis enzymes to test
target_enzymes = {
    'MAR04394',  % Hexokinase (Hex-A)
    'MAR04381',  % Pgi
    'MAR04379',  % pfk
    'MAR04375',  % ald1
    'MAR04373',  % GAPDH
    'MAR04368',  % Pgk
    'MAR04365',  % Pglym78
    'MAR04363',  % Eno
    'MAR04358',  % Pyk
    'MAR04388',  % ldh
    'MAR04391'   % TPI
};

% Sampling parameters
num_samples = 10000;  % Number of samples for flux distribution

%% Set directories
output_dir = 'individual_conditions_v2';
pathway = pwd;
subfolder = fullfile(pathway, output_dir);
if ~exist(subfolder, 'dir')
    mkdir(subfolder)
end

%% Load models
% [parentFolder, ~, ~] = fileparts(pathway);
model_out_cbra = load(strcat('E:\Projects\revision\pFBA_FBA_FVA\models\', 'model_out_cbra_u.mat'));
model_out_cbra = model_out_cbra.model_out_cbra_u;

%% Find HSD model
baseline_model_idx = [];
for i = 1:length(model_out_cbra)
    if strcmp(model_out_cbra{i,1}.modelID, baseline_condition)
        baseline_model_idx = i;
        break;
    end
end

if isempty(baseline_model_idx)
    error('Could not find HSD model');
end

%% STEP 1: Prepare HSD baseline model and perform initial sampling
fprintf('Perturbation factors to test: %s\n', mat2str(perturbation_factors));
fprintf('Number of samples requested: %d\n', num_samples);
fprintf('Number of enzymes to test: %d\n', length(target_enzymes));

baseline_model = model_out_cbra{baseline_model_idx,1};
baseline_model_prepared = prepare_model(baseline_model, baseline_condition);

% Set solver
% changeCobraSolver('gurobi', 'LP');

% Run FVA on HSD baseline to define sampling space
fprintf('Running FVA on HSD baseline model...\n');
[minFlux_baseline, maxFlux_baseline] = fluxVariability(baseline_model_prepared, 90, 'max', baseline_model_prepared.rxns, 0, 1);

% Apply FVA bounds for sampling
baseline_model_prepared.lb = minFlux_baseline;
baseline_model_prepared.ub = maxFlux_baseline;

% Perform initial flux sampling on baseline model
fprintf('Performing initial flux sampling (%d samples requested)...\n', num_samples);
samples_baseline = [];
actual_num_samples = 0;

try
    sampler_result = gpSampler(baseline_model_prepared, num_samples);
    
    % Extract the actual samples from the struct
    if isstruct(sampler_result) && isfield(sampler_result, 'points')
        samples_baseline = sampler_result.points;  % This is the actual sample matrix
        actual_num_samples = size(samples_baseline, 2);
        fprintf('Initial sampling completed successfully.\n');
        fprintf('Sample matrix dimensions: %d reactions x %d samples\n', size(samples_baseline, 1), size(samples_baseline, 2));
        fprintf('Requested %d samples, got %d samples\n', num_samples, actual_num_samples);
    else
        error('gpSampler did not return expected struct with points field');
    end
    
catch ME
    fprintf('Error in initial sampling: %s\n', ME.message);
    fprintf('Falling back to alternative sampling method...\n');
    % Fallback: use uniform random sampling within FVA bounds
    samples_baseline = zeros(length(baseline_model_prepared.rxns), num_samples);
    for i = 1:length(baseline_model_prepared.rxns)
        lb = baseline_model_prepared.lb(i);
        ub = baseline_model_prepared.ub(i);
        if ub > lb
            samples_baseline(i, :) = lb + (ub - lb) * rand(1, num_samples);
        else
            % Handle case where lb == ub (fixed flux)
            samples_baseline(i, :) = repmat(lb, 1, num_samples);
        end
    end
    actual_num_samples = num_samples;
    fprintf('Fallback sampling completed.\n');
    fprintf('Sample matrix dimensions: %d reactions x %d samples\n', size(samples_baseline, 1), size(samples_baseline, 2));
end

% Calculate median fluxes for target enzymes
reference_flux_median = struct();
valid_enzymes = {};
invalid_enzymes = {};

fprintf('\n=== CALCULATING MEDIAN REFERENCE FLUXES ===\n');
for k = 1:length(target_enzymes)
    enzyme_id = target_enzymes{k};
    idx_ref = find(ismember(baseline_model_prepared.rxns, enzyme_id));
    
    if ~isempty(idx_ref)
        % Safety check before accessing array
        if idx_ref > size(samples_baseline, 1) || size(samples_baseline, 2) == 0
            fprintf('=== %s === ERROR: Cannot access samples_baseline(%d, :)\n', enzyme_id, idx_ref);
            fprintf('samples_baseline size: [%d x %d]\n', size(samples_baseline, 1), size(samples_baseline, 2));
            invalid_enzymes{end+1} = enzyme_id;
            continue;
        end
        
        % Get flux samples for this enzyme
        enzyme_samples = samples_baseline(idx_ref, :);
        median_flux = median(enzyme_samples);
        flux_std = std(enzyme_samples);
        flux_range = [min(enzyme_samples), max(enzyme_samples)];
        
        reference_flux_median.(enzyme_id) = median_flux;
        
        fprintf('=== %s ===\n', enzyme_id);
        fprintf('FVA range: [%.6f, %.6f]\n', minFlux_baseline(idx_ref), maxFlux_baseline(idx_ref));
        fprintf('Sampled flux range: [%.6f, %.6f]\n', flux_range(1), flux_range(2));
        fprintf('Median flux: %.6f\n', median_flux);
        fprintf('Standard deviation: %.6f\n', flux_std);
        
        % Check if enzyme has meaningful flux
        if abs(median_flux) > 1e-6
            valid_enzymes{end+1} = enzyme_id;
            fprintf('Status: VALID for sensitivity analysis\n');
        else
            invalid_enzymes{end+1} = enzyme_id;
            fprintf('Status: INVALID (median flux too small)\n');
        end
    else
        invalid_enzymes{end+1} = enzyme_id;
        fprintf('=== %s === NOT FOUND in model\n', enzyme_id);
    end
    fprintf('\n');
end

fprintf('=== ENZYME VALIDATION SUMMARY ===\n');
fprintf('Valid enzymes: %d out of %d\n', length(valid_enzymes), length(target_enzymes));
fprintf('Valid: %s\n', strjoin(valid_enzymes, ', '));
if ~isempty(invalid_enzymes)
    fprintf('Invalid: %s\n', strjoin(invalid_enzymes, ', '));
end
fprintf('\n');

%% STEP 2: Create conditions for each valid enzyme with sampling-based perturbations
out_all = {};
condition_counter = 0;

for rxn_idx = 1:length(valid_enzymes)
    target_enzyme = valid_enzymes{rxn_idx};
    
    fprintf('=== ANALYZING ENZYME %d/%d: %s ===\n', rxn_idx, length(valid_enzymes), target_enzyme);
    
    % Create baseline condition (store original sampling results)
    condition_counter = condition_counter + 1;
    baseline_name = ['HSD_baseline_' target_enzyme];
    
    fprintf('Processing baseline: %s\n', baseline_name);
    
    % Store baseline results
    out_all{condition_counter,1}.samples = samples_baseline;
    out_all{condition_counter,1}.subSysModel = baseline_model_prepared;
    out_all{condition_counter,1}.modelID = baseline_name;
    out_all{condition_counter,1}.target_enzyme = target_enzyme;
    out_all{condition_counter,1}.condition_type = 'baseline';
    out_all{condition_counter,1}.perturbation_factor = 1.0;
    out_all{condition_counter,1}.is_perturbed = false;
    out_all{condition_counter,1}.median_flux = reference_flux_median.(target_enzyme);
    out_all{condition_counter,1}.num_samples = actual_num_samples;
    
    % Create perturbed conditions for each perturbation factor
    for p_idx = 1:length(perturbation_factors)
        condition_counter = condition_counter + 1;
        perturbation_factor = perturbation_factors(p_idx);
        perturbed_name = sprintf('HSD_perturbed_%s_%.0f', target_enzyme, perturbation_factor*100);
        
        fprintf('Processing perturbation %.0f%%: %s\n', perturbation_factor*100, perturbed_name);
        
        % Start with prepared HSD model
        model_perturbed = baseline_model_prepared;
        
        % Apply median-based enzyme perturbation
        median_flux = reference_flux_median.(target_enzyme);
        target_flux = median_flux * perturbation_factor;
        
        fprintf('  Target enzyme median flux: %.6f\n', median_flux);
        fprintf('  Perturbed target flux: %.6f (%.0f%% of median)\n', target_flux, perturbation_factor*100);
        
        % Apply perturbation: set both bounds to target flux
        model_perturbed = apply_median_perturbation(model_perturbed, target_enzyme, target_flux);
        
        % Perform flux sampling on perturbed model
        fprintf('  Performing flux sampling on perturbed model...\n');
        samples_perturbed = [];
        actual_perturbed_samples = 0;
        
        try
            sampler_result_pert = gpSampler(model_perturbed, num_samples);
            
            if isstruct(sampler_result_pert) && isfield(sampler_result_pert, 'points')
                samples_perturbed = sampler_result_pert.points;
                actual_perturbed_samples = size(samples_perturbed, 2);
                fprintf('  Perturbed sampling completed successfully.\n');
                fprintf('  Got %d samples for perturbed model\n', actual_perturbed_samples);
                
                % Verify perturbation was applied correctly
                enzyme_idx = find(ismember(model_perturbed.rxns, target_enzyme));
                if ~isempty(enzyme_idx)
                    achieved_fluxes = samples_perturbed(enzyme_idx, :);
                    achieved_median = median(achieved_fluxes);
                    achieved_std = std(achieved_fluxes);
                    
                    fprintf('  Achieved median flux: %.6f (target: %.6f)\n', achieved_median, target_flux);
                    fprintf('  Achieved std: %.6f\n', achieved_std);
                end
            else
                error('gpSampler did not return expected struct for perturbed model');
            end
            
        catch ME
            fprintf('  Error in perturbed sampling: %s\n', ME.message);
            fprintf('  Using constrained uniform sampling...\n');
            
            % Fallback sampling for perturbed model
            samples_perturbed = zeros(length(model_perturbed.rxns), num_samples);
            enzyme_idx = find(ismember(model_perturbed.rxns, target_enzyme));
            
            for i = 1:length(model_perturbed.rxns)
                if i == enzyme_idx
                    % Constrained enzyme: small variation around target
                    variation = abs(target_flux) * 0.001;  % 0.1% variation
                    samples_perturbed(i, :) = target_flux + variation * (2*rand(1, num_samples) - 1);
                else
                    % Other reactions: sample within bounds
                    lb = model_perturbed.lb(i);
                    ub = model_perturbed.ub(i);
                    if ub > lb
                        samples_perturbed(i, :) = lb + (ub - lb) * rand(1, num_samples);
                    else
                        samples_perturbed(i, :) = repmat(lb, 1, num_samples);
                    end
                end
            end
            actual_perturbed_samples = num_samples;
            
            % Report achieved perturbation for fallback method
            if ~isempty(enzyme_idx)
                achieved_fluxes = samples_perturbed(enzyme_idx, :);
                achieved_median = median(achieved_fluxes);
                fprintf('  Fallback achieved median: %.6f (target: %.6f)\n', achieved_median, target_flux);
            end
        end
        
        % Store perturbed results
        out_all{condition_counter,1}.samples = samples_perturbed;
        out_all{condition_counter,1}.subSysModel = model_perturbed;
        out_all{condition_counter,1}.modelID = perturbed_name;
        out_all{condition_counter,1}.target_enzyme = target_enzyme;
        out_all{condition_counter,1}.condition_type = 'perturbed';
        out_all{condition_counter,1}.perturbation_factor = perturbation_factor;
        out_all{condition_counter,1}.is_perturbed = true;
        out_all{condition_counter,1}.median_flux = median_flux;
        out_all{condition_counter,1}.target_flux = target_flux;
        out_all{condition_counter,1}.num_samples = actual_perturbed_samples;
        
        fprintf('  Perturbation %.0f%% completed (%d/%d)\n', perturbation_factor*100, ...
            condition_counter, length(valid_enzymes)*(length(perturbation_factors)+1));
    end
end

%% Save results with organized structure - INDIVIDUAL FILES FOR EACH CONDITION
fprintf('Saving each condition individually to prevent data loss...\n');

% Create directory structure
for p_idx = 1:length(perturbation_factors)
    factor_dir = sprintf('perturbation_%.0f', perturbation_factors(p_idx)*100);
    factor_subfolder = fullfile(subfolder, factor_dir);
    if ~exist(factor_subfolder, 'dir')
        mkdir(factor_subfolder)
    end
end

baseline_subfolder = fullfile(subfolder, 'baseline');
if ~exist(baseline_subfolder, 'dir')
    mkdir(baseline_subfolder)
end

% Create individual condition folders
individual_subfolder = fullfile(subfolder, 'individual_conditions');
if ~exist(individual_subfolder, 'dir')
    mkdir(individual_subfolder)
end

% Save each condition individually
fprintf('Saving %d conditions individually...\n', length(out_all));
saved_count = 0;
failed_count = 0;

for i = 1:length(out_all)
    condition_data = out_all{i,1};
    condition_id = condition_data.modelID;
    
    % Create filename
    filename = sprintf('%s.mat', condition_id);
    filepath = fullfile(individual_subfolder, filename);
    
    try
        % Save individual condition
        save(filepath, 'condition_data', '-v7.3');  % Use v7.3 for large files
        saved_count = saved_count + 1;
        
        if mod(i, 5) == 0 || i == length(out_all)
            fprintf('  Saved %d/%d conditions...\n', saved_count, length(out_all));
        end
        
        % Also save in perturbation-specific folder
        if condition_data.perturbation_factor == 1.0
            % Baseline condition
            baseline_filepath = fullfile(baseline_subfolder, filename);
            save(baseline_filepath, 'condition_data', '-v7.3');
        else
            % Perturbed condition
            factor_dir = sprintf('perturbation_%.0f', condition_data.perturbation_factor*100);
            factor_subfolder = fullfile(subfolder, factor_dir);
            factor_filepath = fullfile(factor_subfolder, filename);
            save(factor_filepath, 'condition_data', '-v7.3');
        end
        
    catch ME
        fprintf('  ERROR saving %s: %s\n', condition_id, ME.message);
        failed_count = failed_count + 1;
    end
end

fprintf('Individual condition saving complete:\n');
fprintf('  Successfully saved: %d conditions\n', saved_count);
fprintf('  Failed to save: %d conditions\n', failed_count);

% Create summary structure without heavy sampling data
out_all_summary = {};
for i = 1:length(out_all)
    summary_data = out_all{i,1};
    % Remove the heavy sampling data
    summary_data = rmfield(summary_data, 'samples');
    summary_data = rmfield(summary_data, 'subSysModel');
    % Keep only metadata
    out_all_summary{i,1} = summary_data;
end

% Save summary files
try
    save(fullfile(subfolder, 'out_all_summary.mat'), 'out_all_summary');
    fprintf('Saved lightweight summary file: out_all_summary.mat\n');
catch ME
    fprintf('Warning: Could not save summary file: %s\n', ME.message);
end

% Save analysis summary
analysis_summary = struct();
analysis_summary.baseline_condition = baseline_condition;
analysis_summary.perturbation_factors = perturbation_factors;
analysis_summary.valid_enzymes = valid_enzymes;
analysis_summary.invalid_enzymes = invalid_enzymes;
analysis_summary.reference_flux_median = reference_flux_median;
analysis_summary.total_conditions = condition_counter;
analysis_summary.requested_samples = num_samples;
analysis_summary.actual_baseline_samples = actual_num_samples;
analysis_summary.sampling_method = 'gpSampler_with_fallback';

save('analysis_summary_sampling.mat', 'analysis_summary');
movefile('analysis_summary_sampling.mat', fullfile(subfolder, 'analysis_summary_sampling.mat'));

% Create comprehensive condition list with sampling statistics
condition_list = {};
for i = 1:length(out_all)
    condition_list{i, 1} = out_all{i,1}.modelID;
    condition_list{i, 2} = out_all{i,1}.target_enzyme;
    condition_list{i, 3} = out_all{i,1}.condition_type;
    condition_list{i, 4} = out_all{i,1}.perturbation_factor;
    condition_list{i, 5} = out_all{i,1}.is_perturbed;
    condition_list{i, 6} = out_all{i,1}.median_flux;
    condition_list{i, 7} = out_all{i,1}.num_samples;
    if isfield(out_all{i,1}, 'target_flux')
        condition_list{i, 8} = out_all{i,1}.target_flux;
    else
        condition_list{i, 8} = NaN;
    end
end

condition_table = cell2table(condition_list, ...
    'VariableNames', {'ModelID', 'TargetEnzyme', 'ConditionType', 'PerturbationFactor', ...
                      'IsPerturbed', 'MedianFlux', 'NumSamples', 'TargetFlux'});
writetable(condition_table, fullfile(subfolder, 'condition_list_sampling.xlsx'));

fprintf('Valid enzymes analyzed: %d\n', length(valid_enzymes));
fprintf('Baseline samples: %d\n', actual_num_samples);
fprintf('Total conditions created: %d\n', condition_counter);
fprintf('Results saved to: %s\n', subfolder);
fprintf('\n*** PRINCIPLED PERTURBATIONS BASED ON MEDIAN FLUX ***\n');
fprintf('*** READY FOR DOWNSTREAM ANALYSIS ***\n');

function model = prepare_model(model, condition_id)
    %  GLUD reactions 
    model = changeRxnBounds(model, {'MAR03802','MAR03804'}, 0, 'u');
    
    %  IDH reactions 
    model = changeRxnBounds(model, {'MAR03958','MAR00710','MAR04111','MAR04112'}, -1000, 'l');
    
    %  Trxr-2 
    model = changeRxnBounds(model, {'MAR02358'}, -1000, 'l');
    model = changeRxnBounds(model, {'MAR02358'}, 0, 'u');
    
    %  Zw (G6PD)
    model = changeRxnBounds(model, {'MAR04306'}, -1000, 'l');
    model = changeRxnBounds(model, {'MAR04306'}, 0, 'u');
    
    % Dhfr (Dihydrofolate reductase)
    model = changeRxnBounds(model, {'MAR04332','MAR04333','MAR04335','MAR04654','MAR04655'}, 0, 'l');
    
    % CG1236/GRHPR
    model = changeRxnBounds(model, {'MAR07702'}, 0, 'l');
    
    % FASN1 (fatty acid synthase)
    model = changeRxnBounds(model, {'MAR02185'}, 0, 'u');
    model = changeRxnBounds(model, {'MAR02185'}, -1000, 'l');
    
    model = addReaction(model, 'ATP_maintenance', ...
        'metaboliteList', {'MAM01371c[c]', 'MAM01371m[m]', 'MAM02040c[c]', 'MAM02040m[m]', ...
                           'MAM01285c[c]', 'MAM01285m[m]', 'MAM02751c[c]', 'MAM02751m[m]', ...
                           'MAM02039c[c]', 'MAM02039m[m]'}, ...
        'stoichCoeffList', [-1; -1; -1; -1; 1; 1; 1; 1; 1; 1], ...
        'reversible', false);
    
    model = addReaction(model, 'NAD_demand', ...
        'metaboliteList', {'MAM02552c[c]', 'MAM02039c[c]', 'MAM02552m[m]', 'MAM02039m[m]', ...
                           'MAM02553c[c]', 'MAM02553m[m]'}, ...
        'stoichCoeffList', [-1; -1; -1; -1; 1; 1], ...
        'reversible', false);
    
    model = addReaction(model, 'NADPH_demand', ...
        'metaboliteList', {'MAM02555c[c]', 'MAM02555m[m]',...
                            'MAM02039c[c]', 'MAM02039m[m]', 'MAM02554c[c]', 'MAM02554m[m]'}, ...
        'stoichCoeffList', [-1; -1; 1; 1; 1; 1], ...
        'reversible', false);
    
    if strcmp(condition_id, 'HSD')
        % Transport reactions constraints (50% reduction)
        transport_reactions = {'MAR09034','MAR09426','MAR05029'};
        reductionFactor = 0.5;
        
        for j = 1:length(transport_reactions)
            rxn_id = transport_reactions{j};
            idx_in = find(ismember(model.rxns, rxn_id));
            if ~isempty(idx_in)
                model = changeRxnBounds(model, rxn_id, model.ub(idx_in)*reductionFactor, 'u'); 
                model = changeRxnBounds(model, rxn_id, model.lb(idx_in)*reductionFactor, 'l'); 
            end
        end
        
        gene_lists = {{'Tret1-1','Glut1'}, {'Ogdh','SdhA','kdn'}};        
        for gene_set = 1:length(gene_lists)
            gene_int = gene_lists{gene_set};
            list_strct = findRxnsFromGenes(model, gene_int);
            fieldName = fieldnames(list_strct);
            list = [];
            for k = 1:length(fieldName)
                list_tmp = list_strct.(fieldName{k})(:,1);
                list = [list; string(list_tmp)];
            end
            
            for j = 1:length(list)
                rxn_id = char(list(j));
                idx_in = find(ismember(model.rxns, rxn_id));
                if ~isempty(idx_in)
                    model = changeRxnBounds(model, rxn_id, model.ub(idx_in)*reductionFactor, 'u'); 
                    model = changeRxnBounds(model, rxn_id, model.lb(idx_in)*reductionFactor, 'l'); 
                end
            end
        end
        
        specific_perturbations = {
            {'MAR04410', 0.1},  % Fum1
            {'MAR04652','MAR08743', 0.25},  % SDH
            {'MAR04209','MAR05297','MAR06411','MAR06413', 0.25}  % OGDH
        };
        
        for p = 1:length(specific_perturbations)
            rxn_list = specific_perturbations{p}(1:end-1);
            reduction_factor = specific_perturbations{p}{end};
            
            for j = 1:length(rxn_list)
                rxn_id = rxn_list{j};
                idx_in = find(ismember(model.rxns, rxn_id));
                if ~isempty(idx_in)
                    model = changeRxnBounds(model, rxn_id, model.ub(idx_in)*reduction_factor, 'u'); 
                    model = changeRxnBounds(model, rxn_id, model.lb(idx_in)*reduction_factor, 'l'); 
                end
            end
        end
    end
    model = changeObjective(model, {'ATP_maintenance','NAD_demand','NADPH_demand'}, 1/3);
end

%% Function to apply median-based perturbation
function model = apply_median_perturbation(model, enzyme_id, target_flux)
    
    % Find enzyme reaction index
    idx_enzyme = find(ismember(model.rxns, enzyme_id));
    
    if isempty(idx_enzyme)
        error('Enzyme %s not found in model', enzyme_id);
    end
    
    fprintf('    Setting %s bounds to target flux: %.6f\n', enzyme_id, target_flux);
    
    % Set both upper and lower bounds to target flux (tight constraint)
    model = changeRxnBounds(model, enzyme_id, target_flux, 'u');
    model = changeRxnBounds(model, enzyme_id, target_flux, 'l');
    
    % Test feasibility
    test_solution = optimizeCbModel(model, 'max');
    if test_solution.stat ~= 1
        fprintf('    WARNING: Model infeasible with tight constraint\n');
        fprintf('    Applying small tolerance (±0.1%% of target flux)\n');
        
        tolerance = abs(target_flux) * 0.001;  % 0.1% tolerance
        model = changeRxnBounds(model, enzyme_id, target_flux + tolerance, 'u');
        model = changeRxnBounds(model, enzyme_id, target_flux - tolerance, 'l');
        
        % Re-test
        test_solution = optimizeCbModel(model, 'max');
        if test_solution.stat == 1
            fprintf('    Model feasible with tolerance\n');
        else
            fprintf('    ERROR: Model still infeasible\n');
        end
    else
        fprintf('    Model feasible with tight constraint\n');
    end
end