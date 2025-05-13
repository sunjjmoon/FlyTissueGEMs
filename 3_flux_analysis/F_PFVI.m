%% =================== PFI Calculation Script (Updated) ===================
% Purpose: Compute Pathway Flux Index (PFI) and its statistics (Mean, SD, SEM) per subsystem
% ===============================================================

clear all; 
close all; 
clc;

%% ------------------ Set Folder Paths ------------------
% Define pathway and load folder
pathway = pwd;
loadfolder = fullfile(pathway, 'B_fva');

% Load FVA Data (Bounds for Reactions)
if ~exist('out_all_fvaBounded.mat','var')
    load(fullfile(loadfolder, 'out_all_fvaBounded.mat'));
end

% Load Sampling Data for NSD & HSD Conditions
model_sampling_in_name = {'NSD', 'HSD'};
model_sampling_in = cell(1, length(model_sampling_in_name));

for i = 1:length(model_sampling_in_name)
    fileName = fullfile(pathway, 'C_sampling', model_sampling_in_name{i});
    model_sampling_in{i} = load(fileName);
end

% Load Subsystem Definitions
if ~exist('subsysAll_fruitflyGEM.mat','var')
    load(fullfile(pathway, 'subsysAll_fruitflyGEM.mat'));
end

%% ------------------ Create Output Directory ------------------
subfolder = fullfile(pathway, 'F_PFVI_mode');
if ~exist(subfolder, 'dir')
    mkdir(subfolder);
end

%% ------------------ Compute PFI with Statistics ------------------
[modelID, PFI_mean, PFI_SD, PFI_SEM, PFI_COUNT] = compute_pfi_stats(out_all_fvaBounded, subfolder, subsysAll_fruitflyGEM, model_sampling_in);

%% =================== PFI Calculation Function (with Statistics) ===================
function [modelID, pfi_table, pfi_table_sd, pfi_table_sem,pfi_table_count] = compute_pfi_stats(out_all_fvaBounded, subfolder, subsysAll_fruitflyGEM, model_sampling_in)
    close all;

    % Initialize storage for results
    PFI_matrix = zeros(length(subsysAll_fruitflyGEM), length(out_all_fvaBounded));
    PFI_SD_matrix = zeros(length(subsysAll_fruitflyGEM), length(out_all_fvaBounded));
    PFI_SEM_matrix = zeros(length(subsysAll_fruitflyGEM), length(out_all_fvaBounded));
    count_matrix = zeros(length(subsysAll_fruitflyGEM), length(out_all_fvaBounded));

    % Iterate through all defined subsystems
    for j = 1:length(subsysAll_fruitflyGEM)
        % Iterate through each model instance (e.g., NSD, HSD)
        for i = 1:length(out_all_fvaBounded)
            % Extract model FVA bounds
            model = out_all_fvaBounded{i, 1};
            modelID{i, 1} = model.modelID;
            subsys_model = model.subSystems;

            % Standardize subsystem names
            for k = 1:length(subsys_model)
                if iscell(subsys_model{k, 1})
                    subsys_model{k, 1} = subsys_model{k, 1}{1, 1};
                end
            end
            subsys_model = string(subsys_model);
            
            % Identify reactions belonging to the current subsystem
            idxTmp = find(ismember(subsys_model, subsysAll_fruitflyGEM(j)));
            
            % Proceed only if reactions exist in the subsystem
            if ~isempty(idxTmp)
                % Extract reaction names
                rxnNames = model.rxns(idxTmp);
                
                % Extract flux variability bounds (min/max values)
                lb_temp = model.lb(idxTmp);
                ub_temp = model.ub(idxTmp);
                
                % Extract mean flux values from sampling
                mean_flux_values = zeros(length(rxnNames), 1);
                for p = 1:length(rxnNames)
                    model_tmp = model_sampling_in{1, i}.samples;
                    idxFlux = find(ismember(model_tmp.rxns, rxnNames(p)));
                    rounded_flux = round(model_tmp.points(idxFlux, :), 3);  % round to 3 decimals
                    [unique_vals, ~, idx_u] = unique(rounded_flux, 'stable');  % preserve order
                    counts = accumarray(idx_u, 1);  % count frequencies
                    max_count = max(counts);
                    modes = unique_vals(counts == max_count);  % all tied modes
                    mean_flux_values(p, 1) = modes(1);  % pick the first one         
                end
                
                % Compute PFI values per reaction
                individual_pfi_values = abs(ub_temp - lb_temp) .* abs(mean_flux_values);
                
                % Compute PFI statistics
                PFI_matrix(j, i) = mean(individual_pfi_values, 'omitnan');  % Mean PFI
                PFI_SD_matrix(j, i) = std(individual_pfi_values, 0, 'omitnan');  % Standard Deviation
                PFI_SEM_matrix(j, i) = PFI_SD_matrix(j, i) / sqrt(length(idxTmp));  % Standard Error of Mean
                count_matrix(j,i) = length(idxTmp);
            else
                % Assign zero when no reactions exist
                PFI_matrix(j, i) = 0;
                PFI_SD_matrix(j, i) = 0;
                PFI_SEM_matrix(j, i) = 0;
                count_matrix(j,i) = 0;
            end
        end
        disp(subsysAll_fruitflyGEM(j))
    end
    
    % Save PFI values to Excel
    pfi_table = array2table(PFI_matrix, 'RowNames', subsysAll_fruitflyGEM,'VariableNames',modelID);
    writetable(pfi_table, fullfile(subfolder, 'PFVI_results.xlsx'),'WriteRowNames',true);

    pfi_table_sd = array2table(PFI_SD_matrix, 'RowNames', subsysAll_fruitflyGEM,'VariableNames',modelID);
    writetable(pfi_table_sd, fullfile(subfolder, 'PFVI_results_sd.xlsx'),'WriteRowNames',true);

    pfi_table_sem = array2table(PFI_SEM_matrix, 'RowNames', subsysAll_fruitflyGEM,'VariableNames',modelID);
    writetable(pfi_table_sem, fullfile(subfolder, 'PFVI_results_sem.xlsx'),'WriteRowNames',true);

    pfi_table_count = array2table(count_matrix, 'RowNames', subsysAll_fruitflyGEM,'VariableNames',modelID);
    writetable(pfi_table_count, fullfile(subfolder, 'PFVI_results_count.xlsx'),'WriteRowNames',true);

    % Print completion message
    fprintf("✅ PFI Calculation Complete! Results saved to 'PFI_results.xlsx' and 'PFI_results_stats.xlsx'.\n");
end
