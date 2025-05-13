%% Goal:
% - Check normality of NSD vs. HSD conditions
% - Compute skewness + test skewness significance
% - Save results to Excel

close all;
clc;

%% Define paths
pathway = pwd;
save_dirFlux = 'I_normality';
subfolder = fullfile(pathway, save_dirFlux);
if ~exist(subfolder, 'dir')
    mkdir(subfolder);
end

%% Load sampling data
model_sampling_in_name = {'NSD','HSD'};
model_sampling_in = cell(1, length(model_sampling_in_name));

for i = 1:length(model_sampling_in_name)
   fileName = fullfile(pathway, 'C_sampling', [model_sampling_in_name{i}, '.mat']);
   model_sampling_in{i} = load(fileName);
end

%% Extract reactions and points
rxns = model_sampling_in{1}.samples.rxns;
all_rxns_nsd = model_sampling_in{1}.samples.points;
all_rxns_hsd = model_sampling_in{2}.samples.points;

num_rxns = length(rxns);
T_results = table('Size', [num_rxns 8], ...
    'VariableTypes', {'string', 'double', 'double', 'double', 'double', 'double', 'double', 'string'}, ...
    'VariableNames', {'Reaction', 'p_NSD', 'p_HSD', 'Skewness_NSD', 'Skewness_HSD', 'p_Skew_NSD', 'p_Skew_HSD', 'Flag'});

T_results.Reaction = rxns;

%% Loop through reactions and compute normality + skewness + skewness significance
for i = 1:num_rxns
    flux_nsd = all_rxns_nsd(i, :);
    flux_hsd = all_rxns_hsd(i, :);

    % Sample sizes
    n_nsd = sum(~isnan(flux_nsd));
    n_hsd = sum(~isnan(flux_hsd));

    % Skewness
    sk_nsd = skewness(flux_nsd, 0, 2);
    sk_hsd = skewness(flux_hsd, 0, 2);

    % Skewness significance (z-test)
    z_skew_nsd = sk_nsd / sqrt(6 / n_nsd);
    z_skew_hsd = sk_hsd / sqrt(6 / n_hsd);
    p_skew_nsd = 2 * (1 - normcdf(abs(z_skew_nsd)));
    p_skew_hsd = 2 * (1 - normcdf(abs(z_skew_hsd)));

    % Normality test (Lilliefors)
    [H_nsd, p_nsd] = lillietest(flux_nsd);
    [H_hsd, p_hsd] = lillietest(flux_hsd);

    % Store results
    T_results.p_NSD(i) = p_nsd;
    T_results.p_HSD(i) = p_hsd;
    T_results.Skewness_NSD(i) = sk_nsd;
    T_results.Skewness_HSD(i) = sk_hsd;
    T_results.p_Skew_NSD(i) = p_skew_nsd;
    T_results.p_Skew_HSD(i) = p_skew_hsd;

    % Flag logic
    isNonNormal = H_nsd == 1 || H_hsd == 1;
    isSkewed = p_skew_nsd < 0.05 || p_skew_hsd < 0.05;

    if isNonNormal && isSkewed
        T_results.Flag(i) = "Non-normal & Skewed";
    elseif isNonNormal
        T_results.Flag(i) = "Non-normal";
    elseif isSkewed
        T_results.Flag(i) = "Skewed";
    else
        T_results.Flag(i) = "Normal";
    end
end

%% Save detailed results
writetable(T_results, fullfile(subfolder, 'Normality_Skewness_Results_AllRxns.xlsx'));
fprintf("✅ Normality + Skewness check complete. Results saved to Excel.\n");

%% ------------------ Summary Statistics ------------------
num_normal = sum(T_results.Flag == "Normal");
num_skewed = sum(T_results.Flag == "Skewed");
num_nonnormal = sum(T_results.Flag == "Non-normal");
num_both = sum(T_results.Flag == "Non-normal & Skewed");

avg_skew_nsd = mean(T_results.Skewness_NSD, 'omitnan');
avg_skew_hsd = mean(T_results.Skewness_HSD, 'omitnan');
std_skew_nsd = std(T_results.Skewness_NSD, 'omitnan');
std_skew_hsd = std(T_results.Skewness_HSD, 'omitnan');

summary_table = table;
summary_table.Total_Reactions = num_rxns;
summary_table.Num_Normal = num_normal;
summary_table.Num_Skewed = num_skewed;
summary_table.Num_NonNormal = num_nonnormal;
summary_table.Num_NonNormal_And_Skewed = num_both;
summary_table.Perc_Normal = round(100 * num_normal / num_rxns, 2);
summary_table.Avg_Skewness_NSD = avg_skew_nsd;
summary_table.Std_Skewness_NSD = std_skew_nsd;
summary_table.Avg_Skewness_HSD = avg_skew_hsd;
summary_table.Std_Skewness_HSD = std_skew_hsd;

%% Save summary
writetable(summary_table, fullfile(subfolder, 'Normality_Skewness_Summary.xlsx'));
fprintf("Summary table saved to Excel.\n");
