%% Load flagged results from Excel
T_results = readtable(fullfile(pwd, 'I_normality', 'Normality_Skewness_Results_AllRxns.xlsx'));

%% Count flags
flag_all = T_results.Flag;

num_total = height(T_results);
num_normal = sum(flag_all == "Normal");
num_skewed = sum(flag_all == "Skewed");
num_nonnormal = sum(flag_all == "Non-normal");
num_both = sum(flag_all == "Non-normal & Skewed");

%% Compute fractions
frac_normal = num_normal / num_total;
frac_skewed = num_skewed / num_total;
frac_nonnormal = num_nonnormal / num_total;
frac_both = num_both / num_total;

%% Display results
fprintf("📊 Flag Breakdown (%d total reactions):\n", num_total);
fprintf("✔️  Normal only           : %d (%.2f%%)\n", num_normal, 100*frac_normal);
fprintf("↩️  Skewed only           : %d (%.2f%%)\n", num_skewed, 100*frac_skewed);
fprintf("🚫  Non-normal only       : %d (%.2f%%)\n", num_nonnormal, 100*frac_nonnormal);
fprintf("💥  Non-normal & Skewed  : %d (%.2f%%)\n", num_both, 100*frac_both);

%% Extract skewness values
skew_nsd = T_results.Skewness_NSD;
skew_hsd = T_results.Skewness_HSD;

%% Plot histogram of skewness (NSD vs HSD)
figure;
histogram(skew_nsd, 50, 'FaceAlpha', 0.6, 'EdgeColor', 'none','Normalization','probability'); hold on;
histogram(skew_hsd, 50, 'FaceAlpha', 0.6, 'EdgeColor', 'none','Normalization','probability');
legend({'NSD', 'HSD'}, 'Location', 'northeast');
xlabel('Skewness');
ylabel('Probability');
title('Distribution of Skewness across Reactions');
grid on;

pathway = pwd;
save_dirFlux = 'J_skewness';
subfolder = [pathway '\' save_dirFlux];
if ~exist(subfolder, 'dir')
    mkdir(subfolder)
end

saveas(gcf, char(strcat(subfolder,'\skeness.svg')));
saveas(gcf, char(strcat(subfolder,'\skeness.png')));

%% Statistical test: Paired comparison of skewness between NSD and HSD
% Remove NaNs if any
valid_idx = ~isnan(skew_nsd) & ~isnan(skew_hsd);
skew_nsd_clean = skew_nsd(valid_idx);
skew_hsd_clean = skew_hsd(valid_idx);

[p_skew, h_skew, stats] = signrank(skew_nsd_clean, skew_hsd_clean);

fprintf("\n📈 Wilcoxon signed-rank test for skewness (NSD vs HSD):\n");
fprintf("p-value = %.4g\n", p_skew);
if h_skew
    fprintf("Result: Significant difference in skewness distributions ✅\n");
else
    fprintf("Result: No significant difference in skewness distributions ❌\n");
end

diff_skew = skew_hsd_clean - skew_nsd_clean;
median_diff = median(diff_skew);
fprintf("Median skewness shift (HSD - NSD): %.4f\n", median_diff);
