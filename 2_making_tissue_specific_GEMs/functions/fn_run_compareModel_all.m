function fn_run_compareModel_all(model_all, save_dir, tsneGraph, subSysThreshold, subSysThreshold_v2, subSysThreshold_v3, control_col)
% RUN_COMPAREMODEL_ALL
% Compares tissue-specific GEMs using structural similarity, t-SNE, and subsystem coverage.
% 
% Inputs:
%   model_all            - Cell array of models to compare
%   save_dir             - Output directory for results
%   tsneGraph            - Boolean flag (1 = generate t-SNE plot)
%   subSysThreshold      - Threshold (%) for subsystem variation (clustergram 1)
%   subSysThreshold_v2   - Threshold (absolute) for reaction count difference
%   subSysThreshold_v3   - Threshold (%) for relative difference from control
%   control_col          - Column index for control tissue (e.g., 2 = Muscle)
%
% Outputs:
%   Saved figures, cluster assignments, and subsystem difference tables in save_dir

%% Set up
pathway = pwd;
subfolder = fullfile(pathway, save_dir);
if ~exist(subfolder, 'dir')
    mkdir(subfolder)
end
close all;

%% Compare Models
res = compareMultipleModels(model_all);
useModels = res.modelIDs;

%% t-SNE and Clustering (if enabled)
if tsneGraph == 1
    % Hierarchical clustering (raw and colored)
    clustergram(res.structComp, 'Symmetric', false, 'Colormap', 'bone', 'RowLabels', res.modelIDs, 'ColumnLabels', res.modelIDs);
    clustergram(res.structComp, 'Symmetric', false, 'Colormap', redbluecmap, 'RowLabels', res.modelIDs, 'ColumnLabels', res.modelIDs);

    % t-SNE group assignments
    group_names = ["Muscle", "Fat body/Oenocyte", "Gut", "Neurons", "Glia", "Rest"];
    group_colors = [
        1.00, 0.00, 0.00;  % Red
        1.00, 0.00, 1.00;  % Magenta
        1.00, 0.41, 0.16;  % Orange
        0.64, 0.08, 0.18;  % Deep red
        0.47, 0.67, 0.19;  % Green
        0.80, 0.80, 0.80   % Gray
    ];
    shapes = ['o', 's', 'd', '^', 'v', 'h'];

    % Assign groups manually
    results.Group = repmat("Rest", length(res.modelIDs), 1);
    results.Group([26, 21]) = "Muscle";
    results.Group([2, 5]) = "Fat body/Oenocyte";
    results.Group([4, 11]) = "Gut";
    results.Group([6, 9, 19, 22, 23, 25, 30]) = "Neurons";
    results.Group([1, 3, 7, 10, 27, 32]) = "Glia";

    % t-SNE embedding
    rxn2Dmap = tsne(res.reactions.matrix', 'Distance', 'hamming', 'NumDimensions', 2, 'Perplexity', 3);
    figure;
    hold on;
    for i = 1:length(group_names)
        idx = results.Group == group_names(i);
        scatter(rxn2Dmap(idx, 1), rxn2Dmap(idx, 2), 100, group_colors(i, :), shapes(i), 'filled');
    end
    for i = 1:length(res.modelIDs)
        text(rxn2Dmap(i, 1), rxn2Dmap(i, 2), res.modelIDs{i}, 'FontSize', 8, ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');
    end
    xlabel('tSNE1'); ylabel('tSNE2');
    title('t-SNE of Tissue GEMs');
    set(gca, 'FontSize', 12);
    saveas(gcf, fullfile(subfolder, 'tSNE_tissue_clustering.png'));
    saveas(gcf, fullfile(subfolder, 'tSNE_tissue_clustering.svg'));
    disp('t-SNE plot generated and saved.');
end

%% Hierarchical Clustering of Model Structures
distanceMatrix = pdist(res.structComp, 'euclidean');
linkageMatrix = linkage(distanceMatrix, 'Average');
figure;
dendrogram(linkageMatrix, 0, 'Labels', res.modelIDs, 'Orientation', 'top');
xtickangle(45);
title('Hierarchical Clustering of Tissue GEMs');
ylabel('Euclidean Distance');
set(gca, 'FontSize', 12, 'FontWeight', 'bold');
saveas(gcf, fullfile(subfolder, 'tissue_kmeans_clustering.png'));
savefig(fullfile(subfolder, 'tissue_kmeans_clustering.fig'));

%% Assign Clusters (K-means style)
numClusters = 9;
clusterLabels = cluster(linkageMatrix, 'maxclust', numClusters);
clusterTable = table(res.modelIDs, clusterLabels, 'VariableNames', {'Tissue', 'Cluster'});
writetable(clusterTable, fullfile(subfolder, ['tissue_clusters_k' num2str(numClusters) '.xlsx']));
disp('Cluster assignments saved.');

%% Subsystem Coverage Analysis
keep = ismember(res.modelIDs, useModels);
subMat = res.subsystems.matrix(:, keep);

% Calculate differences
subCoverage = (subMat - mean(subMat, 2)) ./ mean(subMat, 2) * 100;
subCoverage_v2 = subMat - subMat(:, control_col);
subCoverage_v3 = (subMat - subMat(:, control_col)) ./ subMat(:, control_col) * 100;

% Identify varying subsystems
inclSub = any(abs(subCoverage) > subSysThreshold, 2);
inclSub_v2 = any(abs(subCoverage_v2) > subSysThreshold_v2, 2);
inclSub_v3 = any(abs(subCoverage_v3) > subSysThreshold_v3, 2);

% Subsystem names
subNames = res.subsystems.ID(inclSub);
subNames_v2 = res.subsystems.ID(inclSub_v2);
subNames_v3 = res.subsystems.ID(inclSub_v3);

% Clustergram 1: % difference
out1 = clustergram(subCoverage(inclSub,:), 'Colormap', redbluecmap, ...
    'DisplayRange', 100, 'RowLabels', subNames, 'ColumnLabels', useModels, 'ShowDendrogram', 'off');
out1Data = array2table(subCoverage(inclSub,:), ...
    'RowNames', subNames, 'VariableNames', useModels);
writetable(out1Data, fullfile(subfolder, 'out1_subsystem_data.csv'), 'WriteRowNames', true);
disp('Subsystem % difference table saved.');

% Heatmap: absolute difference
figure(22);
h = heatmap(useModels, subNames_v2, subCoverage_v2(inclSub_v2,:), ...
    'Colormap', redbluecmap, 'ColorbarVisible', 'on', 'GridVisible', 'on');
s = struct(h); s.XAxis.TickLabelRotation = 90;
maxVal = max(h.ColorLimits); h.ColorLimits = [-maxVal, maxVal];
sorty(h, useModels{control_col}, 'descend');
title('Absolute difference in number of reactions');
print(gcf, fullfile(subfolder, 'rxn.png'));

% Heatmap: % difference
figure(4);
h3 = heatmap(useModels, subNames_v3, subCoverage_v3(inclSub_v3,:), ...
    'Colormap', redbluecmap, 'ColorbarVisible', 'on', 'GridVisible', 'on');
s3 = struct(h3); s3.XAxis.TickLabelRotation = 90;
maxVal = max(h3.ColorLimits); h3.ColorLimits = [-maxVal, maxVal];
sorty(h3, useModels{control_col}, 'descend');
title('% Difference in Subsystem Coverage');
print(gcf, fullfile(subfolder, 'all.png'));

% Heatmap (no cell labels)
figure(44);
h3 = heatmap(useModels, subNames_v3, subCoverage_v3(inclSub_v3,:), ...
    'Colormap', redbluecmap, 'CellLabelColor', 'none', ...
    'ColorbarVisible', 'on', 'GridVisible', 'on');
s3 = struct(h3); s3.XAxis.TickLabelRotation = 90;
h3.ColorLimits = [-maxVal, maxVal];
sorty(h3, useModels{control_col}, 'descend');
title('% Difference in Subsystem Coverage');
print(gcf, fullfile(subfolder, 'all_noLabel.png'));

end
