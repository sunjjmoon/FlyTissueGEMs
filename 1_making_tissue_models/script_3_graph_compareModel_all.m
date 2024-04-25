% Goal:
%   - To compare the models of all the tissues.

% load(model_all.mat);
% Output:   1. out1 - clusterogram. 
%           2. out2 - clusterogram with different structures
%               After making the output, use plot(xxx) and modify the
%               figures there.
%           3. tissue_info:
%               - This will have the updated tissue annotations (e.g., no _
%               but space, etc).

load(strcat(pwd,'\models\','model_all.mat'));

save_dir = '3';
tsneGraph = 1; % 1 is to graph
% Difference by the mean - Original clusterogram
subSysThreshold = 100; % This is to check how much of the difference 
subSysThreshold_v2 = 20; % ;absolute difference
subSysThreshold_v3 = 100; % el difference (%)

control_col = 2; % use the 2nd column (tissue) as the control.

run_compareModel_all(model_all,save_dir,tsneGraph,...
subSysThreshold,subSysThreshold_v2,subSysThreshold_v3,control_col)

function run_compareModel_all(model_all,tsneGraph,subSysThreshold,subSysThreshold_v2,subSysThreshold_v3,control_col)
close all

close all

res = compareMultipleModels(model_all);

useModels = res.modelIDs;

if tsneGraph == 1

    %% Clustering 
    dendo = clustergram(res.structComp, 'Symmetric', false, 'Colormap', 'bone', 'RowLabels', res.modelIDs, 'ColumnLabels', res.modelIDs);
    cg2 = clustergram(res.structComp, 'Symmetric', false, 'Colormap', redbluecmap, 'RowLabels', res.modelIDs, 'ColumnLabels', res.modelIDs);

%% tSNE 
    figure(10)
    rxn2Dmap = tsne(res.reactions.matrix', 'Distance', 'hamming', 'NumDimensions', 2, 'Perplexity', 3);
    scatter(rxn2Dmap(:,1), rxn2Dmap(:,2)); 
    text(rxn2Dmap(:,1), rxn2Dmap(:,2), res.modelIDs);
    xlabel('tSNE1')
    ylabel('tSNE2')
    hold off
end

%% Subsystem Coverage
%% This is to compare the number of reactions in each model
keep = ismember(res.modelIDs, useModels);
subMat = res.subsystems.matrix(:,keep);
subCoverage = (subMat-mean(subMat,2))./mean(subMat,2)*100;
subCoverage_v2 = subMat-subMat(:,control_col); % absolute difference
subCoverage_v3 = (subMat-subMat(:,control_col))./subMat(:,control_col)*100; %rel difference (%)
%%
% Visualize the difference in subsystem coverage with a clustergram,
% including only subsystems with at least a xx% difference in one or more
% GEMs.
% select subsystems to include in plot
inclSub = any(abs(subCoverage)>subSysThreshold,2); % return the logical columns.

subNames = res.subsystems.ID(inclSub);
inclSub_v2 = any(abs(subCoverage_v2)>subSysThreshold_v2,2); % return the logical columns.
subNames_v2 = res.subsystems.ID(inclSub_v2);

% subNames = useModels;
inclSub_v3 = any(abs(subCoverage_v3)>subSysThreshold_v3,2); % return the logical columns.
subNames_v3 = res.subsystems.ID(inclSub_v3);

out1 = clustergram(subCoverage(inclSub,:), 'Colormap', redbluecmap, 'DisplayRange', 100, 'rowLabels', subNames, 'columnLabels', useModels, 'ShowDendrogram', 'OFF');
end
