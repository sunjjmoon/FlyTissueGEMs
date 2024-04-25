%% Set the save and load folder
pathway = pwd;
loadfolder = [pathway '\' 'B_fva_update_v2'];

% load(strcat(loadfolder,'\','out_all_fvaBounded.mat'));
% load(strcat(pwd,'\files\','subsysAll_fruitflyGEM.mat'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter setting
subfolder = [pathway '\' 'D_2_graph'];
if ~exist(subfolder, 'dir')
    mkdir(subfolder)
end

subSysThreshold = 200; % flux difference is greathan 2 times

modelID = graph_box_max_min_metabolicSubsys(out_all_fvaBounded,subfolder,subsysAll_fruitflyGEM,subSysThreshold);

function modelID = graph_box_max_min_metabolicSubsys(out_all_fvaBounded,subfolder,subsysAll_fruitflyGEM,subSysThreshold)
close all

matrix_row_subSys_col_tissue = [];
for j = 1:length(subsysAll_fruitflyGEM)
%% Extract the reaction of interest lb and ub
    lb = zeros(length(out_all_fvaBounded),1);
    ub = zeros(length(out_all_fvaBounded),1);
    pfv = zeros(length(out_all_fvaBounded),1); % store the range of Vmax-Vmin in the pathway
    
    for i = 1:length(out_all_fvaBounded)
        
        model = out_all_fvaBounded{i,1};
        modelID{i,1} = model.modelID;
        subsys_model = model.subSystems;
    
        for k = 1:length(subsys_model)
            if numel(subsys_model{k,1}) ~= 1  % if there is more than 1 cell in the subsystem, use the first one
               subsys_model{k,1} = subsys_model{k,1}{1,1};
            end
        end
        subsys_model = string(subsys_model);
        [Lia,~] = ismember(subsys_model,subsysAll_fruitflyGEM(j));
        idxTmp_orig = find(Lia==1);
        idxTmp =find(ismember(subsys_model,subsysAll_fruitflyGEM(j)));
        
       if isempty(idxTmp) ~= 1 % If the reaction exists
            lb_temp = model.lb(idxTmp);
            ub_temp = model.ub(idxTmp); 
            pfv(i,1) = length(idxTmp).*sum(ub_temp-lb_temp);% Find the median of the range    


       else
            pfv(i,1) = 0;% Find the median            
       end 
    end
matrix_row_subSys_col_tissue{j,1} = pfv';

end
matrix_row_subSys_col_tissue = cell2mat(matrix_row_subSys_col_tissue);
prep = zscore(matrix_row_subSys_col_tissue,1,2); % Standardize by row
RowLabels = subsysAll_fruitflyGEM;

t_tmp = table(RowLabels,matrix_row_subSys_col_tissue);
t_tmp = splitvars(t_tmp);
t_tmp.Properties.VariableNames = ['subsys';modelID];
writetable(t_tmp,strcat(subfolder,'\','tmp.xlsx'),'WriteVariableNames',true);
%% Make the cluster graph 
% Specify a custom colormap
rng(1); % for reproducibility
Y = tsne(prep');
Z = linkage(prep','ward');
dendrogram(Z)
T = cluster(Z, 'criterion', 'inconsistent', 'cutoff', 1);
figure
h = gscatter(Y(:,1), Y(:,2), T);
set(h, 'MarkerSize', 30); % Adjust the size (50 in this example)
hold on;
text(Y(:, 1), Y(:, 2), modelID, 'FontSize', 9, 'Color', 'k', 'HorizontalAlignment', 'right');
title('Comparison of pathway flux variability');
xlabel('tSNE1');
ylabel('tSNE2');
tmptmp = table(modelID,T);
% Add new legen
legend('Location', 'eastoutside');
legend boxoff
writetable(tmptmp,strcat(subfolder,'\','model_cluster.xlsx'));

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);

%% tsne for reactions
rng(1); % for reproducibility
Y = tsne(prep);
Z = linkage(prep,'ward');
dendrogram(Z)
% Use the inconsistency criterion to automatically determine the number of clusters
T = cluster(Z, 'criterion', 'inconsistent', 'cutoff', 1.1);
% T = clusterdata(Z,'Linkage','average','Maxclust',4);
figure
% scatter(Y(:,1), Y(:,2), 50, T, 'filled'); % Adjust the size (50 in this example)
h = gscatter(Y(:,1), Y(:,2), T);
set(h, 'MarkerSize', 30); % Adjust the size (50 in this example)
% Add confidence ellipses around cluster centroids
hold on;
text(Y(:, 1), Y(:, 2), RowLabels, 'FontSize', 10, 'Color', 'k', 'HorizontalAlignment', 'right');
title('Comparison of metabolic flux variability');
xlabel('tSNE1');
ylabel('tSNE2');
tmptmp = table(RowLabels,T);
legend('Location', 'eastoutside');
legend boxoff

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
% print -dpdf -painters tsne_FVA

% Make sure to close the clusterogram
RowLabels = subsysAll_fruitflyGEM;
cg = clustergram(matrix_row_subSys_col_tissue,...
        'RowLabels',RowLabels,...
        'ColumnLabels',modelID,...
        'Colormap',jet,...
        'Standardize','Row',...
        'Cluster','Row',...
        'ColumnPdist','correlation',...
        'ShowDendrogram','off',... 
        'ImputeFun', @knnimpute);  % DendogramShooff- I first said 'on' and then investigated manually to find the groups
RowLabels = subsysAll_fruitflyGEM;
cg = clustergram(matrix_row_subSys_col_tissue,...
        'RowLabels',RowLabels,...
        'ColumnLabels',modelID,...
        'Colormap',jet,...
        'Standardize','Row',...
        'Cluster','Column',...
        'ColumnPdist','correlation',...
        'ShowDendrogram','off',... 
        'ImputeFun', @knnimpute);  % DendogramShooff- I first said 'on' and then investigated manually to find the groups
cg.Dendrogram = 1;
cgax.FontSize = 10;
cgax.XTickLabelRotation = 30;
cg_plot = plot(cg);
cg_plot.FontSize = 12;
cg_plot.Title.String = ['Pathway flux variability (n_{path}\Sigma(V_{i,max}-V_{i,min}))'];

cg_plot.Title.FontSize = 12;
axes(cg_plot);
% cg_plot.Title = 
%% Once finished up to the above, pause, and manually add the title etc. Make sure to clik graph user to set.
cg_plot.Title.String = ['Pathway flux variability (n_{path}\Sigma(V_{i,max}-V_{i,min}))'];
cg_plot.Title.FontSize = 12;
axes(cg_plot);
% Create a separate subplot for the colorbar
colorbar('Location', 'southoutside');
filename = strcat(subfolder, '\', 'all');
% print(gcf,strcat(filename,'.jpeg'), '-djpeg', '-r300')  
saveas(gcf,strcat(filename,'.fig'))  
% exportgraphics(gcf,strcat(filename,'.pdf'))  
% 
% int_idx = [72,104,103,132,101,...
%             118,3,63,31,106,...
%                 62,23,15,17,...
%                 16,18,54,36,...
%                 92,70,9,1,102,...
%                 44,116,117,95,11,...
%                 137,130,45,124];
% int_idx = flip(int_idx);
% interest_label = {'Glycolysis','Pentose Phsophate Pathway','Glucuronate','TCA','Oxidative Phosphorylation',...
%                     'Pyrurvate','Alanine-Aspartate-Glutamate','Folate','Butanoate','Phenylalanine',...
%                     'Fatty acid Oxidation','FAO (n-9,mito)','FAO (even-chain,mito)','FAO (odd-chain, mito)',...
%                     'FAO (even-chain,peroxi)','FAO (odd-chain, peroxi)','Fatty acid Synthesis','Carnite shuttle (mito)',...
%                         'Lysine','Glycerophosphlipid','Arginine-Proline','Acyl-CoA hydrolysis','Co-A biosynthesis',...
%                         'CoA synthesis','Purine','Pyrimidine','N-glycan','Ascorbate',...
%                         'BCAA','Transport','methionine','Starch and sucrose'};
% interest_label = flip(interest_label);
% h2 = clustergram(matrix_row_subSys_col_tissue(int_idx,:),...
%         'RowLabels',interest_label,...
%         'ColumnLabels',modelID,...
%         'Colormap',jet,...
%         'Standardize','Row',...
%         'Showdendrogram','off',...
%         'Cluster','Row',...
%         'ColumnPdist','correlation',...
%         'ImputeFun', @knnimpute)
% 
% cg_plot = plot(h2)
% cg_plot.FontSize = 12;
% cg_plot.Title.String = ['Pathway flux variability (n_{path}\Sigma(V_{i,max}-V_{i,min}))'];
% 
% cg_plot.Title.FontSize = 12;
% axes(cg_plot);
% colorbar('Location', 'southoutside');
% 
% filename = strcat(subfolder, '\', 'all_sel');
% print(gcf,strcat(filename,'.jpeg'), '-djpeg', '-r300')  
% saveas(gcf,strcat(filename,'.fig'))  
% % exportgraphics(gcf,strcat(filename,'.fig'))  
% 
% %% Save subsystem with deviation larger than certain percentage
% % select subsystems to include in plot
% subCoverage = (matrix_row_subSys_col_tissue-mean(matrix_row_subSys_col_tissue,2))./mean(matrix_row_subSys_col_tissue,2)*100;
% inclSub = any(abs(subCoverage)>subSysThreshold,2); % return the logical columns.
% 
% h3 = clustergram(matrix_row_subSys_col_tissue(inclSub,:),...
%         'RowLabels',RowLabels(inclSub),...
%         'ColumnLabels',modelID,...
%         'Colormap',jet,...
%         'Standardize','Row',...
%         'Showdendrogram','off',...
%         'Cluster','Row',...
%         'ColumnPdist','correlation',...
%         'ImputeFun', @knnimpute)
% %         'Standardize','row',...
% % colorbar
% % colorbar
% cg_plot = plot(h3)
% cg_plot.FontSize = 12;
% cg_plot.Title.String = ['Pathway flux variability (n_{path}\Sigma(V_{i,max}-V_{i,min}))'];
% cg_plot.Title.FontSize = 12;
% cg.DisplayRatio = 0.1;
% 
% axes(cg_plot);
% colorbar('Location', 'southoutside');
% 
% filename = [subfolder, '\', 'all_sel_th_' num2str(subSysThreshold)];
% print(gcf,strcat(filename,'.jpeg'), '-djpeg', '-r300')  
% saveas(gcf,strcat(filename,'.fig'))  
% % exportgraphics(gcf,strcat(filename,'.fig'))  

end
