%% Goal
%   - To make heatmaps 
%   input:
%       - From the 'script_run_outerjoin_all.m'
%       - data_2
%           - data that contains
%               row: rxns,
%               col: models, 
%               elements: flux values.
%           - Final results that is 1 x m subsystems. It contains
%           comparisons of fluxes
%       - subsysout
%           - the subsystem names. Must match the size with the data_2
%       -Save_dir
%           - the parent folder is results in flux analysis
%           -e.g.,) save_dir = heatmap;
%% Load the model
% load('model_out_cbra2.mat') % cbra2 is where the model id is updated
save_dirFlux = '5_heatmap';
write = 1;
load(strcat(pwd,'\4_outerjoin\','data_store3.mat')) 
subsysOut = {'Glycolysis','Pentose phsophate pathway',...
        'Tricarboxylic acid cycle (TCA)', 'Nucleotide metabolism',...
        'Pyruvate metabolism','Alanine, aspartate and glutamate metabolism',...
        'Folate metabolism','Amino sugar metabolism','transport'...
        'Fatty acid oxidation','bouf7','bouf9','boef','boof','bobcf','boduf',...
        'boefp','boofp','bopp','boufp7','boufp9','Fatty acid biosynthesis',...
        'Fatty acid biosynthesis (even)','Cholesterol metabolism',...
        'Glycerolipid metabolism','Sphingolipid metabolism',...
        'Arginine and proline metabolism','CoA synthesis','Pantothenate and CoA biosynthesis',...
        'NADPH (cyto)','NADPH (mito)'};

run_heatmap_subsys(data_store3,subsysOut,save_dirFlux)

function run_heatmap_subsys(data_store3,subsysOut,save_dirFlux)

%% Set the save folder
pathway = pwd;
subfolder = [pathway '\' save_dirFlux];
if ~exist(subfolder, 'dir')
    mkdir(subfolder)
end

for i = 1:length(subsysOut)
    data_store3{1,i};
    rxnIDs = table2array(data_store3{1,i}(:,1));
    modelNames = data_store3{1,i}.Properties.VariableNames(2:end)';
    valsCell = table2array(data_store3{1,i}(:,2:end));
    
    % valCell normalize - normalize to the maximum range e.g., 1000. to
    % make all ranging from -1 to 1.
    valsCell = valsCell./1000;

filename = strcat(subfolder, '\', string(subsysOut{i}), '_all.jpeg');
figure(1)
h2 = heatmap(valsCell,'GridVisible','on','ColorbarVisible','on');
h2.XLabel = 'models';
h2.YLabel = ['Reactions'];
h2.XData = modelNames;
h2.YData = rxnIDs;
h2.ColorLimits = [-1 1];
h2.Colormap = jet;
h2.Title = subsysOut{i};
if i == 9 || i == 10 %transport reaction or FAO
    grid off % for transport reaction, it gets dark color
else
    grid on
end

%% Now, will take only those with some differences
diff = (valsCell-mean(valsCell,2))./mean(valsCell,2)*100;
inclSub = any(abs(diff)>0,2); % return the logical columns.
val_sel = valsCell(inclSub,:);

filename = strcat(subfolder, '\', string(subsysOut{i}), '_sel.jpeg');
figure(4)
h3 = heatmap(val_sel,'XLabel','models','YLabel',['Normalized fluxes' newline '(mmol g-DCW^{-1} hr^{-1})'],...
    'XData',modelNames,'YData',rxnIDs(inclSub),'ColorLimits',[-1 1],...
    'Colormap',jet,'GridVisible','on','ColorbarVisible','on');
h3.YLabel = ['Reactions'];
h3.Title = subsysOut{i};
if i == 9 || i == 10 %transport reaction or FAO
    grid off % for transport reaction, it gets dark color
else
    grid on
end

%% Figure 5 - Graph without the rows that contain only NaN and 0
rowsToContain = ~all(isnan(valsCell) | valsCell == 0 ,2); % get the rows that does not have only NaN and 0.
val_sel_update = valsCell(rowsToContain,:);

%% 
% Enforce NaN to 0. and Z-score. Then, remove all 0 reactions.
valsCell_u = valsCell;
nan_idx = isnan(valsCell_u);
valsCell_u(nan_idx) = 0;
valsCell_u_uu = zscore(valsCell_u,1,2);
rowsToContain = ~all((valsCell_u_uu) == 0 ,2); % get the rows that does not have only NaN and 0.
val_sel_update = valsCell_u_uu(rowsToContain,:);

filename = strcat(subfolder, '\', string(subsysOut{i}), '_sel_noNaN_0.jpeg');
figure(5)
h3 = heatmap(val_sel_update,'YLabel',['Normalized fluxes' newline '(mmol g-DCW^{-1} hr^{-1})'],...
    'XData',modelNames,'YData',rxnIDs(rowsToContain),'ColorLimits',[-3 3],...
    'Colormap',jet,'GridVisible','on','ColorbarVisible','on');
ax = gca;
ax.YDisplayLabels = nan(size(ax.YDisplayData));
h3.YLabel = ['Reaction rates'];
h3.Title = subsysOut{i};
if i == 9 || i == 10 %transport reaction or FAO
    grid off % for transport reaction, it gets dark color
else
    grid on
end
print(gcf,filename, '-djpeg', '-r300')
% saveas(gcf,strcat(subfolder, '\', string(subsysOut{i}), '_sel_noNaN_0.fig'));
% exportgraphics(gcf,strcat(subfolder, '\', string(subsysOut{i}), '_sel_noNaN_0.pdf'));

%% Figures normalized to control
%   This is done to make it easier to draw arrows - 3x, 1x, etc.
valsCell_norm_0 =(valsCell-valsCell(:,1))./valsCell(:,1)*100; % ctr is 0
valsCell_norm_1 =valsCell./valsCell(:,1); % ctr i 1

valsCell_norm_0_sel =(val_sel-val_sel(:,1))./val_sel(:,1)*100; % ctr is 0
valsCell_norm_0_sel(:,1) =abs(valsCell_norm_0_sel(:,1)); % ctr is 0
valsCell_norm_1_sel =val_sel./val_sel(:,1); % ctr i 1

end
end