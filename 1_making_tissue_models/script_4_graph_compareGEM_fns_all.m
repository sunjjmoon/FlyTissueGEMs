% Goal:
%   - To compare the models functionalities.
%  Compare only the specific models

save_dir = 'fn_analz';

% tissue_info = 'tissue_info';

%% Load the data
info = strcat(pwd,'\files\','tissue_info','.xlsx');
tissue_info = readtable(info,'ReadVariableNames',true);
tissue_info = table2cell(tissue_info(:,1));

save_in = 1;

out = run_compareGEM_fns_all(model_all,tissue_info,save_in,save_dir)

function out = run_compareGEM_fns_all(model_all,tissue_info,save_in,save_dir)
%% Define the save folder
pathway = pwd;
subfolder = [pathway '\' save_dir];
if ~exist(subfolder, 'dir')
    mkdir(subfolder)
end

close all

%set the location of the tasks files
full_path = strcat(pathway,'\files\','metabolicTasks_Full_e.xlsx');
res_func = compareMultipleModels(model_all, false, false, [], true, full_path);

useModels = tissue_info;
res.modelIDs = useModels;

%% Identify which tasks differed among the GEMs (i.e., not all passed or all failed).
% obtain the rows that are different, excluding not all pass or failure.
isDiff = ~all(res_func.funcComp.matrix == 0, 2) & ~all(res_func.funcComp.matrix == 1, 2); % 2 is to return column (by after working for raws)
diffTasks = res_func.funcComp.tasks(isDiff);

%% Count the numbers
count = sum(res_func.funcComp.matrix(isDiff,:),1)';
[y2,idx] = sort(count,1,'descend');
tissue_sorted = useModels(idx);
y2_norm_zscore = normalize(y2);
T = table(tissue_sorted,y2,y2_norm_zscore);
write(T,strcat(subfolder,'\','availabel_fns_all.xlsx'));
save(strcat(subfolder,'\','T.mat'),'T');

%% visualize the matrix
out = plot(1);
spy(res_func.funcComp.matrix(isDiff,:), 30);
% apply some formatting changes
set(gca, 'XTick', 1:numel(useModels), 'XTickLabel', useModels, 'XTickLabelRotation', 90, ...
    'YTick', 1:numel(diffTasks), 'YTickLabel', diffTasks, 'YAxisLocation', 'right');
xlabel(gca, '');
legend('pass','Location','southeast')
legend('boxoff')
% title('Functional comparison')

if save_in ==1
%     print(gcf,strcat(subfolder,'\tmp.tiff'), '-dtiff', '-r600')
    saveas(gcf,strcat(subfolder,'\tmp.fig'))
end

%% Heatmap graph
 figure(5)
    h = heatmap(y2_norm_zscore','Colormap', redbluecmap,...
    'ColorbarVisible','off','GridVisible','on','CellLabelColor','none');
    h.XData = tissue_sorted; 
    h.YData = '.';
    h.YLabel = 'z-score';
%     h.YTickLabelRotation = 90;
    h.Title = 'Metabolic functional comparison';
    colorbar
    
 if save_in ==1
%     print(gcf,strcat(subfolder,'\count.tiff'), '-dtiff', '-r600')
    saveas(gcf,strcat(subfolder,'\count.fig'))

 end

end