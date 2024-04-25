% For A_fva_bound_results_v2
% I did the analysis in O2 because it took quite a time if I do it in my laptop.

%% Defind the savefolder
dir = 'B_fva_update_v2';
pathway = pwd;
subfolder = [pathway '\' dir];
if ~exist(subfolder, 'dir')
    mkdir(subfolder)
end

%% Load the files
dir = 'A_fva_bounds_results_v2';
loadfolder = [pathway '\' dir];

load(strcat(loadfoler,'out_all.mat'));
load(strcat(loadfoler,'subSysModel.mat'));

% Perform flux variability analysis
out_all_fvaBounded = {};
for i = 1:length(out_all)
    all_fvaBounded = out_all{i, 1}.subSysModel ;
    all_fvaBounded = changeRxnBounds(all_fvaBounded,all_fvaBounded.rxns,out_all{i,1}.models(:,1),'l');
    all_fvaBounded = changeRxnBounds(all_fvaBounded,all_fvaBounded.rxns,out_all{i,1}.models(:,2),'u');
    out_all_fvaBounded{i,1} = all_fvaBounded;
end
    save(strcat(subfolder,'\','out_all_fvaBounded.mat'),'out_all_fvaBounded');

