%% Step 1: Get the tissue information
% goal:
%   - to have an excel file that has
%       1. tissue info
%       2. tissue info for graph
pathway = pwd;
load(strcat(pathway,'\models\','model_all.mat'))
[tissue_info_ref, tissue_info_graph,all_subsys] = run_getTissueInfo(model_all);
save('tissue_info.mat','tissue_info_ref')
save('all_subsys.mat','all_subsys')

function [tissue_info_ref, tissue_info_graph,all_subsys] = run_getTissueInfo(model_all)
for i = 1:length(model_all)
    tissue_info_ref{i,1} = model_all{1,i}.id;
end
% This will be changing depending on the tissue information annotation.
% tissue_info_graph = extractAfter(tissue_info_graph,'_');
tissue_info_ref = strrep(tissue_info_ref,'_',' ');
tissue_info_ref = strrep(tissue_info_ref,'adult','');
tissue_info_graph = tissue_info_ref;

%get the all subsystems
model_all{1, 2}.subSystems{3474, 1}=model_all{1, 2}.subSystems{3474, 1}{1, 1};
all_subsys = unique(string(model_all{1, 2}.subSystems));

tissue_info_ref_name = 'tissu_info.xlsx';
writecell(tissue_info_ref,tissue_info_ref_name);

end