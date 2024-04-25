%% Simplify the ID
%  Start the script 1 first and then do this.
pathway = pwd;
load(strcat(pathway,'\models\','model_out_cbra.mat'))
load('tissue_info.mat');

%% Since some of the names are too long, manually change
% tissue_info_ref{1,1} = 'CNSglialCell';
% tissue_info_ref{6,1} = 'PeriphernalNervSys';
tissue_info_ref{7,1} = 'ReticularNeuroAssoGlialCell';
tissue_info_ref{23,1} = 'legTasteBristleChemoNeuron';
tissue_info_ref{29,1} = 'prefollicleStalkCell';

[model_out_cbra2, modelID_names] = run_modelIDchange(model_out_cbra,tissue_info_ref);

save('tissue_info.mat','tissue_info_ref')


function [model_out_cbra2, modelID_names] = run_modelIDchange(model_out_cbra,tissue_info)
% goal:
%   -to change the modelID such that I may be able to easily make a folder
%   for each tissues, etc
%   - model_in
%           - The cobraformate
%   tissue_info = 'tissue_info_simplified'
%   model_in = model_all
%% Load the data
% path = pwd();
% info = [path, '\', tissue_info '.xlsx'];
% tissue_info_ref
% tissu_info = readtable(info,'ReadVariableNames',false);
% tmp = genes_GEM(contains(genes_GEM.Rank,'high'),:);
% tissue_info = table2cell(tissue_info_ref(:,1));
model_out_cbra2 = model_out_cbra;
for i = 1:length(model_out_cbra)
   model_out_cbra2{i,1}.modelID = tissue_info{i,1};
end

saveName = strcat(pwd,'\models\','model_out_cbra2.mat');
save(saveName,'model_out_cbra2')

modelID_names = tissue_info;


end
