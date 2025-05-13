% Goal:
%   - To combine tissue dtaa
%   - For details, check the structuring_data.xlsx and transcriptomics data
%   strategy file.
% clear
mat = dir('*.mat');
model_all_struct = [];
for i = 1:length(mat)
   tmp = load(mat(i).name).model_out;
   model_all_struct = [model_all_struct;tmp'];
   model_all{i} = tmp; 
end

save_name = 'model_all.mat';
save(save_name,'model_all')