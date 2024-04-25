% Goal:
% % to make a data struct that is
% compatible for running the tINIT2.
% %  Input:
%       transcripts_file is of excel file, which
%       will be made into a data_struct file.
pathway = pwd;
save_dir = 'tmp';
saveFolder = [pathway '\' save_dir];
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder)
end

filename = 'FCA_pseudobulk_table_31_clusters_u.xlsx';
out = run_makeTranscrptIn_v2(filename);
movefile('transcripts_dataStruct.mat',saveFolder)

function out = run_makeTranscrptIn_v2(filename)


current_folder = pwd;
tmp2 = strcat(current_folder, '\files\', filename);
gtex_data = readtable(tmp2);
data_struct.tissues = gtex_data.Properties.VariableNames(2:end)';  % sample (tissue) names
data_struct.tissues = convertStringsToChars(string(data_struct.tissues));
data_struct.genes = gtex_data.gene;  % gene names
data_struct.levels = table2array(gtex_data(:, 2:end)); % This is the gene expression level data
data_struct.threshold = nanmean(data_struct.levels,1)'; % Used for the average value of all the gene expressions.
save('transcripts_dataStruct.mat', 'data_struct')
out = data_struct;
end