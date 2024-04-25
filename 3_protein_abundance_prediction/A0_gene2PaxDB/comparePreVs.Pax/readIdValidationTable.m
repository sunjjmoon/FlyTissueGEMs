function id_validation_gene_out = readIdValidationTable(file_IdVa)
% Read the tsv files obtained from the flybase batch downloads
% Goal:
% Input:
%   1. file_IdVa = 'id_validation_table_flyFieldsGenes';

fileload = strcat(file_IdVa,'.txt');
id_validation_gene_out = readtable(fileload,'delimiter', '\t', 'ReadVariableNames',true);

%% Multiple genes may be mapped. The first mapped gene is the one to be exported.
[~,ia,~] = unique(table2cell(id_validation_gene_out(:,1)),'stable');
id_validation_gene_out = id_validation_gene_out(ia,:);
end