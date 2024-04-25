function [gene2paxDB_T_raw, gene2paxDBsel] = linkGeneToPaxdb(flybase_fields_out_c_u_t,id_validation_gene_out,write)
% goal:
%   - To link the genes to the original PaxDB file.
%% Var Names
f_varNames = flybase_fields_out_c_u_t.Properties.VariableNames;
i_varNames = id_validation_gene_out.Properties.VariableNames;

%% Load the refPaxDB
flybase_fields_out_c_u_t_cell = table2cell(flybase_fields_out_c_u_t);
id_validation_gene_out_cell = table2cell(id_validation_gene_out);
%% Loop the Paxdb proteins to find the match in flybase_fields and concetanate
data_tmp = strings(length(flybase_fields_out_c_u_t_cell),3);
for i = 1:length(flybase_fields_out_c_u_t_cell)
    [Lia,~] = ismember(id_validation_gene_out_cell(:,1),flybase_fields_out_c_u_t_cell(i,2));
    if sum(Lia) ~= 0
        data_tmp(i,:) = string(id_validation_gene_out_cell(Lia,:)); % abundance col
    end
end
%% Concetanate
gene2paxDB = [flybase_fields_out_c_u_t_cell, data_tmp];

%% Remove empty string
tmp = find(gene2paxDB(:,2) == ""); % find "" string
gene2paxDB(tmp,:) = []; % remove the row;

gene2paxDB_T_raw = table(gene2paxDB);
gene2paxDB_T_raw = splitvars(gene2paxDB_T_raw);
gene2paxDB_T_raw.Properties.VariableNames = [f_varNames, i_varNames];

gene2paxDBsel = [gene2paxDB_T_raw(:,{'current_symbol'}), gene2paxDB_T_raw(:,{'abundance'})];

if write ==1
    filename = ['results_gene2PaxDBabundance.xlsx'];
    writetable(gene2paxDBsel,filename);    
end
end