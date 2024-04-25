%% Define the save folder
pathway = pwd;
save_dir = 'D';
subfolder = [pathway '\' save_dir];
if ~exist(subfolder, 'dir')
    mkdir(subfolder)
end

%% Load the input
load(strcat(pathway,'\C\','out_empty2.mat'));
load(strcat(pathway,'\C\','fruitfly2_ECupd.mat'));

%% Load the Gene-EC inro
fileIn=strcat(pathway,'\extra\','Sun_Jin_EC_numbers');

%% Run the functions
[out_empty3,gem_uu,T2] = run_genelistOfEmptyEC(fileIn,out_empty2,gem_u);
movefile('out_empty3.mat',subfolder);
movefile('fruitfly3.mat',subfolder);
writetable(T2,strcat(subfolder,'\updatedECINfo.xlsx'))

function [out_empty3,gem_uu,T2] = run_genelistOfEmptyEC(fileIn,out_empty2,gem_u)
% Goal:
%   - to link EC numbers to reactions that are 1) are not linked to EC
%   based on the KEGG-ID, but 2) have the genes, which may have the EC
%   number associated with them.
% 
%   Input:
%       fileIn = 'Sun_Jin_EC_numbers';
%       out_empty2 = out_empty2; %string that contains reactions without EC
%       and, some ECs have been updated through KeggID link.
% 
%% Load the EC data and parse it for GEM format eccodes (e.g. xxx;xxxx)
file_tmp = strcat(fileIn, '.xlsx');
df = string(table2cell(readtable(file_tmp,'ReadVariableNames', false)));
df(:,5) = strrep(string(df(:,5)),'ec:','');
% df(:,{'EC number'}) = strrep(df(:,{'EC number'}),
df_gene = df(:,1); % 1st is gene_list in a format that match with GEM
df_ECnumb = df(:,5); % 5th is ec number

% concetanate ec to (e.g., x.x.x.x;y.y.y.y) for those with multiple ec.
[~, idxa, ~] = unique(df_gene,'stable'); %stable keeps the order.
df_unique = df(idxa,:);
for i = 1:length(df_unique)
    [Lia,~] = ismember(df_gene,df_unique(i)); % How many original genes are member of ith unique gene
    if sum(Lia) ~= 1% meaning that gene(Lia) have multiple
        ec_tmp = string(df_ECnumb(Lia,1)); % must be a vector [Lia x 1]
        if contains(ec_tmp,'.') == 1 % some has multiple genes but no EC
            df_unique(i,5) = strjoin(ec_tmp,';');
        end
    end
end

% Since only genes that are associated with EC will be linked to
% reacctions, first neglect the genes that are not associated with EC
not_empty_logic = ~cellfun('isempty',cellstr(df_unique(:,5)));
df_unique_u = df_unique(not_empty_logic,:);
disp([newline '# of genes with eccodes: ' num2str(length(df_unique_u)), ' of ' num2str(length(df_unique)) ' genes'])

%% Find reactions that contain the genes in the gene list above.
% retrieve idx of the reactions that have genes but no EC
tmp = cellstr(out_empty2(:,4));
tmp_logic = [~cellfun('isempty',tmp), cellfun('isempty',cellstr(out_empty2(:,5)))]; % we want 1 (not empty) and 1 (empty)
tmp_idx = find(all(tmp_logic,2));
df_orig_geneWOec = out_empty2(tmp_idx,:);

% loop the query genes to match with gene list
df_orig_geneWOec_u = df_orig_geneWOec;
numb_of_rxns_updated = [];
ec_col = strings(length(df_orig_geneWOec_u),1);
for i = 1:length(df_unique_u)
    gene = df_orig_geneWOec(:,4);
    Lia = contains(gene,df_unique_u(i,1));
    if sum(Lia)~=0
        idx = find(Lia);
%         disp(df_unique_u(i,1))
        for j = 1:length(idx)
            if contains(ec_col(idx(j),1),'.') == 1 % if the array contains other EC, then strjoin
                ec_col(idx(j),1) = strcat(ec_col(idx(j),1),';',df_unique_u(i,5));
            else   
                ec_col(idx(j),1) = df_unique_u(i,5);
            end
        end
        numb_of_rxns_updated = [numb_of_rxns_updated; idx];
    end
end

numb_of_rxns_updated_unique = unique(numb_of_rxns_updated,'stable');
disp([newline '# of reactions updated with EC by presence of genes: ' num2str(length(numb_of_rxns_updated_unique)) ,newline,...
        '# of total reactions that were not EC linked but genes: ' num2str(length(ec_col)), newline,...
        '% of reactions for updated EC: ' num2str(length(numb_of_rxns_updated_unique)) ' / ' num2str(length(ec_col)) ' = ' ,...
        num2str(length(numb_of_rxns_updated_unique)/length(ec_col))])
df_orig_geneWOec_u(:,5) = ec_col;

out_empty3 = out_empty2;
for i = 1:length(df_orig_geneWOec_u)
    [Lia,~] = ismember(out_empty2(:,1),df_orig_geneWOec_u(i));
    out_empty3(Lia,5) = df_orig_geneWOec_u(i,5);
end

%% update the EC to the model
non_empty_logic = ~cellfun('isempty',cellstr(out_empty3(:,5))); % 5th col is EC info col.
non_empty_idx = str2double(out_empty3(non_empty_logic,1));% 1st is idx of the original gem associated with reactions.
ec_to_add = cellstr(out_empty3(non_empty_logic,5));

gem_uu = gem_u;
rxns_info = [];
genes_info = [];
for i = 1:length(ec_to_add)
    gem_uu.eccodes(non_empty_idx(i),1) = ec_to_add(i);
    rxns_info{i} = gem_uu.rxns(non_empty_idx(i),1);
    genes_info{i} = gem_uu.grRules(non_empty_idx(i),1);
end
disp([newline 'total # of reactions updated (kegg+gene): ' num2str(length(ec_to_add))])
    
gem_uu.id = strcat(gem_u.id, '_ECupd_fromGene');
% save([strcat(gem_uu.id,'.mat')],'gem_uu');
genes_info = vertcat(genes_info{:});

T2 = table(rxns_info',genes_info,ec_to_add);

T = table(out_empty3);
filename = 'genes_assoc_EC_basedOn_genes.xlsx';
% write(T,filename);
save('out_empty3.mat','out_empty3')
save('fruitfly3.mat','gem_uu')
end