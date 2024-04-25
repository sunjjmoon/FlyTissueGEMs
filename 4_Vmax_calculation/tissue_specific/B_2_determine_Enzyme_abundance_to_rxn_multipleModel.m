% Goal:
%   To include the enzyme abundance to the individual models
%       using the rxnGeneMat matrix. 
%% 0. Set the path
path = pwd;

%% Make the output folder
output = 'B_2_rxn2protein_IntoModel';
output_folder = strcat(path,'\',output);
if ~exist(output_folder,'dir')
    mkdir(output_folder)
end

%% Load the input file
fileName_ref = 'model_out';
fileName_ref_in = strcat(path,'\B_1_gene2protein_IntoModel\',fileName_ref);

model_out = link_enzyme2rxn(fileName_ref_in,path);
% save
save(strcat(output_folder,'\','model_out.mat'),'model_out');

function out = link_enzyme2rxn(fileName_ref_in,path)
% Use the rxnGeneMat to make a vector that contains the apparent enzyme
% abundance of the reaction. These basically indiate the grRules. It does
% not differentiate the on or off.
%% Extract the model and rxnGeneMat
% Reference is for the original model. 
df_ref_tmp = load(strcat(fileName_ref_in,'.mat'));
df_ref_tmp = df_ref_tmp.model_out;

%% Multiple models (2023.10.29)
for p = 1:length(df_ref_tmp)
df_ref = df_ref_tmp{1,p};

df_rxnGeneMat = df_ref.rxnGeneMat;
b = df_rxnGeneMat;

%% Extract the protein abundance info - a
a = df_ref.genes;
a_val = df_ref.protein;

%% Initial analysis
total_rxns = length(df_ref.rxns);
total_genes = length(df_ref.genes);
total_protein_without_0 = sum(a_val~=0);
disp([newline 'Initial analysis']);
disp(['# of total rxns: ' num2str(total_rxns) newline ...
         '# of total genes: ' num2str(total_genes) newline ...
          '# of total proteins with value: ' num2str(total_protein_without_0)...
            ' (' num2str(total_protein_without_0./total_genes*100) ' %, out of ' num2str(total_genes) ')']);
        
rxn_without_genes = sum(cellfun('isempty',df_ref.grRules));
rxn_with_genes = sum(~cellfun('isempty',df_ref.grRules));
disp(['# of rxns with gene linked: ' num2str(rxn_with_genes)...
        ' (' num2str(rxn_with_genes/total_rxns*100) ' %, out of ' num2str(total_rxns) ')' newline ...
       '# of rxns without gene linked: ' num2str(rxn_without_genes)...
        ' (' num2str(rxn_without_genes/total_rxns*100) ' %, out of ' num2str(total_rxns) ')'])

%% Output
% apparent enzyme abundance for rxn. If 1 gene associated, use the protein
% abundance. If multiple, use the maximum. If none, use 0. 
df_ref.AppProtein4rxns = []; % apparent enzyme abundance - If 

%% Determine the genes involved in this reaction
for i = 1:size(df_rxnGeneMat,1) % loop the rxns
    idx_gene_associated_to_rxn = find(b(i,:)==1);

    % add the protein abundance to the reactions
    df_tmp = [];
    for j = 1:length(idx_gene_associated_to_rxn)
        protein_val = df_ref.protein(idx_gene_associated_to_rxn,1);
        df_tmp = [df_tmp;protein_val];
    end
    df = max(df_tmp); % if there are more than 1 protein associated, use the maximum values.
    if isempty(df) ~=1 % is gene is associated to reaction
        df_ref.AppProtein4rxns(i,1) = df;
    else % if gene is not associated to reaction, just have 0 value.
        df_ref.AppProtein4rxns(i,1) = 0;
    end
end
df_ref.update_B2 = 'apparent protein abundance field is made';
out{1,p} = df_ref;

%% Results:
total_rxns = length(df_ref.rxns);
total_genes = length(df_ref.genes);
total_App_protein_without_0 = sum(df_ref.AppProtein4rxns~=0);
disp([newline 'Results']);
disp(df_ref.id)
disp(['# of total rxns: ' num2str(total_rxns) newline ...
          '# of apparent protein Abundance added: ' num2str(total_App_protein_without_0)...
            ' (' num2str(total_App_protein_without_0./total_rxns*100) ' %, out of ' num2str(total_rxns) ')']);

disp(['# of total rxns with gene linked: ' num2str(rxn_with_genes) newline ...
        '# of apparent protein Abundance added: ' num2str(total_App_protein_without_0)...
            ' (' num2str(total_App_protein_without_0./rxn_with_genes*100) ' %, out of ' num2str(rxn_with_genes) ')']);


%% Export
df_ref.update_B2 = 'protein abundance is added to the corresponding gene(mmol/gDCW)';

out{1,p} = df_ref;
end

end