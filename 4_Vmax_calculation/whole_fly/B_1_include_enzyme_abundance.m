% Goal:
%   To include the enzyme abundance to the model
%   Simply, match the gene and corresponding protein abundance.
%% 0. Set the path
path = pwd;
%% Make the output folder
output = 'B_1_gene2protein_IntoModel';
output_folder = strcat(path,'\',output);
if ~exist(output_folder,'dir')
    mkdir(output_folder)
end

currentFolder = pwd();  % Get the current folder
addpath(genpath(currentFolder));

%% 1. Load the model
load(strcat(currentFolder,'\models\','fluitfly3.mat'));
%% 2. Specify input
model_in = gem_uu;
fileName_ref = 'all';
model_out = gene2protein_intoModel(model_in,fileName_ref,path);
% save
save(strcat(output_folder,'\','model_out.mat'),'model_out');

function out = gene2protein_intoModel(model_in,fileName_ref,path)
% to append the protein number extracted from the previous A function and
% make include it to the model.
%% Extract the data of protein abundance to the corresponding genes 
filename = strcat(fileName_ref,'.xlsx');
sheetName = 'mmol_gDCW';
df_all = readtable(filename,'Sheet',sheetName);
df_all_var = df_all.Properties.VariableNames;

a = df_all{:,1};
a_val = df_all{:,2:end};

rxns = model_in.rxns;
b = model_in.genes; % Ref is where going to be updated.

%% Determine whether a is in b -- a is the reference, b is the gene in model
df_u = []; % The first column is b, 2nd column is a.
for i = 1:length(b) % loop the model genes
[Lia,~] = ismember(a(:,1),b(i,1)); % Find whethter a is part of element in b. 
%     count = count +sum(Lia);
%     disp(count)
    if sum(Lia) > 1 % Then, multiple match
        idx = find(Lia);
        tmp = a(idx,1); 
        tmp_df_cell = repmat(b(i,:),length(idx),1);
        df_tmp = [tmp_df_cell, tmp];
        df_u = [df_u; df_tmp];
    elseif sum(Lia) == 1  % That means, there are 1 to 1 match
        idx = find(Lia); % index that matches to a.
        df_tmp = [b(i,:), a(idx,1),a_val(idx,1)];
        df_u = [df_u; df_tmp];
        model_in.protein(i,1)= a_val(idx,1); %protein abundance is updated.
    else
        model_in.protein(i,1) = 0;
        continue % for 0, do not update but continue.
    end
end

%% Results
total_elements = length(model_in.protein);
empty_indices = sum(model_in.protein==0);
numb_of_filled_protein = total_elements-empty_indices;
percentage_of_rxn_with_EC = numb_of_filled_protein./total_elements*100;
disp([newline '# of filled_protein: ' num2str(numb_of_filled_protein) ' out of ' num2str(total_elements) ' ('...
        num2str(percentage_of_rxn_with_EC) ' %)' newline ...
        '# of empty protein values: ' num2str(sum(empty_indices))]);

model_in.update_B = 'protein abundance is added to the corresponding gene(mmol/gDCW)';
out = model_in;
%% Regression model to fill the missing values
% Since only 710 proteins were included in the 1846 gene sets (38.5%), I
% decided use the regression pipeline to include the rest. To do so, use
% the following. These values were determined from the protein abundance
% pipeline
slope = 77.5689;
y_int = 1.0269*10^3;
% Load the mRNA data
filename = strcat(fileName_ref,'.xlsx');
sheetName = 'original_mrna';
df_all_original_mrna = readtable(filename,'Sheet',sheetName);
df_all_var_orig = df_all_original_mrna.Properties.VariableNames;
a = df_all_original_mrna{:,1};
a_val = df_all_original_mrna{:,2:end};

% 1. Find the empty indices
idx_empty = find(model_in.protein == 0);
b = model_in.genes(idx_empty);
%% Determine whether a is in b -- a is the reference, b is the gene in model
df_u = []; % The first column is b, 2nd column is a.
for i = 1:length(b) % loop the model genes
    idx_at_original = idx_empty(i);
    [Lia,~] = ismember(a(:,1),b(i,1)); % Find whethter a is part of element in b. 
%     count = count +sum(Lia);
%     disp(count)
    if sum(Lia) > 1 % Then, multiple match
        idx = find(Lia);
        tmp = a(idx,1); 
        tmp_df_cell = repmat(b(i,:),length(idx),1);
        df_tmp = [tmp_df_cell, tmp];
        df_u = [df_u; df_tmp];
    elseif sum(Lia) == 1  % That means, there are 1 to 1 match
        idx = find(Lia); % index that matches to a.
        df_tmp = [b(i,:), a(idx,1),a_val(idx,1)];
        df_u = [df_u; df_tmp];
        
        % Calculate the ppm values
        protein_ppm = (slope.*(a_val(idx,1))+y_int);
        % convert to mmol/gDCW referring to the A scritp
        % This is from the average weight of the 2872 proteins from the
        % parameter dataset.
        maa = 67.12*1000 ; %The average molecular mass of an amino acid in Daltons (or g/mol)
        protein_mmol_gDCW = protein_ppm./1000./maa;
        model_in.protein(idx_at_original,1)= protein_mmol_gDCW; %protein abundance is updated. 
    else
        model_in.protein(idx_at_original,1) = 0;
        continue % for 0, do not update but continue.
    end
end
%% Results
total_elements = length(model_in.protein);
empty_indices = sum(model_in.protein==0);
numb_of_filled_protein = total_elements-empty_indices;
percentage_of_rxn_with_EC = numb_of_filled_protein./total_elements*100;
disp([newline '# of filled_protein: ' num2str(numb_of_filled_protein) ' out of ' num2str(total_elements) ' ('...
        num2str(percentage_of_rxn_with_EC) ' %)' newline ...
        '# of empty protein values: ' num2str(sum(empty_indices))]);

model_in.update_B = 'protein abundance is added to the corresponding gene(mmol/gDCW)';
out = model_in;
%% Fill the 0 values with the regression
% Now 1788 values were filled with the regression. 58 values are still missing.
% Find the 0 values and use the same number.
%rsq = 0.0627, cor = 0.3;l
slope = 77.5689;
y_int = 1.0269*10^3;
% 1. Find the empty indices
idx_empty = find(model_in.protein == 0);
b = model_in.genes(idx_empty);
%% Determine whether a is in b -- a is the reference, b is the gene in model
for i = 1:length(b) % loop the model genes
    idx_at_original = idx_empty(i);
    protein_ppm = (slope.*0+y_int);
    maa = 67.12*1000 ; %The average molecular mass of an amino acid in Daltons (or g/mol)
    protein_mmol_gDCW = protein_ppm./1000./maa;
    model_in.protein(idx_at_original,1)= protein_mmol_gDCW; %protein abundance is updated. 
end
%% Results
total_elements = length(model_in.protein);
empty_indices = sum(model_in.protein==0);
numb_of_filled_protein = total_elements-empty_indices;
percentage_of_rxn_with_EC = numb_of_filled_protein./total_elements*100;
disp([newline '# of filled_protein: ' num2str(numb_of_filled_protein) ' out of ' num2str(total_elements) ' ('...
        num2str(percentage_of_rxn_with_EC) ' %)' newline ...
        '# of empty protein values: ' num2str(sum(empty_indices))]);

model_in.update_B = 'protein abundance is added to the corresponding gene(mmol/gDCW)';
out = model_in;
end