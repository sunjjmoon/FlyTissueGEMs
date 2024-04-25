% Goal:
%   - To link the kcat value to genes. Then, I will connect the genes to
%   reactions
%% 0. Set the path
path = pwd;

% Make the output folder
output = 'C_con_kcatEC_rxn';
output_folder = strcat(path,'\',output);
if ~exist(output_folder,'dir')
    mkdir(output_folder)
end

%% Load the input file
fileName_kcat = 'brenda_drosophila_kcat';
fileName_kcat_in = strcat(path,'\input_data_kinetics\',fileName_kcat);

fileName_ref = 'model_out';
fileName_ref_in = strcat(path,'\B_2_rxn2protein_IntoModel\',fileName_ref);


%% Set the parameter values
write = 0;
write_all = 1; % Write all

%% Calculate the protein abundance from gene
model_out = con_kcatEC_rxn(fileName_ref_in,fileName_kcat_in,write,path,write_all);
% save
save(strcat(output_folder,'\','model_out.mat'),'model_out');

function out = con_kcatEC_rxn(fileName_ref_in,fileName_kcat_in,write,path,write_all)

%% Extract the reference
% Reference is for the original model. 
% df_ref = load(strcat(fileName_ref_in,'.mat'));
% df_ref = df_ref.model_out;
df_ref_tmp = load(strcat(fileName_ref_in,'.mat'));
df_ref_tmp = df_ref_tmp.model_out;

%% Multiple models (2023.10.29)
for p = 1:length(df_ref_tmp)
df_ref = df_ref_tmp{1,p};

df_ref_EC = df_ref.eccodes;
df_ref_grRules = df_ref.grRules;

% kcat file, (to be appended, B)
filename = strcat(fileName_kcat_in,'.xlsx');
df_kcat = readtable(filename);
df_kcat_var = df_kcat.Properties.VariableNames;

%% Loop the kcat (a) and see if kcat is in the reference (or the EC numbers of the model)
a = df_kcat{:,1}; % EC number in the kcat file -- data to be added to b.
a_val = df_kcat{:,2};

ref = df_ref_EC;
b = ref;

%% Determine whether a is in b -- a is the data to be added to b, b is the gene in model
df_u = []; % The first column is b, 2nd column is a.
kcat_store =[];
idx_not_filled = [];
for i = 1:length(b) % loop the rxn list
    b_tmp = strsplit(b{i,1},';')'; % to check whether the rxn has multiple EC.
    
    if length(b_tmp) > 1 % If there are more than 1 EC numbers to the corresponding rxns, extract the max kcat, and add only the highest num to that reaction
        kcat_tmp_tmp = [];
        for k = 1:length(b_tmp)
            [Lia,~] = ismember(a(:,1),b_tmp(k,1)); % Find whethter a is part of element in b. 
            idx = find(Lia); % index that matches to a.
            kcat_val = a_val(idx,1);
            kcat_tmp_tmp = [kcat_tmp_tmp; kcat_val];
        end
            df_tmp = max(kcat_tmp_tmp); %get the maximum kcat value.
            if sum(df_tmp) ~=0 % values are present
                df_ref.kcat(i,1) = df_tmp;
                idx_not_filled(i,1) = 0; % because this is filled. we say 0.
            else
               df_ref.kcat(i,1) = 0;
               idx_not_filled(i,1) = i; % because this is filled. we say 0.
            end
    else
        [Lia,~] = ismember(a(:,1),b(i,1)); % Find whethter a is part of element in b. 
        idx = find(Lia); % index that matches to a.
        
        if sum(Lia) == 1  % if exist, fill.
            df_tmp = [b(i,:), a(idx,1),a_val(idx,1)];
            df_u = [df_u; df_tmp];
            df_ref.kcat(i,1)= a_val(idx,1); %protein abundance is updated.     
            idx_not_filled(i,1) = 0; % because this is filled. we say 0.            
        else
            df_ref.kcat(i,1) = 0;
            idx_not_filled(i,1) = i; % index that is not filled.
%             continue % for 0, do not update but continue.
        end
    end
end

idx_not_filled = sum(idx_not_filled ~=0); % Determine the index that is not filled.
%% results
empty_indices = cellfun('isempty',df_ref.eccodes);
non_empty_indices = length(df_ref.kcat)-sum(empty_indices);
numb_of_kcat_without_0 = sum(df_ref.kcat~=0);

disp([newline '# of kcat updte: ' num2str(numb_of_kcat_without_0) ' out of ' num2str(length(df_ref.kcat)),...
        ' (' num2str((numb_of_kcat_without_0/length(df_ref.kcat)*100)) '%) from Drosophila, ' num2str(idx_not_filled),...
        ' not filled'])
disp(['# of kcat updte: ' num2str(numb_of_kcat_without_0) ' out of ' num2str(non_empty_indices),...
        ' (' num2str((numb_of_kcat_without_0/non_empty_indices*100)) '%, out of non-empty EC) from Drosophila, '...
        num2str(non_empty_indices-numb_of_kcat_without_0),...
        ' not filled'])
    
%% Now check with other organism and fill the EC numbers without kcat values from Drosophila analysis
fileName_kcat_name = {'brenda_Caenorhabditis elegans_kcat','brenda_Danio rerio_kcat',...
                        'brenda_Mus_musculus_kcat','brenda_homoSapiens_kcat'};
numb_of_kcat_filled_now = [numb_of_kcat_without_0];
for j = 1:length(fileName_kcat_name)
    % kcat 
    fileName_kcat_in = string(fileName_kcat_name(j));
    filename = strcat(path,'\input_data_kinetics\',fileName_kcat_in,'.xlsx');
    df_kcat = readtable(filename);
    df_kcat_var = df_kcat.Properties.VariableNames;

    %% Loop the kcat (a) and see if kcat is in the reference (or the EC numbers of the model)
    a = df_kcat{:,1}; % EC number in the kcat file -- data to be added to b.
    a_val = df_kcat{:,2};

    %% determine EC array that is not filled
    idx_not_filled_yet = find(df_ref.kcat ==0);
    b_diffOrgnaism = ref(idx_not_filled_yet,1);
    %% Determine whether a is in b -- a is the data to be added to b, b is the gene in model
    df_u = []; % The first column is b, 2nd column is a.
    idx_not_filled_diffOrganism = [];
    
for i = 1:length(idx_not_filled_yet) % loop the rxn list
    b_tmp = strsplit(b_diffOrgnaism{i,1},';')'; % to check whether the rxn has multiple EC.
    
    if length(b_tmp) > 1 % If there are more than 1 EC numbers to the corresponding rxns, extract the max kcat, and add only the highest num to that reaction
        kcat_tmp_tmp = [];
        for k = 1:length(b_tmp)
            [Lia,~] = ismember(a(:,1),b_tmp(k,1)); % Find whethter a is part of element in b. Lia is the index from a.
            idx = find(Lia); % index that matches to a.
            kcat_val = a_val(idx,1);
            kcat_tmp_tmp = [kcat_tmp_tmp; kcat_val];
        end
            df_tmp = max(kcat_tmp_tmp); %get the maximum kcat value.
            if isempty(df_tmp) ~=1 % values are present
                df_ref.kcat(idx_not_filled_yet(i),1) = df_tmp;
            else
               df_ref.kcat(idx_not_filled_yet(i),1) = 0;
               idx_not_filled_diffOrganism = [idx_not_filled_diffOrganism; idx_not_filled_yet(i)]; % 
            end
    else % only 1 EC per reaction
        [Lia,~] = ismember(a(:,1),b_diffOrgnaism(i,1)); % Find whethter a is part of element in b. 
        idx = find(Lia); % index that matches to a, which is the kcat values of diff animals
        
        if sum(Lia) == 1  % if exist, fill.
            df_tmp = [b_diffOrgnaism(i,:), a(idx,1),a_val(idx,1)];
            df_u = [df_u; df_tmp];
            df_ref.kcat(idx_not_filled_yet(i),1)= a_val(idx,1); %protein abundance is updated.
        else
            df_ref.kcat(idx_not_filled_yet(i),1) = 0;
            idx_not_filled_diffOrganism = [idx_not_filled_diffOrganism; idx_not_filled_yet(i)]; % 
        end
    end
end

idx_filled_tmp = idx_not_filled_diffOrganism(idx_not_filled_diffOrganism == 0);

%% results
numb_of_kcat_without_0 = sum(df_ref.kcat~=0);
numb_of_kcat_filled_now = [numb_of_kcat_filled_now; numb_of_kcat_without_0];

empty_indices = cellfun('isempty',df_ref.eccodes);
non_empty_indices = length(df_ref.kcat)-sum(empty_indices);

numb_of_kcat_filled_this_organism = numb_of_kcat_filled_now(end,1)-numb_of_kcat_filled_now(end-1,1);
disp([newline 'Loop ' num2str(j) ', ' char(fileName_kcat_in) ' : ' num2str(length(numb_of_kcat_filled_now)) ,...
    ' filled out of ' num2str([length(idx_not_filled_diffOrganism)+length(idx_filled_tmp)])]);
disp(['Accumulated # of kcat updte: ' num2str(numb_of_kcat_without_0) ' out of ' num2str(length(df_ref.kcat)),...
        ' (' num2str((numb_of_kcat_without_0/length(df_ref.kcat)*100)) '%), ' num2str(length(idx_not_filled_diffOrganism)),...
        ' not filled'])
disp(['Accumulated of kcat updte: ' num2str(numb_of_kcat_without_0) ' out of ' num2str(non_empty_indices),...
        ' (' num2str((numb_of_kcat_without_0/non_empty_indices*100)) '%, out of non-empty EC) from Drosophila, '...
         num2str(non_empty_indices-numb_of_kcat_without_0),...
        ' not filled '])
end
%% 3. Find the empty indices for the EC number
empty_indices = cellfun('isempty',df_ref.eccodes);
numb_of_non_empty_EC_rxns = length(df_ref.eccodes)-sum(empty_indices);
percentage_of_rxn_with_EC = numb_of_non_empty_EC_rxns./length(df_ref.eccodes)*100;
disp([newline 'number of EC annotated to rxns: ' num2str(numb_of_non_empty_EC_rxns) ' out of ' num2str(length(df_ref.eccodes)) newline ...
        'not annotated with EC: ' num2str(sum(empty_indices)) newline ...
        '% of EC annotated: ' num2str(percentage_of_rxn_with_EC) ' %']);
disp(df_ref.id)

df_ref.update_C = 'kcat value added to corresponding EC(1/s)';

out{1,p} = df_ref;
end
end


