% Goal:
%   - To fill the EC with loose constraints.
%% 0. Set the path
path = pwd;
% Make the output folder
output = 'E_fill_the_EC_gap_otherOrganism';
output_folder = strcat(path,'\',output);
if ~exist(output_folder,'dir')
    mkdir(output_folder)
end

%% Load the input file
fileName_kcat = 'brenda_Caenorhabditis elegans_kcat';
fileName_kcat_in = strcat(path,'\input_data_kinetics\',fileName_kcat);
% 
fileName_kcat = 'brenda_Danio rerio_kcat';
fileName_kcat_in = strcat(path,'\input_data_kinetics\',fileName_kcat);

fileName_kcat = 'brenda_Mus_musculus_kcat';
fileName_kcat_in = strcat(path,'\input_data_kinetics\',fileName_kcat);

% 
fileName_kcat = 'brenda_homoSapiens_kcat';
fileName_kcat_in = strcat(path,'\input_data_kinetics\',fileName_kcat);

fileName_ref = 'model_out';
fileName_ref_in = strcat(path,'\D_fill_the_EC_gap\',fileName_ref);

%% Set the parameter values
write = 0;
write_all = 1; % Write all

%% Calculate the protein abundance from gene
model_out = fill_the_EC_gap_otherOrganism(fileName_ref_in,fileName_kcat_in,write,path,write_all);
% save
save(strcat(output_folder,'\','model_out.mat'),'model_out');
% save('gene2protein.mat', 'gene2protein');

function out = fill_the_EC_gap_otherOrganism(fileName_ref_in,fileName_kcat_in,write,path,write_all)

%% Extract the reference
% Reference is for the original model. 
df_ref_tmp = load(strcat(fileName_ref_in,'.mat'));
df_ref_tmp = df_ref_tmp.model_out;

%% Multiple models (2023.10.29)
for p = 1:length(df_ref_tmp)
df_ref = df_ref_tmp{1,p};

df_ref_EC = df_ref.eccodes;
% kcat file, (to be appended, B)
filename = strcat(fileName_kcat_in,'.xlsx');
df_kcat = readtable(filename);
df_kcat_var = df_kcat.Properties.VariableNames;

%% Loop the kcat (a) and see if kcat is in the reference (or the EC numbers of the model)
a = df_kcat{:,1}; % EC number in the kcat file -- data to be added to b.
a_val = df_kcat{:,2};

% Initial conditions
empty_indices = cellfun('isempty',df_ref.eccodes);
non_empty_indices = length(df_ref.kcat)-sum(empty_indices);

disp([newline 'Original, kcat filled: ' num2str(sum(df_ref.kcat~=0)) ' out of ' num2str(length(df_ref.kcat)) ' (',...
        num2str(sum(df_ref.kcat~=0)/length(df_ref.kcat)*100), ' %), missing kcat: ',...
        num2str(sum(df_ref.kcat==0))]);    
disp(['Original, kcat filled: ' num2str(sum(df_ref.kcat~=0)) ' out of ' num2str(non_empty_indices) ' (',...
        num2str(sum(df_ref.kcat~=0)/non_empty_indices*100), ' %), missing kcat: ',...
        num2str(sum(df_ref.kcat==0)-sum(empty_indices))]);

%% 1. extract EC numbers of the first three values. - If duplicates exist, max values are used. 1.1.1.x
% Use cellfun and regexp to extract the first three parts from each element in a
EC_first_three = cellfun(@(str) regexp(str, '\d+\.\d+\.\d+', 'match', 'once'), a, 'UniformOutput', false);
% Find unique EC number and the corresponding values
% Iterate through unique values and find the maximum corresponding value in c
EC_first_three_unique = unique(EC_first_three);
a_val_unique = zeros(length(EC_first_three_unique),1);
for i = 1:length(EC_first_three_unique)
    EC_first_three_tmp = EC_first_three_unique{i};
    [Lia,~] = ismember(EC_first_three, EC_first_three_tmp);
    if sum(Lia) > 1 % if there are duplicates, get the only the maxiimum value
        idx_tmp = find(string(EC_first_three) == EC_first_three_tmp);
        a_val_unique(i,1) = max(a_val(idx_tmp));
    else
        a_val_unique(i,1) = a_val(Lia);
    end
end
%% Determine whether EC numbers from kcat val (a) is in the original model EC (b)
idx_not_filled_yet = find(df_ref.kcat ==0);
b = df_ref_EC(idx_not_filled_yet,1);

df_u = []; % The first column is b, 2nd column is a.
idx_not_filled = [];
for i = 1:length(b) % loop the original EC 
    % Use cellfun and regexp to extract the first three parts from each element in a
    b_tmp = strsplit(b{i,1},';')'; % to check whether the rxn has multiple EC.
    if ~isempty(b_tmp{1})
        if length(b_tmp) > 1 % If there are more than 1 EC numbers to the corresponding rxns, extract the max kcat, and add only the highest num to that reaction
            b_tmp_tmp = [];
            b_tmp = cellfun(@(str) regexp(str, '\d+\.\d+\.\d+', 'match', 'once'), b_tmp, 'UniformOutput', false);
            for k = 1:length(b_tmp)
                [Lia,~] = ismember(EC_first_three_unique(:,1),b_tmp(k,1)); % Find whethter a is part of element in b. 
                idx = find(Lia); % index that matches to a.
                kcat_val = a_val_unique(idx,1);
                b_tmp_tmp = [b_tmp_tmp; kcat_val];
            end
                df_tmp = max(b_tmp_tmp); %get the maximum kcat value.
                if sum(df_tmp) ~=0 % values are present
                    df_ref.kcat(idx_not_filled_yet(i),1) = df_tmp;
                else
                   df_ref.kcat(idx_not_filled_yet(i),1) = 0;
                   idx_not_filled= [idx_not_filled; idx_not_filled_yet(i)]; % because this is filled. we say 0.
                end
        else
            pattern = '^\d+\.\d+\.\d+'; %% \d+: \d is a shorthand for matching any digit (0-9), and the + means to match one or more occurrences of the preceding element. So, \d+ matches one or more digits.
                                     %   \.: This is used to match a literal period ('.') character. Since the period is a special character in regular expressions, it needs to be escaped with a backslash to match a literal period.
             b_tmp = regexp(b_tmp, pattern, 'match');

             if any(cellfun('isempty', b_tmp)) == 0 % contains no empty cell
                [Lia,~] = ismember(string(EC_first_three_unique(:,1)),string(b_tmp)); % Find whethter a is part of element in b. 
                idx = find(Lia); % index that matches to a.

                if sum(Lia) == 1  % if exist, fill.
                    df_tmp = [b_tmp, EC_first_three(idx,1),a_val_unique(idx,1)];
                    df_u = [df_u; df_tmp];
                    df_ref.kcat(idx_not_filled_yet(i),1)= a_val_unique(idx,1); %protein abundance is updated.     
                else
                    df_ref.kcat(idx_not_filled_yet(i),1) = 0;
                    idx_not_filled= [idx_not_filled; idx_not_filled_yet(i)]; % because this is filled. we say 0.
                end
             end
        end
    else
        continue
    end
end
idx_not_filled = sum(df_ref.kcat ==0); % Determine the index that is not filled.
idx_not_filled = sum(idx_not_filled ~=0); % Determine the index that is not filled.

%% results
numb_of_kcat_without_0 = sum(df_ref.kcat~=0);
disp([newline '1.1.1.x, # of kcat updte: ' num2str(numb_of_kcat_without_0) ' out of ' num2str(length(df_ref.kcat)),...
        ' (' num2str((numb_of_kcat_without_0/length(df_ref.kcat)*100)) '%) from ohter Organism, ' num2str(idx_not_filled),...
        ' not filled'])
disp(['1.1.1.x, # of kcat updte: ' num2str(numb_of_kcat_without_0) ' out of ' num2str(non_empty_indices),...
        ' (' num2str((numb_of_kcat_without_0/non_empty_indices*100)) '%) from ohter Organism, ' num2str(non_empty_indices-numb_of_kcat_without_0),...
        ' not filled'])        
%% 2. Now, move on to 2 digits -- e.g.,) 1.1.x.x
% Use cellfun and regexp to extract the first three parts from each element in a
EC_first_three = cellfun(@(str) regexp(str, '\d+\.\d+', 'match', 'once'), a, 'UniformOutput', false);
% Find unique EC number and the corresponding values
% Iterate through unique values and find the maximum corresponding value in c
EC_first_three_unique = unique(EC_first_three);
a_val_unique = zeros(length(EC_first_three_unique),1);
for i = 1:length(EC_first_three_unique)
    EC_first_three_tmp = EC_first_three_unique{i};
    [Lia,~] = ismember(EC_first_three, EC_first_three_tmp);
    if sum(Lia) > 1 % if there are duplicates, get the only the maxiimum value
        idx_tmp = find(string(EC_first_three) == EC_first_three_tmp);
        a_val_unique(i,1) = max(a_val(idx_tmp));
    else
        a_val_unique(i,1) = a_val(Lia);
    end
end
%% Determine whether EC numbers from kcat val (a) is in the original model EC (b)
idx_not_filled_yet = find(df_ref.kcat ==0);
b = df_ref_EC(idx_not_filled_yet,1);

df_u = []; % The first column is b, 2nd column is a.
idx_not_filled = [];
for i = 1:length(b) % loop the original EC 
    % Use cellfun and regexp to extract the first three parts from each element in a
    b_tmp = strsplit(b{i,1},';')'; % to check whether the rxn has multiple EC.
    if ~isempty(b_tmp{1})
        if length(b_tmp) > 1 % If there are more than 1 EC numbers to the corresponding rxns, extract the max kcat, and add only the highest num to that reaction
            b_tmp_tmp = [];
            b_tmp = cellfun(@(str) regexp(str, '\d+\.\d+', 'match', 'once'), b_tmp, 'UniformOutput', false);
            for k = 1:length(b_tmp)
                [Lia,~] = ismember(EC_first_three_unique(:,1),b_tmp(k,1)); % Find whethter a is part of element in b. 
                idx = find(Lia); % index that matches to a.
                kcat_val = a_val_unique(idx,1);
                b_tmp_tmp = [b_tmp_tmp; kcat_val];
            end
                df_tmp = max(b_tmp_tmp); %get the maximum kcat value.
                if sum(df_tmp) ~=0 % values are present
                    df_ref.kcat(idx_not_filled_yet(i),1) = df_tmp;
                else
                   df_ref.kcat(idx_not_filled_yet(i),1) = 0;
                   idx_not_filled= [idx_not_filled; idx_not_filled_yet(i)]; % because this is filled. we say 0.
                end
        else
            pattern = '^\d+\.\d+'; %% \d+: \d is a shorthand for matching any digit (0-9), and the + means to match one or more occurrences of the preceding element. So, \d+ matches one or more digits.
                                     %   \.: This is used to match a literal period ('.') character. Since the period is a special character in regular expressions, it needs to be escaped with a backslash to match a literal period.
             b_tmp = regexp(b_tmp, pattern, 'match');
             if any(cellfun('isempty', b_tmp)) == 0 % contains no empty cell
                [Lia,~] = ismember(string(EC_first_three_unique(:,1)),string(b_tmp)); % Find whethter a is part of element in b. 
                idx = find(Lia); % index that matches to a.

                if sum(Lia) == 1  % if exist, fill.
                    df_tmp = [b_tmp, EC_first_three(idx,1),a_val_unique(idx,1)];
                    df_u = [df_u; df_tmp];
                    df_ref.kcat(idx_not_filled_yet(i),1)= a_val_unique(idx,1); %protein abundance is updated.     
                else
                    df_ref.kcat(idx_not_filled_yet(i),1) = 0;
                    idx_not_filled= [idx_not_filled; idx_not_filled_yet(i)]; % because this is filled. we say 0.
                end
             end
        end
    else
        continue
    end
end
idx_not_filled = sum(df_ref.kcat ==0); % Determine the index that is not filled.
idx_not_filled = sum(idx_not_filled ~=0); % Determine the index that is not filled.

%% results
numb_of_kcat_without_0 = sum(df_ref.kcat~=0);
disp([newline '1.1.x.x, # of kcat updte: ' num2str(numb_of_kcat_without_0) ' out of ' num2str(length(df_ref.kcat)),...
        ' (' num2str((numb_of_kcat_without_0/length(df_ref.kcat)*100)) '%) from ohter Organism, ' num2str(idx_not_filled),...
        ' not filled'])
disp(['1.1.x.x, # of kcat updte: ' num2str(numb_of_kcat_without_0) ' out of ' num2str(non_empty_indices),...
        ' (' num2str((numb_of_kcat_without_0/non_empty_indices*100)) '%) from ohter Organism, ' num2str(non_empty_indices-numb_of_kcat_without_0),...
        ' not filled'])        
%% 3. Now, move on to 2 digits -- e.g.,) 1.x.x.x
% Use cellfun and regexp to extract the first three parts from each element in a
EC_first_three = cellfun(@(str) regexp(str, '\d+', 'match', 'once'), a, 'UniformOutput', false);
% Find unique EC number and the corresponding values
% Iterate through unique values and find the maximum corresponding value in c
EC_first_three_unique = unique(EC_first_three);
a_val_unique = zeros(length(EC_first_three_unique),1);
for i = 1:length(EC_first_three_unique)
    EC_first_three_tmp = EC_first_three_unique{i};
    [Lia,~] = ismember(EC_first_three, EC_first_three_tmp);
    if sum(Lia) > 1 % if there are duplicates, get the only the maxiimum value
        idx_tmp = find(string(EC_first_three) == EC_first_three_tmp);
        a_val_unique(i,1) = max(a_val(idx_tmp));
    else
        a_val_unique(i,1) = a_val(Lia);
    end
end
%% Determine whether EC numbers from kcat val (a) is in the original model EC (b)
idx_not_filled_yet = find(df_ref.kcat ==0);
b = df_ref_EC(idx_not_filled_yet,1);

df_u = []; % The first column is b, 2nd column is a.
idx_not_filled = [];
for i = 1:length(b) % loop the original EC 
    % Use cellfun and regexp to extract the first three parts from each element in a
    b_tmp = strsplit(b{i,1},';')'; % to check whether the rxn has multiple EC.
    if ~isempty(b_tmp{1})
        if length(b_tmp) > 1 % If there are more than 1 EC numbers to the corresponding rxns, extract the max kcat, and add only the highest num to that reaction
            b_tmp_tmp = [];
            b_tmp = cellfun(@(str) regexp(str, '\d+', 'match', 'once'), b_tmp, 'UniformOutput', false);
            for k = 1:length(b_tmp)
                [Lia,~] = ismember(EC_first_three_unique(:,1),b_tmp(k,1)); % Find whethter a is part of element in b. 
                idx = find(Lia); % index that matches to a.
                kcat_val = a_val_unique(idx,1);
                b_tmp_tmp = [b_tmp_tmp; kcat_val];
            end
                df_tmp = max(b_tmp_tmp); %get the maximum kcat value.
                if sum(df_tmp) ~=0 % values are present
                    df_ref.kcat(idx_not_filled_yet(i),1) = df_tmp;
                else
                   df_ref.kcat(idx_not_filled_yet(i),1) = 0;
                   idx_not_filled= [idx_not_filled; idx_not_filled_yet(i)]; % because this is filled. we say 0.
                end
        else
            pattern = '^\d+'; %% \d+: \d is a shorthand for matching any digit (0-9), and the + means to match one or more occurrences of the preceding element. So, \d+ matches one or more digits.
                                     %   \.: This is used to match a literal period ('.') character. Since the period is a special character in regular expressions, it needs to be escaped with a backslash to match a literal period.
             b_tmp = regexp(b_tmp, pattern, 'match');
             if any(cellfun('isempty', b_tmp)) == 0 % contains no empty cell
                [Lia,~] = ismember(string(EC_first_three_unique(:,1)),string(b_tmp)); % Find whethter a is part of element in b. 
                idx = find(Lia); % index that matches to a.

                if sum(Lia) == 1  % if exist, fill.
                    df_tmp = [b_tmp, EC_first_three(idx,1),a_val_unique(idx,1)];
                    df_u = [df_u; df_tmp];
                    df_ref.kcat(idx_not_filled_yet(i),1)= a_val_unique(idx,1); %protein abundance is updated.     
                else
                    df_ref.kcat(idx_not_filled_yet(i),1) = 0;
                    idx_not_filled= [idx_not_filled; idx_not_filled_yet(i)]; % because this is filled. we say 0.
                end
             end
        end
    else
        continue
    end
end
idx_not_filled = sum(df_ref.kcat ==0); % Determine the index that is not filled.
idx_not_filled = sum(idx_not_filled ~=0); % Determine the index that is not filled.

%% results
numb_of_kcat_without_0 = sum(df_ref.kcat~=0);
disp([newline '1.x,x,x, # of kcat updte: ' num2str(numb_of_kcat_without_0) ' out of ' num2str(length(df_ref.kcat)),...
        ' (' num2str((numb_of_kcat_without_0/length(df_ref.kcat)*100)) '%) from ohter Organism, ' num2str(idx_not_filled),...
        ' not filled' ])
disp(['1.x.x.x, # of kcat updte: ' num2str(numb_of_kcat_without_0) ' out of ' num2str(non_empty_indices),...
        ' (' num2str((numb_of_kcat_without_0/non_empty_indices*100)) '%) from ohter Organism, ' num2str(non_empty_indices-numb_of_kcat_without_0),...
        ' not filled'])
    
df_ref.kcat = df_ref.kcat; %convert 1/s -> 1/hr

df_ref.update_E = 'Higher organism EC considered';
out{1,p} = df_ref;
  
%% 3. Find the empty indices for the EC number
empty_indices = cellfun('isempty',df_ref.eccodes);
numb_of_non_empty_EC_rxns = length(df_ref.eccodes)-sum(empty_indices);
percentage_of_rxn_with_EC = numb_of_non_empty_EC_rxns./length(df_ref.eccodes)*100;
disp([newline 'number of EC annotated to rxns: ' num2str(numb_of_non_empty_EC_rxns) newline...
        '% of EC annotated: ' num2str(percentage_of_rxn_with_EC) ' %']);
disp(df_ref.id)
    
end
end
    