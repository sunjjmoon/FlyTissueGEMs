%% Goal
%   - to concetanate the ppm values. 
%       a : reference (x-axis later)
%       b : to be appended. (y-axis later)
%% Define the save folder
pathway = pwd;
save_dir = 'C';
subfolder = [pathway '\' save_dir];
if ~exist(subfolder, 'dir')
    mkdir(subfolder)
end

%% Load inputs
a_in = strcat(pathway,'\A0_gene2PaxDB\comparePreVs.Pax\','results_gene2PaxDBabundance'); % v2 is for the whole fly only
b_in = strcat(pathway,'\A\','results_gene2param_T2uuu'); % v2 is for the whole fly only
write = 1;

%% Run
a2b = concetanateAtoB_protein(a_in,b_in,write);
T = tableJoin(a_in,b_in,write);

%% Move files
concat.a2b = a2b; concat.T = T;
save(strcat(subfolder,'\','concat.mat'),'concat');
movefile('outerJoin_b_v2.xlsx',subfolder);

function a2b = concetanateAtoB_protein(a_in,b_in,write)
%% Extract the data, df_param (Ref)
filename = strcat(a_in,'.xlsx');
a = readtable(filename);
a_var = a.Properties.VariableNames;

% mrna file, (to be appended, B)
filename = strcat(b_in,'.xlsx');
b = readtable(filename);
b_var = b.Properties.VariableNames;
b_var(1) = {'genes'};

%% Loop the df_param to append the associated mRNA abundances
a_ref = table2cell(a);
b_cell = table2cell(b);
%% Loop the param (Ref) to find the matched array for query (a, mRNA)
df_u = [];
for i = 1:length(a_ref)
    [Lia,~] = ismember(b_cell(:,1),a_ref(i,1));
%     count = count +sum(Lia);
%     disp(count)
    if sum(Lia) > 1 % Then, multiple match
        idx = find(Lia);
        tmp = b_cell(idx,:); 
        tmp_df_cell = repmat(a_ref(i,:),length(idx),1);
        df_tmp = [tmp_df_cell, tmp];
        df_u = [df_u; df_tmp];
    elseif sum(Lia) == 1  % That means, there are 1 to 1 match
        idx = find(Lia);
        df_tmp = [a_ref(i,:), b_cell(idx,:)];
        df_u = [df_u; df_tmp];
    else
        continue % for 0, do not update but continue.        
    end
end
a2b = string(df_u);
a2b_T_raw = table(a2b);
a2b_T_raw = splitvars(a2b_T_raw);
a2b_T_raw.Properties.VariableNames = [a_var, b_var];

if write ==1
    filename = strcat(a_in,'_concatA2B.xlsx');
    writetable(a2b_T_raw,filename);    
end           
end

function T = tableJoin(a_in,b_in,write)
%% Extract the data, df_param (Ref)
filename = strcat(a_in,'.xlsx');
a = readtable(filename);
a.Properties.VariableNames(1) = {'symbol'};

% mrna file, (to be appended, B)
filename = strcat(b_in,'.xlsx');
b = readtable(filename);
b.Properties.VariableNames(1) = {'symbol'};

[T,ileft,iright] = outerjoin(a,b);

idx = ~cellfun(@isnan,table2cell(T(:,4)));
T = T(idx,:);
idx = ~cellfun(@isnan,table2cell(T(:,2))); % no abundace for pax
T = T(idx,:);

if write ==1
    filename = strcat('outerJoin_b_v2.xlsx');
    writetable(T,filename);    
end

end