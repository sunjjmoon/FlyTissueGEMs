%% Goal
%   - to calculate protein based on mRNA expression levels and the Schwan
%   parameter values

%% Define the save folder
pathway = pwd;
save_dir = 'A';
subfolder = [pathway '\' save_dir];
if ~exist(subfolder, 'dir')
    mkdir(subfolder)
end

%% Load inputs
file_param = strcat(pathway,'\files\','results_human2flyGenematched_param_uuu'); % v2 is for the whole fly only
file_mrna = strcat(pathway,'\files\','mRNA_v2'); % v2 is for the whole fly only
write = 1;

%% Run
[gene2param, var_out] = calProteinFromgene_v4(file_param,file_mrna,write);
save(strcat(subfolder,'\','gene2param.mat'),'gene2param');
save(strcat(subfolder,'\','var_out.mat'),'var_out');

movefile('results_gene2param_T2uuu.xlsx',subfolder);

function [gene2param, var_out] = calProteinFromgene_v4(file_param,file_mrna,write)
%% Extract the data, df_param (Ref)
filename = strcat(file_param,'.xlsx');
df_param = readtable(filename);
df_param_var = df_param.Properties.VariableNames;

% mrna file, (to be appended, B)
filename = strcat(file_mrna,'.xlsx');
df_mran = readtable(filename);
df_mran_var = df_mran.Properties.VariableNames;

%% Loop the df_param to append the associated mRNA abundances
ref = table2cell(df_param);
a = table2cell(df_mran);
%% Loop the param (Ref) to find the matched array for query (a, mRNA)
df_u = [];
for i = 1:length(ref)
    [Lia,~] = ismember(a(:,1),ref(i,1));
%     count = count +sum(Lia);
%     disp(count)
    if sum(Lia) > 1 % Then, multiple match
        idx = find(Lia);
        tmp = a(idx,end); 
        tmp_df_cell = repmat(ref(i,:),length(idx),1);
        df_tmp = [tmp_df_cell, tmp];
        df_u = [df_u; df_tmp];
    elseif sum(Lia) == 1  % That means, there are 1 to 1 match
        idx = find(Lia);
        df_tmp = [ref(i,:), a(idx,end)];
        df_u = [df_u; df_tmp];
    else
        continue % for 0, do not update but continue.        
    end
end

% Get only the unique rows.
[~,ia,~] = unique(df_u(:,1),'stable');
df_u = df_u(ia,:);

rowNames = df_u(:,1);
df_u(:,1:3) =[];
df_u = string(df_u);
df_u = str2double(df_u);

avo = 6.022e23;
total_nucleotides = 2.6e-15;
rpk_readsTranscript_per_readsTot_perTranscriptLength = df_u(:,end)./10^6./10^3;
mRNA_copyNumb = rpk_readsTranscript_per_readsTot_perTranscriptLength.*total_nucleotides.*avo;
protein_copyNumb = df_u(:,1)./df_u(:,2).*mRNA_copyNumb;

r = 5e-6; % 5 um for drosophila in average
vol_s2 = 4/3*pi*(r).^3.*10^3; % L
avg_aaMass = 110; % Da per a.a
protein_mass_hela = 150*10^-12;
protein_ppm = protein_copyNumb.*df_u(:,end-2).*avg_aaMass./avo./protein_mass_hela*10^6; % multiply by length of aa and mass, which will give Da/protein molecules

gene2param = [df_u, mRNA_copyNumb, protein_copyNumb, protein_ppm];

gene2param_T_raw = table(gene2param,'RowNames',rowNames);
gene2param_T_raw = splitvars(gene2param_T_raw);
gene2param_T_raw.Properties.VariableNames = [df_param_var(4:end), df_mran_var(end),...
    {'mRNA_copyNumb'}, {'protein_copyNumb'}, {'protein_ppm'}];

var_out = gene2param_T_raw.Properties.VariableNames;

if write ==1
    filename = ['results_gene2param_T2uuu.xlsx'];
    writetable(gene2param_T_raw,filename,'WriteRowNames',true);    
end
end