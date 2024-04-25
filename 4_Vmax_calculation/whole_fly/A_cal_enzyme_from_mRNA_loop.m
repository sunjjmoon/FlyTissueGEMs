currentFolder = pwd();  % Get the current folder
addpath(genpath(currentFolder));

%% 0. Set the path
path = pwd;

%% Load the input file
folder = strcat(path,'\','input_data');
fileName_mrna = 'mRNA_v2';
fileName_mran_in = strcat(folder,'\',fileName_mrna);
fileName_ref = 'results_human2flyGenematched_param_uuu';
fileName_ref_in = strcat(folder,'\',fileName_ref);

%% Set the parameter values
write = 1;
write_all = 1; % Write all

%% Calculate the protein abundance from gene
gene2protein = calProteinFromgene_loop(fileName_ref_in,fileName_mran_in,write,path,write_all);

save(strcat(path,'\A_gene2protein_loop\','gene2protein.mat'), 'gene2protein');

function gene2protein = calProteinFromgene_loop(file_param,file_mrna,write,path,write_all)
%% Make the output folder
output = 'A_gene2protein_loop';
output_folder = strcat(path,'\',output);
if ~exist(output_folder,'dir')
    mkdir(output_folder)
end


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

df_protein_copy = [];
df_protein_ppm = [];
df_protein_mmolGdcw = [];

gene2protein = struct();

for j = 2:size(a,2) %Loop each column, except the first one, which is the gene_id.
%% Loop the param (Ref) to find the matched array for query (a, mRNA)
df_u = [];
for i = 1:length(ref)
    [Lia,~] = ismember(a(:,1),ref(i,1));
%     count = count +sum(Lia);
%     disp(count)
    if sum(Lia) > 1 % Then, multiple match
        idx = find(Lia);
        tmp = a(idx,j); 
        tmp_df_cell = repmat(ref(i,:),length(idx),1);
        df_tmp = [tmp_df_cell, tmp];
        df_u = [df_u; df_tmp];
    elseif sum(Lia) == 1  % That means, there are 1 to 1 match
        idx = find(Lia);
        df_tmp = [ref(i,:), a(idx,j)];
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
correction_factor = 10; % Since the values are not rpkm or fpkm, it is not technically the rpkm values. but using this makes in the same order of magnitude.
rpk_readsTranscript_per_readsTot_perTranscriptLength = df_u(:,end)./10^6./10^3.*correction_factor;
mRNA_copyNumb = rpk_readsTranscript_per_readsTot_perTranscriptLength.*total_nucleotides.*avo;
protein_copyNumb = df_u(:,1)./df_u(:,2).*mRNA_copyNumb.*correction_factor;
r = 5e-6; % 5 um for drosophila in average
vol_s2 = 4/3*pi*(r).^3.*10^3; % L
avg_aaMass = 110; % Da per a.a
protein_mass_hela = 150*10^-12;
protein_ppm = protein_copyNumb.*df_u(:,end-2).*avg_aaMass./avo./protein_mass_hela*10^6; % multiply by length of aa and mass, which will give Da/protein molecules

%% Convert ppm to mmol/gDCW - based on Lewis' method
%Determine the molecular weight of the protein: Obtain the molecular weight of the protein of interest in g/mol. This information is typically available in protein databases or literature.
%Convert ppm to mg/g: Since ppm represents the ratio of the protein mass to the total mass, you need to convert it to milligrams per gram (mg/g). This can be done by dividing the ppm value by 1,000.
%Convert mg/g to mmol/g: To convert mass (mg) to moles (mmol), divide the mass in milligrams by the molecular weight of the protein in grams/mole.
%Consider the dry cell weight (DCW): If you want to express the concentration in mmol per gram of dry cell weight (gDCW), you need to account for the dry cell weight of the sample. The dry cell weight represents the weight of the cells without the presence of water. It can be determined experimentally through cell harvesting and drying techniques.
%Perform the final conversion: Divide the concentration in mmol/g by the dry cell weight in grams to obtain the concentration in mmol/gDCW.
maa = df_u(:,8).*1000 ; %The average molecular mass of an amino acid in Daltons (or g/mol)
protein_mmol_gDCW = protein_ppm./1000./maa;
gene2param = [df_u, mRNA_copyNumb, protein_copyNumb, protein_ppm];
gene2param_T_raw = table(gene2param,'RowNames',rowNames);
gene2param_T_raw = splitvars(gene2param_T_raw);
gene2param_T_raw.Properties.VariableNames = [df_param_var(4:end), df_mran_var(end),...
    {'mRNA_copyNumb'}, {'protein_copyNumb'}, {'protein_ppm'}];
var_out = gene2param_T_raw.Properties.VariableNames;
if write ==1   
    filename = strcat(output_folder,'\',char(df_mran_var(j)),'_2protein.xlsx');
    writetable(gene2param_T_raw,filename,'WriteRowNames',true);    

end
%% add up the protein information
df_protein_copy = [df_protein_copy, protein_copyNumb];
df_protein_ppm = [df_protein_ppm, protein_ppm];
df_protein_mmolGdcw = [df_protein_mmolGdcw,protein_mmol_gDCW];   
end  
%% Write excel file for protein all information in each sheet
% Write the original
sheetName = 'original_mrna';
if write_all ==1   
    filename = strcat(output_folder,'\','all.xlsx');
    writetable(df_mran,filename,'Sheet',sheetName,'WriteRowNames',true);    
end
% Write the protein copy number
gene2param_T1_raw = table(df_protein_copy,'RowNames',rowNames);
gene2param_T1_raw = splitvars(gene2param_T1_raw);
gene2param_T1_raw.Properties.VariableNames = df_mran_var(2:end);
sheetName = 'protein_copy';
if write_all ==1   
    filename = strcat(output_folder,'\','all.xlsx');
    writetable(gene2param_T1_raw,filename,'Sheet',sheetName,'WriteRowNames',true);    
end
% Write the protein ppm
gene2param_T_raw = table(df_protein_ppm,'RowNames',rowNames);
gene2param_T_raw = splitvars(gene2param_T_raw);
gene2param_T_raw.Properties.VariableNames = df_mran_var(2:end);
sheetName = 'ppm';
if write_all ==1   
    filename = strcat(output_folder,'\','all.xlsx');
    writetable(gene2param_T_raw,filename,'Sheet',sheetName,'WriteRowNames',true);    
end
% Write the protein mmol_gDCW
gene2param_T2_raw = table(df_protein_mmolGdcw,'RowNames',rowNames);
gene2param_T2_raw = splitvars(gene2param_T2_raw);
gene2param_T2_raw.Properties.VariableNames = df_mran_var(2:end);
sheetName = 'mmol_gDCW';
if write_all ==1   
    filename = strcat(output_folder,'\','all.xlsx');
    writetable(gene2param_T2_raw,filename,'Sheet',sheetName,'WriteRowNames',true);    
end

%% Structure out
gene2protein.original = df_mran;
gene2protein.protein_copy = gene2param_T1_raw;
gene2protein.protein_ppm = gene2param_T_raw;
gene2protein.protein_mmolGdcw = gene2param_T2_raw;

end