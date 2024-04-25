function flybase_fields_out_c_u_t = readFlybaseFields(file_flybase,refPaxDB)
% Read the tsv files obtained from the flybase batch downloads
% Goal:
% Input:
%   1. file_flybase = 'FlyBase_Fields_geneFromPaxDBString';
%   2. refPaxDB = '7227-WHOLE_ORGANISM-integrated-paxdb';
%       - this is the refence paxDB file used to get the flybase data. The
%       resulting gene names will be concetanated to the original files.
%% Load flybase Fields data (based on the refPaxDB)
fileload = strcat(file_flybase,'.txt');
flybase_fields_out = readtable(fileload,'delimiter', '\t', 'ReadVariableNames',true);
var_flybase = flybase_fields_out.Properties.VariableNames;
flybase_fields_out_c = table2cell(flybase_fields_out);
% flybase_fields_out = string(extractAfter(df_fbpp,'.'));

%% Load the refPaxDB
fileload = strcat(refPaxDB,'.txt');
df = readtable(fileload,'delimiter', '\t', 'ReadVariableNames',true);
df_cell = table2cell(df);
df_cell(:,2) = extractAfter(df_cell(:,2),'.');

%% Loop the Paxdb proteins to find the match in flybase_fields and concetanate
abundance = [];
for i = 1:length(flybase_fields_out_c)
    [Lia,~] = ismember(df_cell(:,2),flybase_fields_out_c(i,1));
    abundance(i,1) = string(df_cell(Lia,3)); % abundance col
end

%% Concetanate
flybase_fields_out_c_u = [flybase_fields_out_c, num2cell(abundance)];
T = table(flybase_fields_out_c_u);
T = splitvars(T);
T.Properties.VariableNames = [var_flybase, {'abundance'}];
flybase_fields_out_c_u_t = T;

%% It is important that one should check "associated_gene" to ID check to convert the symbols into the alphabetic.
end