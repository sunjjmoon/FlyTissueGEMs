function paxDB_out = readPaxData(filePaxDB)
% Goal:
%   - to extract the PaxDB in a format to use the flyGene matched
% Input:
%   1. filePaxDB = '7227-WHOLE_ORGANISM-integrated-paxdb';

fileload = strcat(filePaxDB,'.txt');
df = readtable(fileload,'delimiter', '\t', 'ReadVariableNames',true);
VarNames = df.Properties.VariableNames;
df_cell = table2cell(df);
df(:,2) = extractAfter(df_cell(:,2),'.');
paxDB_out = df;
% paxDB_out = cell2table(df,'VariableNames',VarNames);
% paxDB_out = string(extractAfter(df_fbpp,'.'));

end