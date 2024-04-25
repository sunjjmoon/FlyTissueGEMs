%% Define the save folder
pathway = pwd;
save_dir = 'A';
subfolder = [pathway '\' save_dir];
if ~exist(subfolder, 'dir')
    mkdir(subfolder)
end

%% Load the files 
load(strcat(pathway,'\models\','fruitfly1.mat'))
GEM_model = fruitflyGEM;

%% Find the empty arrays
[out_empty,uniqKEGGid,T] = run_findEmptyArray_v2(GEM_model);
movefile('out_empty.mat',subfolder);

function [out_empty,uniqKEGGid,T] = run_findEmptyArray_v2(GEM_model)
% v2: to get the rxn associated KEGG id as well.
%  GEM_model = fruitflyGEM;
out_empty =[];
for i = 1:length(GEM_model.eccodes)
    tmp = find(isempty(GEM_model.eccodes{i,1}));
    if tmp == 1
        idx_empty = i;
        reaction = GEM_model.rxns{idx_empty,1};
        keggID = GEM_model.rxnKEGGID{idx_empty,1};
        grRules = GEM_model.grRules{idx_empty,1};
        tmp_store = [{idx_empty}, {reaction}, {keggID}, {grRules}];
        out_empty = [out_empty;tmp_store];
    else
        continue
    end
end

% Determine the number of reactions without EC numbers.
out_empty = string(out_empty);
disp([newline, '# of rxns without EC: ' num2str(length(out_empty)), newline,...
        'tot of rxns : ' num2str(length(GEM_model.rxns)), newline,...
        '% of rxns wo EC: ' num2str(length(out_empty)./length(GEM_model.rxns)*100)]);

T = table(out_empty); T = splitvars(T);
T.Properties.VariableNames = {'idx_empty', 'reaction', 'keggID', 'grRules'};
filename = 'raw_rxns_wo_ec.xlsx';
% write(T,filename);

%% Extract the unique keggID
a = out_empty(:,3); % obtain the KEGG-ID.
b = a(~cellfun('isempty',a)); % obtain non-empty kegg-ID. Or reactions with EC and Kegg. 161 reactions.
uniqKEGGid = unique(b,'stable'); % uniq KEGG id 
disp([newline, '# of rxns with kegg-ID: ' num2str(length(b)), newline,...
       '# of rxns wo kegg-ID: ' num2str([length(a)-length(b)]) newline,...
        'tot of rxns without EC : ' num2str(length(out_empty)) newline,...
        '% of rxns that will be updated with EC: ' num2str(length(b)./length(a)*100)]);

save('out_empty.mat','out_empty')
end