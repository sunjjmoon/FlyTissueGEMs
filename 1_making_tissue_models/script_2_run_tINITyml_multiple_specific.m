model_out = run_tINITyml_multiple_specific(model,data_struct);

function model_out = run_tINITyml_multiple_specific(model,data_struct)
%% Add the boundary Metabolites - If does not exist
pattern = 'x';
replace = 'p'; 
model.mets = strrep(model.mets,pattern,replace);
model.comps = strrep(model.comps,pattern,replace);
savePath = pwd;
%% Check parseTaskList
setRavenSolver('gurobi')
essentialTasks = parseTaskList(strcat(pwd,'\files\','metabolicTasks_Essential.txt'));

for i = 1:length(data_struct.tissues)

data_struct_in.tissues = data_struct.tissues(i);  % sample (tissue) names
data_struct_in.genes = data_struct.genes(:,1);  % gene names
data_struct_in.levels = data_struct.levels(:,i); % This is the gene expression level data
data_struct_in.threshold = 1; % 

refModel = model;  % the reference model from which the GEM will be extracted
tissue = char(data_struct_in.tissues);  % must match the tissue name in data_struct.tissues
celltype = [];  % used if tissues are subdivided into cell type, which is not the case here
hpaData = [];  % data structure containing protein abundance information (not used here)
arrayData = data_struct_in;  % data structure with gene (RNA) abundance information
metabolomicsData = [];  % list of metabolite names if metabolomic data is available
removeGenes = true;  % (default) remove lowly/non-expressed genes from the extracted GEM
taskFile = [];  % we already loaded the task file, so this input is not required
useScoresForTasks = true;  % (default) use expression data to decide which reactions to keep
printReport = true;  % (default) print status/completion report to screen
taskStructure = essentialTasks;  % metabolic task structure (used instead "taskFile")
params = [];  % additional optimization parameters for the INIT algorithm
paramsFT = [];  % additional optimization parameters for the fit-tasks algorithm

model_out = getINITModel2(refModel, tissue, celltype, hpaData, arrayData, metabolomicsData, removeGenes, taskFile, useScoresForTasks, printReport, taskStructure, params, paramsFT);
out_name = strcat(data_struct.tissues(i), '.mat');
out_name2 = char(strcat(savePath,'/',...
    out_name));
model_out.id = char(data_struct.tissues(i));
save(out_name2, 'model_out')

end
end