%% Make the tissue-specific model
path = pwd;
load(strcat(path,'\models\','fruitfly3.mat'));
model = fruitfly_v3;
data_struct = load(strcat(path,'\1_transcripts\','transcripts_dataStruct.mat'));
data_struct = data_struct.data_struct;
model_out = run_tINITyml_multiple_specific(model,data_struct);

function model_out = run_tINITyml_multiple_specific(model, data_struct)
% RUN_TINITYML_MULTIPLE_SPECIFIC
% Builds tissue-specific GEMs using the tINIT2 algorithm with a YAML-based
% Fruitfly GEM and pseudobulk transcriptomics input.
%
% INPUTS:
%   model        - Base metabolic model structure (e.g., parsed from YAML)
%   data_struct  - Struct containing gene expression data with fields:
%                  .genes, .tissues, .levels, .threshold (optional)
%
% OUTPUT:
%   model_out    - The final model created in the last iteration (others saved to disk)
%
% For each tissue, this script runs tINIT and saves the resulting model as a .mat file.
path = pwd;

%% Prepare Model (Fix compartments and metabolite formatting)
model.mets = strrep(model.mets, 'x', 'p');      % Replace 'x' with 'p' in metabolite IDs
model.comps = strrep(model.comps, 'x', 'p');    % Replace 'x' with 'p' in compartment IDs

%% Set up solver and task file
setRavenSolver('gurobi');

% Load predefined metabolic tasks
taskFilePath = strcat(path,'\files\metabolicTasks_Essential.txt');
essentialTasks = parseTaskList(taskFilePath);

%% Output directory
outPat = '2_tissue_specific_gems';
saveFolder = fullfile(path, outPat);
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder)
end

%% Loop through tissues (adjust range if needed)
for i = 1:length(data_struct.tissues)
    i=26
    % Prepare tissue-specific data input
    data_struct_in.tissues = data_struct.tissues(i);
    data_struct_in.genes   = data_struct.genes;
    data_struct_in.levels  = data_struct.levels(:, i);
    data_struct_in.threshold = 1;  % Fixed expression threshold : Mean-SD ~ 1.02.

    % Define tINIT input parameters
    refModel         = model;
    tissue           = char(data_struct_in.tissues);
    celltype         = [];
    hpaData          = [];
    arrayData        = data_struct_in;
    metabolomicsData = [];
    removeGenes      = true;
    taskFile         = [];
    useScoresForTasks = true;
    printReport      = true;
    taskStructure    = essentialTasks;
    params           = [];
    paramsFT         = [];

    % Run tINIT2
    model_out = getINITModel2(refModel, tissue, celltype, hpaData, ...
        arrayData, metabolomicsData, removeGenes, taskFile, ...
        useScoresForTasks, printReport, taskStructure, params, paramsFT);

    % Assign model ID and save
    model_out.id = tissue;
    outFileName = fullfile(outPat, [tissue, '.mat']);
    save(outFileName, 'model_out');

    fprintf('Saved model for tissue: %s\n', tissue);
end

end
