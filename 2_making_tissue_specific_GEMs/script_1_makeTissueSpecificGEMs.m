%% Script: Generate Transcriptomic Data Structure for tINIT2
% This script prepares a transcriptomic data structure compatible with tINIT2.
% It loads a pseudobulk transcriptomics dataset (e.g., from Fly Cell Atlas),
% processes the gene expression matrix, and saves it in a standardized format.

clc;
clear;
close all;

%% Set File and Folder Paths
pathway = pwd;
save_dir = '1_transcripts';
saveFolder = fullfile(pathway, save_dir);

if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end

% Input transcriptomics file (located in /files)
filename = 'FCA_pseudobulk_table_31_clusters_u.xlsx';
input_file = fullfile(pathway, 'files', filename);

% Generate data structure
out = run_makeTranscriptIn_v2(input_file);

% Move output to save folder
movefile('transcripts_dataStruct.mat', saveFolder);

disp(['Transcriptomic data structure saved to: ' saveFolder]);

%% Function: Create tINIT-Compatible Transcriptomic Data Structure
function out = run_makeTranscriptIn_v2(filepath)
% RUN_MAKETRANSCRIPTIN_V2
% Converts a gene expression table into a data_struct compatible with tINIT2.
%
% INPUT:
%   filepath - Full path to the Excel file (.xlsx) containing the transcriptomics data.
%              The file should have gene names in the first column and expression values in the rest.
%
% OUTPUT:
%   out - Struct with fields:
%         - genes:    Gene names
%         - tissues:  Sample/tissue names
%         - levels:   Expression levels (genes × tissues)
%         - threshold: Mean expression per tissue (used as threshold)
%   The struct is saved as 'transcripts_dataStruct.mat'.

% Read data
gtex_data = readtable(filepath);

% Parse fields
data_struct.genes     = gtex_data.gene;                              % Gene names
data_struct.tissues   = convertStringsToChars(string(gtex_data.Properties.VariableNames(2:end)))'; % Tissue names
data_struct.levels    = table2array(gtex_data(:, 2:end));            % Expression values
data_struct.threshold = nanmean(data_struct.levels, 1)';             % Mean expression per tissue

% Save struct
save('transcripts_dataStruct.mat', 'data_struct');

% Return output
out = data_struct;
end
