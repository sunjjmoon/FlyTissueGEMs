%% Compare Metabolic Models Across All Tissues
% This script compares the metabolic models of all tissues using t-SNE 
% visualization and subsystem-level analysis.

clc;
clear;
close all;

load(strcat(pwd,'\2_tissue_specific_gems\','model_all.mat'));

%% Define Save Directory
% Update this path to reflect the output location
save_dir = '3_analysis_network_subsystem';  % 

%% Plot Settings
tsneGraph = 1;  % Set to 1 to generate t-SNE plots

%% Subsystem Difference Thresholds
subSysThreshold     = 100;  % Relative difference threshold for standard clusterogram
subSysThreshold_v2  = 20;   % Absolute difference threshold (for alternate graph)
subSysThreshold_v3  = 100;  % Relative difference threshold (alternate graph)

%% Tissue Column Index
control_col = 2;  % Use the second column for tissue annotations (e.g., 'Muscle', 'Gut', etc.)

%% Run Comparison Function
fn_run_compareModel_all(model_all, save_dir, tsneGraph, ...
    subSysThreshold, subSysThreshold_v2, subSysThreshold_v3, control_col);

% After running, you may use `plot(...)` commands in the output directory 
% to further customize figures (e.g., adjust colorbars, titles).

% Example:
% colorbar('northoutside');
% title('Comparison of Metabolic Network Structures');
