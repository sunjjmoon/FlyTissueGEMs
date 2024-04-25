% Goal:
%   To calculat the Vmax = kcat * ApparentEnzyme
%% 0. Set the path
path = pwd;

%% Make the output folder
output = 'F';
output_folder = strcat(path,'\',output);
if ~exist(output_folder,'dir')
    mkdir(output_folder)
end

%% Load the input file
fileName_ref = 'model_out';
fileName_ref_in = strcat(path,'\E_fill_the_EC_gap_otherOrganism\',fileName_ref);

[model_out gene_Vmax] = cal_Vmax(fileName_ref_in,output_folder);
% save
save(strcat(output_folder,'\','model_out.mat'),'model_out');

%% Extract the top 150 genes and their gene ontology
% This is to figure out what kind of biological functions these genes are
% related to
% Sample cell array containing expressions
expressions = gene_Vmax.grRules;
% Define a regular expression pattern to remove 'and', 'or', and special characters
pattern = '\s*(and|or|[(){}])\s*';
% Remove the specified elements using regexprep
flattened = string(regexprep(expressions, pattern, ' '));
% writetable(table(flattened),strcat(path,'\F_2_GOanalysis\','gene_Vmax.txt'));
writetable(table(gene_Vmax.Vmax),strcat(path,'\F\','Vmax.txt'));
b = [string(flattened) gene_Vmax.Vmax];
writetable(table(b),strcat(path,'\F\','gene_Vmax.txt'));

%%  Remove the nan and zeros
expressions = gene_Vmax.df_genes_sorted_noNanZero;
% Define a regular expression pattern to remove 'and', 'or', and special characters
pattern = '\s*(and|or|[(){}])\s*';
% Remove the specified elements using regexprep
gene_flattened = string(regexprep(expressions, pattern, ' '));
% writetable(table(flattened),strcat(path,'\F_2_GOanalysis\','gene_Vmax.txt'));
writetable(table(gene_Vmax.Vmax_sorted_noNanZero),strcat(path,'\F\','Vmax_regIncl.txt'));
b = [string(gene_flattened) gene_Vmax.Vmax_sorted_noNanZero];
writetable(table(b),strcat(path,'\F\','gene_Vmax_regIncl.txt'));

T = table(gene_flattened,gene_Vmax.Vmax_sorted_noNanZero);
writetable(T,strcat(path,'\F\','gene_Vmax_regIncl.xlsx'));

function [df_ref gene_Vmax] = cal_Vmax(fileName_ref_in,output_folder)
% Calculate by multiplying kcat and Apparent Enzyme.
% The maximum and minimum values will be determined, and the 0 values will
% be filled with the maximum values.

%% Extract the model and rxnGeneMat
% Reference is for the original model. 
df_ref = load(strcat(fileName_ref_in,'.mat'));
df_ref = df_ref.model_out;
df_kcat = df_ref.kcat;
df_genes = df_ref.grRules;
a = df_kcat.*3600; % s -> hr.

%% Extract the protein abundance info - a
b = df_ref.AppProtein4rxns;

%% Vmax Calculation
Vmax = b.*a;
vmax_max = max(Vmax);
vmax_min = min(Vmax);

%% Results
total_rxns = length(df_ref.rxns);
total_genes = length(df_ref.genes);
total_Vmax_without_0 = sum(Vmax~=0);
disp([newline 'Initial analysis']);
disp(['# of total rxns: ' num2str(total_rxns) newline ...
          '# of total proteins with value: ' num2str(total_Vmax_without_0)...
            ' (' num2str(total_Vmax_without_0./total_rxns*100) ' %, out of ' num2str(total_rxns) ')']);
disp(['Vmax val (max): '  num2str(vmax_max)]);
disp(['Vmax val (min): '  num2str(vmax_min)]);
disp(['Vmax val (med): '  num2str(nanmedian(Vmax))]);
disp(['Number of 0 vales: '  num2str(sum(Vmax==0)) ' ( ' num2str(sum(Vmax==0)/total_rxns*100) ' %, out of ' num2str(total_rxns) ')']);

out_of_1000 = sum(Vmax>1000);
disp([newline 'Out of max: '  num2str(out_of_1000) ' ( ' num2str(out_of_1000/total_rxns*100) ' %, out of ' num2str(total_rxns) ')']);

% %% 2023. 11. 16.
Vmax_excludes_0 = Vmax(Vmax~=0);
disp(['Number of Vmax without 0 : ' num2str(length(Vmax_excludes_0)) ' out of ',...
        num2str(total_rxns) ' (' num2str(length(Vmax_excludes_0)/total_rxns*100) ')']);
disp(['Vmax val (max): '  num2str(max(Vmax_excludes_0))]);
disp(['Vmax val (min): '  num2str(min(Vmax_excludes_0))]);
disp(['Vmax val (med): '  num2str(nanmedian(Vmax_excludes_0))]);
out_of_1000 = sum(Vmax_excludes_0>1000);
disp(['Out of max: '  num2str(out_of_1000) ' ( ' num2str(out_of_1000/total_rxns*100) ' %, out of ' num2str(total_rxns) ')']);

%% Output
% apparent enzyme abundance for rxn. If 1 gene associated, use the protein
% abundance. If multiple, use the maximum. If none, use 0. For 0, we will
% later use the maximum vmax to make loose boundary.
%% Graph
x = 1:1:length(Vmax);
% Sort the Vmax values in ascending order
[Vmax_sorted, order] = sort(Vmax);
% gene_sorted = 
log_Vmax_sorted = log10(Vmax_sorted);
log_Vmax_sorted_nonInfNan = log_Vmax_sorted(~isinf(log_Vmax_sorted));
log_Vmax_sorted_nonInfNan = log_Vmax_sorted_nonInfNan(~isnan(log_Vmax_sorted_nonInfNan));
% Create a scatter plot with ascending order
figure(1)
scatter([1:length(log_Vmax_sorted_nonInfNan)], log_Vmax_sorted_nonInfNan, 'filled');
yline(3,'--','Default max','Linewidth',1)
% Set axis labels and title
xlabel('reaction numbers');
ylabel('Vmax (log_{10} mmol g-dCW^{-1} hr^{-1})');
grid on;

% Save the scatter plot as a JPEG image
filename = strcat(output_folder,'\','scatter_plot.jpg');  % Specify the desired file name
set(gca, 'FontSize',15);
saveas(gcf, filename, 'jpg');

% Historgram
figure(2)% Sample data
% Specify the bin width
binwidth = 0.5;
hist = histogram(log_Vmax_sorted, 'BinWidth', binwidth,'Normalization','probability');
xline(3,'--','Default max','Linewidth',1.5)
xlabel('log_{10}Vmax (mmol g-dCW^{-1} hr^{-1})');
ylabel('Frequency');
grid on;
filename = strcat(output_folder,'\','histogram.jpg');  % Specify the desired file name
set(gca, 'FontSize',15);
saveas(gcf, filename, 'jpg');
exportgraphics(gcf,strcat(output_folder,'\','histogram.pdf'))% saveas(gcf,'shortVSlong.tiff');

% cumulative distribution function (cdf) plot
figure(3)
h = cdfplot(log_Vmax_sorted);
h.LineWidth = 2;
xlabel('Vmax (log_{10} mmol g-dCW^{-1} hr^{-1})');
ylabel('Fraction of reactions');
xline(3,'--','Default max','Linewidth',1)
grid on;
filename = strcat(output_folder,'\','cdf_plot.jpg');  % Specify the desired file name
set(gca, 'FontSize',15);
saveas(gcf, filename, 'jpg');

%% Redefine the Vmax values
% Define the threshold for values greater than 10^3
threshold = 10^3; % This is the default max setting

vmax_tmp = Vmax;
vmax_tmp(vmax_tmp > threshold) = threshold; % for those higher, set the max
vmax_tmp(vmax_tmp == 0) = threshold; % for unbounded, set the max

df_ref.ub = vmax_tmp;

% For the minimum, already the directionality was considered
idx = find(df_ref.lb ~= 0);
df_ref.lb(idx) = -1.*vmax_tmp(idx);

%% Export
df_ref.Vmax = Vmax;
out = df_ref;
df_ref.update_F = 'Vmax_calculated and updated';

[Vmax_sorted, order] = sort(Vmax,'descend');
gene_Vmax.grRules = df_genes(order);
gene_Vmax.Vmax = Vmax_sorted;

df_genes_sorted = df_genes(order);
noZeroIdx = find(Vmax_sorted~=0);
Vmax_sorted_noNanZero = Vmax_sorted(noZeroIdx);
df_genes_sorted_noNanZero = df_genes_sorted(noZeroIdx);

noNanIdx = find(~isnan(Vmax_sorted_noNanZero));
Vmax_sorted_noNanZero = Vmax_sorted_noNanZero(noNanIdx);
df_genes_sorted_noNanZero = df_genes_sorted_noNanZero(noNanIdx);

gene_Vmax.df_genes_sorted_noNanZero = df_genes_sorted_noNanZero;
gene_Vmax.Vmax_sorted_noNanZero = Vmax_sorted_noNanZero;
end