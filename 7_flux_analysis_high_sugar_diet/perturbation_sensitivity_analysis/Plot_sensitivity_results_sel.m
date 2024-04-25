clc
clear all
%% Set the save and load folder
pathway = pwd;
% Get the current working directory
% Extract the parent folder
dirContents = dir(pathway);
% Initialize an empty cell array to store folder names
folderNames = {};
% Loop through each item in dirContents
for i = 1:numel(dirContents)
    % Check if the item is a directory and not '.' or '..'
    if dirContents(i).isdir && ~strcmp(dirContents(i).name, '.') && ~strcmp(dirContents(i).name, '..')
        % Add the folder name to the folderNames cell array
        folderNames{end+1} = dirContents(i).name;
    end
end
% Filter folder names that start with a number followed by an underscore
filteredFolderNames = folderNames(startsWith(folderNames, {'1_', '3_', '5_', '9_','10_'}));

%% Extract the names
% Initialize a cell array to store the extracted words
extractedWords = cell(size(filteredFolderNames));

% Loop through each folder name and extract the words between underscores
for i = 1:numel(filteredFolderNames)
    % Apply regular expression to extract words between underscores
    tokens = regexp(filteredFolderNames{i}, '_([^_]+)_', 'tokens');
    
    % Store the extracted words in the cell array
    if ~isempty(tokens)
        extractedWords{i} = tokens{1}{1};
    else
        extractedWords{i} = ''; % Handle case where no match is found
    end
end

% Display the extracted words
disp('Extracted words:');
disp(extractedWords');


%% Load the data
sens = zeros(length(filteredFolderNames),1);
gene_names = [];
for i = 1:length(filteredFolderNames)
    name = filteredFolderNames{1,i};
     
    loadfolder = [pathway  '\' name];
    df = readtable(strcat(loadfolder,'\summary_data.xlsx'),'Sheet','Summary');
    sens(i,1) = df.sensitivity;
    gene = extractBetween(name,'_','_');
    gene_names{i,1} = gene;

end
%% Sort by largest decrease to lower
[sens_ordered, I] = sort(sens,'descend');
gene_names_ordered = gene_names(I);
gene_names_ordered = string(gene_names_ordered);
gene_names_ordered_cat = categorical(gene_names_ordered);
gene_names_ordered_cat_reorder = reordercats(gene_names_ordered_cat,gene_names_ordered);

%% Plot the graph
xlabText = ['Perturbation Effect' newline '(\Deltasummed rates_{glyc}/\Deltarates_{glyc,i,upper})'];
figure(1)
h= barh(gene_names_ordered_cat_reorder,sens_ordered);
xlabel(xlabText);
set(gca, 'FontSize', 18,'FontWeight','bold'); % Set font size for tick labels
grid on 
grid minor

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
% print('-dpdf', '-painters', 'Sensitivity_results_sel.pdf');
print(gcf, 'Sensitivity_results_sel.tiff', '-dtiff', '-r300');


t = table(gene_names_ordered_cat_reorder,sens_ordered);
writetable(t, 'Sensitivity_out_sel.xlsx');
