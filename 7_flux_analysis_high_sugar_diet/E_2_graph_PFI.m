close all
clear all
%% Load the files
text_draw = 8; % draw top 5 text
bottom_draw = 5; % draw the bottom 5.
text_font_size = 10;

pathway = pwd;
save_dir = 'E_2_PFI';
save_folder = strcat(pathway,'\',save_dir);

if ~exist(save_folder,'dir')
    mkdir(save_folder)
end

filename = 'output';
sheetName = 'flux';
df = readtable(strcat(pathway,'\E_1_PFI_table\',filename,'.xlsx'),'Sheet',sheetName);

sheetName = 'SubSysOrig';
Subsystem = table2cell(readtable(strcat(pathway,'\E_1_PFI_table\',filename,'.xlsx'),'Sheet',sheetName));

flux_NSD = table2array(df(:,2))';
flux_HSD = table2array(df(:,3))';

relChange = flux_HSD./flux_NSD;
relChange_tmp = flux_HSD-flux_NSD;
relChange_tmp = flux_HSD./flux_NSD;

relChange = (relChange_tmp-mean(relChange_tmp))./std(relChange_tmp);
relChange = relChange_tmp;

[relChange_sorted,idx] = sort(relChange,'ascend');
subsystem_sorted = Subsystem(idx);
subsystem_sorted = Subsystem(idx);
sz = 30;

x = log10(flux_NSD);
x = 1:1:length(relChange);
y = relChange_sorted;

%% Plot
%% Report the excelfile
% Remove the unrealistic output - the tryptophan pathway is wrongly annotated with melatonin degradation not exisiting in muscle
idx = find(string(subsystem_sorted) == 'Tryptophan metabolism');
subsystem_sorted(idx) = [];
relChange_sorted(idx) = [];
df_report = table(subsystem_sorted,relChange_sorted');
writetable(df_report,strcat(save_folder,'\','output_fvs_graph.xlsx'))

%% Without the Inf and NaN
[sorted_data_v2,subsystem_sorted_v2,flux_NSD_v2] = sort_excluding_inf_nan(relChange_sorted,subsystem_sorted,flux_NSD)
x = 1:1:length(subsystem_sorted_v2);
y = sorted_data_v2;
%% 
% y = log2(y);
y = log2(y);
%% Plot
figure(2)
% Example with a color map (parula)
cmap = jet(length(y));
scatter(x, y, sz, cmap, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
% scatter(x,y,sz,'filled')
hold on;  % Keep the current plot while adding text
% Add text to the first three points
length_idx = length(y);
end_point = length_idx-bottom_draw+1;
% Add text to the first three points
for i = 1:text_draw
    text(x(i)+1, y(i), [subsystem_sorted(i)], 'FontSize', text_font_size, 'Color', 'b');
end
hold on
for i = length_idx:-1:end_point
    text(x(i)+1, y(i), [subsystem_sorted_v2(i)], 'FontSize', text_font_size, 'Color', 'r');
end
% Draw the least affected
yline(0)
xlabel(['Metabolic subsystems'])
ylabel(['Relative flux potential index' newline '(log2(FVS_{HSD}/FVS_{NSD}))'])
set(gca,'FontSize',13)
grid on
grid minor
hold off


save_fig = 'rel_fvs_v2';
% saveas(gcf,strcat(save_folder,'\',save_fig,'.tiff'))
saveas(gcf,strcat(save_folder,'\',save_fig,'.fig'))

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
% print('-dpdf', '-painters', char(strcat(save_folder,'\fva.pdf')));

%% Report the excelfile
df_report = table(subsystem_sorted,relChange_sorted');
writetable(df_report,strcat(save_folder,'\','output_fvs_graph.xlsx'),'Sheet','withoutNaNInf')


function [sorted_data_v2,subsystem_sorted_v2,flux_NSD_v2] = sort_excluding_inf_nan(data,subsystem_sorted,flux_NSD)
  %  Sorts a list or array, excluding values with inf or nan.
%  Args:   data: A list or array of numbers.
%  Returns:   A sorted list or array, excluding values with inf or nan.

  % Find indices of non-finite values (inf or nan)
  non_finite_indices = isnan(data) | isinf(data);

  % Sort the data excluding non-finite values
  sorted_data_v2 = data(~non_finite_indices);
  subsystem_sorted_v2 = subsystem_sorted(~non_finite_indices);
  flux_NSD_v2 = flux_NSD(~non_finite_indices);

  end

