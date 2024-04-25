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

[model_out, gene_Vmax,Vmax_order_tissues] = cal_Vmax(fileName_ref_in,output_folder);
% save
save(strcat(output_folder,'\','model_out.mat'),'model_out');

function [out gene_Vmax,Vmax_order_tissues] = cal_Vmax(fileName_ref_in,output_folder)
close all
% Calculate by multiplying kcat and Apparent Enzyme.
% The maximum and minimum values will be determined, and the 0 values will
% be filled with the maximum values.

%% Extract the model and rxnGeneMat
% Reference is for the original model. 
df_ref_tmp = load(strcat(fileName_ref_in,'.mat'));
df_ref_tmp = df_ref_tmp.model_out;

%% Multiple models (2023.10.29)
for p = 1:length(df_ref_tmp)
df_ref = df_ref_tmp{1,p};

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
disp(['Vmax val (med): '  num2str(median(Vmax))]);
disp(['# of 0 vales: '  num2str(sum(Vmax==0)) ' ( ' num2str(sum(Vmax==0)/total_rxns*100) ' %, out of ' num2str(total_rxns) ')']);

out_of_1000 = sum(Vmax>1000);
disp(['Out of max: '  num2str(out_of_1000) ' ( ' num2str(out_of_1000/total_rxns*100) ' %, out of ' num2str(total_rxns) ')']);

%% Output
% apparent enzyme abundance for rxn. If 1 gene associated, use the protein
% abundance. If multiple, use the maximum. If none, use 0. 
%% Redefine the Vmax values
% Define the threshold for values greater than 10^3
threshold = 10^3; % This is the default max setting

vmax_tmp = Vmax;
vmax_tmp(vmax_tmp > threshold) = threshold; % for those higher, set the max
vmax_tmp(vmax_tmp == 0) = threshold; % for unbounded, set the max

out{1,p} = df_ref;
out{1,p}.ub = vmax_tmp;

% For the minimum, already the directionality was considered
idx = find(out{1,p}.lb ~= 0);
out{1,p}.lb(idx) = -1.*vmax_tmp(idx);

%% Export
out{1,p}.Vmax = Vmax;
out{1,p}.update_F = 'Vmax_calculated and updated';

[Vmax_sorted, order] = sort(Vmax,'descend');
log_Vmax_sorted = log10(Vmax_sorted);

gene_Vmax{1,p}.grRules = df_genes(order);
gene_Vmax{1,p}.Vmax = Vmax_sorted;

gene_Vmax{1,p}.logVmax = log10(sort(Vmax,'ascend'));
gene_Vmax{1,p}.id = df_ref.id;

end

%% Re-order the tissues by having higher median values for the Vmax
for p = 1:length(gene_Vmax)
    val_tmp = gene_Vmax{1,p}.Vmax;
    val_tmp = val_tmp(val_tmp~=0);
    vmax_med(p,1) = median(val_tmp,'omitnan');
    vmax_mean(p,1) = mean(val_tmp,'omitnan');
    tmp_id{p,1} = gene_Vmax{1,p}.id;
    
    val_tmp_sorted = sort(val_tmp,'ascend');
    diff = abs(val_tmp_sorted-1000);
    [tmp_mix,idx_min] = min(diff);
    number_that_reach_max(p,1) = (idx_min);

end
% vmax_med = vmax_mean;
[Vmax_med_nonZero_sorted, order] = sort(vmax_med,'descend');
% Normalize the numeric values using zscore
Vmax_med_nonZero_sorted_normalized = zscore(Vmax_med_nonZero_sorted);
tmp_id_sorted = tmp_id(order);

t = table(tmp_id_sorted,Vmax_med_nonZero_sorted);
writetable(t,strcat(output_folder,'\','med_Vmax_sorted_noZscore.xlsx'));

Vmax_med_nonZero_sorted = normalize(Vmax_med_nonZero_sorted);

t = table(tmp_id_sorted,Vmax_med_nonZero_sorted);
writetable(t,strcat(output_folder,'\','med_Vmax_sorted.xlsx'));

% extract the colors
n = jet(length((tmp_id_sorted)));  % number of colors
n = flip(n);

%%  Graph heatmap
% Create a heatmap
figure(10);
h = heatmap(Vmax_med_nonZero_sorted,'Colormap', flip(n),...
    'ColorbarVisible','on','GridVisible','on','CellLabelColor','none');
h.XData = '.'; 
h.YData = tmp_id_sorted;
filename = strcat(output_folder,'\','heatmap');  % Specify the desired file name
saveas(gcf, strcat(filename,'.jpg'), 'jpg');
savefig(gcf, strcat(filename,'.fig'));

t = table(tmp_id,number_that_reach_max);
writetable(t,strcat(output_folder,'\','number_of_reactions_hit_vMax_1000.xlsx'));

%% Graph
for k = 1:length(gene_Vmax)
        i = order(k); % graph from the highest to lowest
        log_Vmax_sorted = gene_Vmax{1,i}.logVmax;
        log_Vmax_sorted = log_Vmax_sorted(~isinf(log_Vmax_sorted));

    x = 1:1:length(log_Vmax_sorted);
    cz = ones(length(x),1).*n(k,:);
    figure(1)
    s = scatter(x', log_Vmax_sorted,5,cz,'filled');
    xlabel('reactions','FontWeight','bold');
    ylabel('log_{10}Vmax (mmol g-dCW^{-1} hr^{-1})','FontWeight','bold');
    % Adjust the appearance of the chart as needed
    grid on;
    id_tmp{1,k} = gene_Vmax{1,i}.id;
    hold on    
end
yline(3,'--','Default max','Linewidth',1)
set(gca, 'FontSize',12);
l = legend(id_tmp,'Location','eastoutside','FontSize',8);
legend boxoff

% Save the scatter plot as a JPEG image
filename = strcat(output_folder,'\','scatter_plot');  % Specify the desired file name
saveas(gcf, strcat(filename,'.jpg'), 'jpg');
savefig(gcf, strcat(filename,'.fig'));

set(gcf,'Units','inches');
screenposition = get(gcf,'Position');
set(gcf,...
    'PaperPosition',[0 0 screenposition(3:4)],...
    'PaperSize',[screenposition(3:4)]);
% print -dpdf -painters scatter_all
print('-dpdf', '-painters', strcat(output_folder,'\','scatter_all.pdf'));

 Vmax_order_tissues = struct();
Vmax_order_tissues.Vmax_med_nonZero_sorted = Vmax_med_nonZero_sorted ;
Vmax_order_tissues.id_tmp = id_tmp';

for k = 1:length(gene_Vmax)
        i = order(k); % graph from the highest to lowest
        log_Vmax_sorted = gene_Vmax{1,i}.logVmax;
        log_Vmax_sorted = log_Vmax_sorted(~isinf(log_Vmax_sorted));

    x = 1:1:length(log_Vmax_sorted);
    % Sort the Vmax values in ascending order
    % Create a scatter plot with ascending order
    cz = ones(length(x),1).*[rand,rand,rand];
    
    figure(3)
    h = cdfplot(log_Vmax_sorted);
    h.LineWidth = 2;
    h.Color = n(k,:);
%     h.Color = [rand rand rand];
    xlabel('log_{10}Vmax (mmol g-dCW^{-1} hr^{-1})');
    ylabel('Fraction of reactions');
    grid on;
    hold on
    id_tmp{1,k} = gene_Vmax{1,i}.id;

end
xline(3,'--','Default max','Linewidth',1)
set(gca, 'FontSize',13);
legend(id_tmp,'Location','eastoutside','FontSize',8);
legend boxoff

% Save the scatter plot as a JPEG image
filename = strcat(output_folder,'\','cdf_plot');  % Specify the desired file name
saveas(gcf, strcat(filename,'.jpg'), 'jpg');
savefig(gcf, strcat(filename,'.fig'));

%% CDF - cal 
for k = 1:length(gene_Vmax)
  
    i = order(k); % graph from the highest to lowest
    log_Vmax_sorted = gene_Vmax{1,k}.logVmax;
    log_Vmax_sorted = log_Vmax_sorted(~isinf(log_Vmax_sorted));

%Empirical cumulative distribution function (cdf) plot
    [f,x,flo,fup] = ecdf(log_Vmax_sorted);
    
    [minDifference, index] = min(abs(f-0.5)); % find the idex that is closes to the 0.5
    log_fity_val = x(index);
    log_fity_val_store(k,1) = log_fity_val;
    log_fity_val_store_id{k,1} = gene_Vmax{1,k}.id;

end
[~, order_2] = sort(log_fity_val_store,'descend');

Vmax_order_tissues.log_fity_val_store = log_fity_val_store;
Vmax_order_tissues.log_fity_val_store_id = log_fity_val_store_id;

%% CDF graph
for k = 1:length(gene_Vmax)
  
    i = order_2(k); % graph from the highest to lowest
    log_Vmax_sorted = gene_Vmax{1,i}.logVmax;
    log_Vmax_sorted = log_Vmax_sorted(~isinf(log_Vmax_sorted));

    [f,x,flo,fup] = ecdf(log_Vmax_sorted);

    figure(4)
    h = stairs(x,f);   
    h.LineWidth = 2;
%     h.Color = [rand,rand,rand];
    h.Color = n(k,:);

    xlabel('Vmax (log_{10} mmol g-dCW^{-1} hr^{-1})');
    ylabel('Fraction of reactions');
    grid on;
    hold on;
    id_tmp{1,k} = gene_Vmax{1,i}.id; 
end   
    xline(3,'--','Linewidth',1)
    filename = strcat(output_folder,'\','cdf_plot_v2');  % Specify the desired file name
    set(gca, 'FontSize',13);
    legend(id_tmp,'Location','eastoutside','FontSize',8);
    legend boxoff

    % Define the y-position for the text label
yPosition = 0.1;  % Adjust this value as needed
% Define the rotation angle in degrees
rotationAngle = 90;  % Adjust this value as needed
% Add a rotated text label at the specified position
text(3.2, yPosition, 'Default max', 'HorizontalAlignment',...
        'center', 'BackgroundColor', 'none', 'Rotation', rotationAngle);
    
saveas(gcf, strcat(filename,'.jpg'), 'jpg');
savefig(gcf, strcat(filename,'.fig'));
end


