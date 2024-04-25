%% Set the save and load folder
pathway = pwd;
loadfolder = [pathway '\' 'B_fva_update_v2'];

load(strcat(loadfolder,'\','out_all_fvaBounded.mat'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter setting - Type interest and save_dir
interest = 'Reactions for subsystem glycolysis_gluconeogenesis';
save_dirFlux = 'glycolysis_gluconeogenesis';

% save folder define
subfolder = [pathway '\' 'C_graph_fva' '\' save_dirFlux];
if ~exist(subfolder, 'dir')
    mkdir(subfolder)
end

modelID = graph_box_max_min(out_all_fvaBounded,subfolder,interest);

function modelID = graph_box_max_min(out_all_fvaBounded,subfolder,interest)
close all
%% Load the files
pathway = pwd;
loadfolder = [pathway '\rxn_path_files\'];
filename = strcat(loadfolder, interest, '.tsv');
opts = detectImportOptions(filename, 'FileType', 'text');
opts.Delimiter = '\t';
opts.Encoding = 'UTF-8';
tmp1 = readtable(filename, opts);

reactions = string(tmp1.ReactionID);
reactionsID = string(tmp1.Genes);
trimmedStr = strtrim(reactionsID);
eqn = string(tmp1.Equation);

%% Combine the reaction with the gene
rxn_gene = cellfun(@(r, g) strcat(r, ' (', g, ')'), reactions, trimmedStr, 'UniformOutput', false);

for j = 1:length(reactions)
list = reactions(j);

%% Extract the reaction of interest lb and ub
    lb = zeros(length(out_all_fvaBounded),1);
    ub = zeros(length(out_all_fvaBounded),1);
    for i = 1:length(out_all_fvaBounded)
        modelID{i,1} = out_all_fvaBounded{i, 1}.modelID;
        idxTmp = find(ismember(out_all_fvaBounded{i,1}.rxns,list));

       if isempty(idxTmp) ~= 1 % If the reaction exists
            lb(i,1) = out_all_fvaBounded{i,1}.lb(idxTmp);
            ub(i,1) = out_all_fvaBounded{i,1}.ub(idxTmp);    
           
       else
           lb(i,1) = 0;
           ub(i,1) = 0;
       end 
    end
    
    %scale
    lb = (lb)./1000; ub = (ub)./1000;

    mean_val = (abs(ub)+abs(lb))./2+lb;
%     mean_val = ub-lb;
    range_mean = mean(mean_val);
    
%% Graph
g = categorical(modelID);
g = reordercats(g);
% Create a bar graph with error bars
% figure(1);
% b = barh(g,lb,'FaceColor','none','EdgeColor','k','LineWidth',2);
% hold on;
% b2 = barh(g,ub,'FaceColor','none','EdgeColor','k','LineWidth',2);
% xlim([-1.1 1.1]);
idx_oo_store = [];
lower_end_store = [];
for oo = 1:length(ub)
    if (ub(oo) > 0) && (lb(oo) > 0)
        upper(oo,1) = ub(oo)-lb(oo); % since the stacked values will be greater than 1.
        idx_oo_store(oo,1) = 1; %  %both are positive
    elseif (ub(oo) < 0) && (lb(oo) < 0)
        upper(oo,1) = ub(oo)-lb(oo); % since the stacked values will be greater than 1.        
        idx_oo_store(oo,1) = -1; % 
%         lower(oo,1) = lb(oo)-ub(oo);
        lower(oo,1) = upper(oo,1).*-1;

    else
        upper(oo,1) = ub(oo);
        lower(oo,1) = lb(oo);
%         upper(oo,1) = ub(oo)-lb(oo); % since the stacked values will be greater than 1.        
        idx_oo_store(oo,1) = 0; % no problem 
    end
end

% Contains both negative and positive
idx_oo_store_v2 = [];
dir_store = zeros(length(ub),1);
for oo = 1:length(ub)
    multiple = ub(oo).*lb(oo);
    dir = ub(oo) >= 0;
    if multiple < 0
        idx_oo_store_v2(oo,1) = -1; % range spans to both directions
    else 
        idx_oo_store_v2(oo,1) = 0; % no problem 
        if dir == 1 % both are positive direction 
            dir_store(oo,1) = 1;
        else
            dir_store(oo,1) = -1;
        end
    end
end

    
figure(1)
if sum(ismember(idx_oo_store_v2,0)) == length(lb) % if lower limit is positive and upper limit is positive
    if dir_store(1) == -1  % meaning that both upper and lower is negative direction
         h = barh(g,[ub lower],'stacked','BaseValue',max(ub),'FaceColor',[0.5 0.5 0.5]);
        for u = 1:numel(h)
            h(u).FaceColor = [0.5 0.5 0.5];
            if u == 1
                h(u).Visible = 'off';
            end
        end  
    else 
        h = barh(g,[lb upper],'stacked','BaseValue',min(lb),'FaceColor',[0.5 0.5 0.5]);

        for u = 1:numel(h)
            h(u).FaceColor = [0.5 0.5 0.5];
            if u == 1
                h(u).Visible = 'off';
            end
        end
    end
    
else % If the flux range spans in both direction
    if (min(lb) < 0 && min(ub) >= 0)  
    h = barh(g,[lb upper],'stacked','BaseValue',min(lb),'FaceColor',[0.5 0.5 0.5]);
        for u = 1:numel(h)
            h(u).FaceColor = [0.5 0.5 0.5];
            if u == 1
                h(u).Visible = 'off';
            end
        end
            hold on
             h3 = barh(g,[lb],'stacked','BaseValue',0);
             h3.FaceColor = 'flat';
%              h3.EdgeColor =  [1 1 1];
             h3.CData = [1 1 1].*ones(length(lb),1);
%              h3.EdgeColor = [1 1 1].*ones(length(lb),1);

             idx_min = find(lb < 0);
             h3.CData(idx_min,:) = [0.5 0.5 0.5].*ones(length(idx_min),1);
    
    elseif (min(lb) < 0 && min(ub) < 0) 
    h = barh(g,[ub lower],'stacked','BaseValue',max(ub),'FaceColor',[0.5 0.5 0.5]);
        for u = 1:numel(h)
            h(u).FaceColor = [0.5 0.5 0.5];
            if u == 1
                h(u).Visible = 'off';
            end
        end
            hold on
             h3 = barh(g,[ub],'stacked','BaseValue',0);
             h3.FaceColor = 'flat';
             h3.CData = [1 1 1].*ones(length(lb),1);
             hold on 
             idx_min = find(ub > 0);
             h3.CData(idx_min,:) = [0.5 0.5 0.5].*ones(length(idx_min),1);
    else
        h = barh(g,[lb],'stacked','BaseValue',min(lb),'FaceColor',[0.5 0.5 0.5]);
        hold on
        h2 = barh(g,[ub],'stacked','BaseValue',0,'FaceColor',[0.5 0.5 0.5]);
    
    end
end
xline(0,'--','LineWidth',1)
% Set the limit
if max(ub)>0
    xlim_range = [min(lb)*1.03,max(ub)*1.03];
    if sum(xlim_range) ~= 0 
        set(gca,'XLim',[min(lb)*1.03,max(ub)*1.03]);
    end
else
    xlim_range = [min(lb)*1.03,max(ub)/1.03];
    if sum(xlim_range) ~= 0 
        set(gca,'XLim',[min(lb)*1.03,max(ub)/1.103]);
    end
end

xlabel('Flux range (max - min, normalized)')
xline(range_mean,'--','LineWidth',0.5);

%% PDF file does not accept the -->. change.
title_tmp = eqn(j);
% Replace "⇒" with "->"
title_tmp = strrep(title_tmp, '⇒', '->');
% Replace "⇔" with "<->"
title_tmp = strrep(title_tmp, '⇔', '<->');

title(title_tmp,'FontSize',6)
grid on
grid minor

hold off


%     filename = strcat(subfolder, '\', reactions(j), '.jpeg');
%     print(gcf,filename, '-djpeg', '-r300')
    exportgraphics(gcf,strcat(subfolder, '\', reactions(j), '.pdf'))
tmp = [lb ub];
end
end
