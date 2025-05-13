%% Set the save and load folder
pathway = pwd;

loadfolder = [pathway '\' '2_FBA_FVA_rxn_bound_updated'];
load(strcat(pathway, '\2_FBA_FVA_rxn_bound_updated\', 'out_all_fvaBounded.mat'));

% interest file
interest = 'Reactions for subsystem glycolysis_gluconeogenesis';
save_dirFlux = 'glycolysis_gluconeogenesis';

% Define save folder
subfolder = [pathway '\' '3_sensitivity'];
if ~exist(subfolder, 'dir')
    mkdir(subfolder)
end
sens_int = {'MAR04363'};

graph = 0;
modelID = graph_box_max_min(out_all_fvaBounded, subfolder, interest, graph,sens_int);

function modelID = graph_box_max_min(out_all_fvaBounded, subfolder, interest, graph,sens_int)
close all
pathway = pwd;
[parentFolder, ~, ~] = fileparts(pathway);
loadfolder = [parentFolder '\rxn_path_files\'];

filename = strcat(loadfolder, interest, '.tsv');
opts = detectImportOptions(filename, 'FileType', 'text');
opts.Delimiter = '\t';
opts.Encoding = 'UTF-8';
tmp1 = readtable(filename, opts);

reactions = string(tmp1.ReactionID);
reactionsID = string(tmp1.Genes);
eqn = string(tmp1.Equation);
trimmedStr = strtrim(reactionsID);
rxn_gene = cellfun(@(r, g) strcat(r, ' (', g, ')'), reactions, trimmedStr, 'UniformOutput', false);

for j = 1:length(reactions)
    list = reactions(j);
    lb = zeros(length(out_all_fvaBounded), 1);
    ub = zeros(length(out_all_fvaBounded), 1);

    for i = 1:length(out_all_fvaBounded)
        modelID{i, 1} = out_all_fvaBounded{i, 1}.modelID;
        idxTmp = find(ismember(out_all_fvaBounded{i, 1}.rxns, list));
        if ~isempty(idxTmp)
            lb(i, 1) = out_all_fvaBounded{i, 1}.lb(idxTmp);
            ub(i, 1) = out_all_fvaBounded{i, 1}.ub(idxTmp);
        else
            lb(i, 1) = 0;
            ub(i, 1) = 0;
        end
    end

    lb = lb ./ 1000;
    ub = ub ./ 1000;
    mean_val = (abs(ub) + abs(lb)) ./ 2 + lb;
    range_mean = mean(mean_val);

    if graph == 1
        g = categorical(modelID);
        g = reordercats(g);

        idx_oo_store = [];
        upper = zeros(length(ub), 1);
        lower = zeros(length(ub), 1);

        for oo = 1:length(ub)
            if ub(oo) > 0 && lb(oo) > 0
                upper(oo) = ub(oo) - lb(oo);
                idx_oo_store(oo) = 1;
            elseif ub(oo) < 0 && lb(oo) < 0
                upper(oo) = ub(oo) - lb(oo);
                lower(oo) = -upper(oo);
                idx_oo_store(oo) = -1;
            else
                upper(oo) = ub(oo);
                lower(oo) = lb(oo);
                idx_oo_store(oo) = 0;
            end
        end

        idx_oo_store_v2 = zeros(length(ub), 1);
        dir_store = zeros(length(ub), 1);

        for oo = 1:length(ub)
            if ub(oo) * lb(oo) < 0
                idx_oo_store_v2(oo) = -1;
            else
                idx_oo_store_v2(oo) = 0;
                dir_store(oo) = ub(oo) >= 0;
                if dir_store(oo) == 0
                    dir_store(oo) = -1;
                end
            end
        end

        figure(1)
        if all(idx_oo_store_v2 == 0)
            if dir_store(1) == -1
                h = barh(g, [ub lower], 'stacked', 'BaseValue', max(ub));
            else
                h = barh(g, [lb upper], 'stacked', 'BaseValue', min(lb));
            end

            for u = 1:numel(h)
                h(u).FaceColor = [0.5 0.5 0.5];
                if u == 1
                    h(u).Visible = 'off';
                end
            end
        else
            if min(lb) < 0 && min(ub) >= 0
                h = barh(g, [lb upper], 'stacked', 'BaseValue', min(lb));
            elseif min(lb) < 0 && min(ub) < 0
                h = barh(g, [ub lower], 'stacked', 'BaseValue', max(ub));
            else
                h = barh(g, [lb], 'stacked', 'BaseValue', min(lb));
                h2 = barh(g, [ub], 'stacked', 'BaseValue', 0);
            end

            for u = 1:numel(h)
                h(u).FaceColor = [0.5 0.5 0.5];
                if u == 1
                    h(u).Visible = 'off';
                end
            end

            hold on
            h3 = barh(g, [lb], 'stacked', 'BaseValue', 0);
            h3.FaceColor = 'flat';
            h3.CData = ones(length(lb), 1);
            idx_min = find(lb < 0);
            h3.CData(idx_min, :) = 0.5 * ones(length(idx_min), 1);
        end

        xline(0, '--', 'LineWidth', 1)

        if max(ub) > 0
            set(gca, 'XLim', [min(lb)*1.03, max(ub)*1.03]);
%         else
%             set(gca, 'XLim', [min(lb)*1.03, max(ub)/1.103]);
        end

        xlabel('Flux range (max - min, normalized)')
        xline(range_mean, '--', 'LineWidth', 0.5);
        title(eqn(j), 'FontSize', 6)
        grid on; grid minor; hold off

        filename = strcat(subfolder, '\', reactions(j), '.png');
        print(gcf, filename, '-djpeg', '-r300');
    end
end

%% Summary extraction
summaryData = cell(length(reactions), 5);
[lb_nsd, ub_nsd, mean_nsd] = deal(zeros(length(reactions),1));
[lb_hsd, ub_hsd, mean_hsd] = deal(zeros(length(reactions),1));

for j = 1:length(reactions)
    for i = 1:length(out_all_fvaBounded)
        idxTmp = find(ismember(out_all_fvaBounded{i,1}.rxns, reactions(j)));
        if i == 1
            if ~isempty(idxTmp)
                lb_nsd(j) = out_all_fvaBounded{i,1}.lb(idxTmp);
                ub_nsd(j) = out_all_fvaBounded{i,1}.ub(idxTmp);
                mean_nsd(j) = (lb_nsd(j) + ub_nsd(j)) / 2;
            end
        else
            if ~isempty(idxTmp)
                lb_hsd(j) = out_all_fvaBounded{i,1}.lb(idxTmp);
                ub_hsd(j) = out_all_fvaBounded{i,1}.ub(idxTmp);
                mean_hsd(j) = (lb_hsd(j) + ub_hsd(j)) / 2;
            end
        end
    end
end

rxn_int_idx = find(ismember(reactions, sens_int));
rxn_int_mean_nsd = mean_nsd(rxn_int_idx);
rxn_int_mean_hsd = mean_hsd(rxn_int_idx);
del_rxn_int = rxn_int_mean_hsd - rxn_int_mean_nsd;
abs_del_rxn_int_norm = (del_rxn_int ./ rxn_int_mean_nsd);

del_sum_abs_change_norm = abs(nanmean((mean_hsd - mean_nsd) ./ mean_nsd));
mean_relChange = mean_hsd ./ mean_nsd;
mean_relChange_log2 = log2(mean_relChange);

sum_nsd_mean = sum(abs(mean_nsd));
sum_hsd_mean = sum(abs(mean_hsd));
mean_relChange_sum = sum_hsd_mean / sum_nsd_mean;
mean_relChange_log2_sum = log2(mean_relChange_sum);
sensitivity = del_sum_abs_change_norm / abs_del_rxn_int_norm;

summaryTable = table(reactions, rxn_gene, lb_nsd, ub_nsd, lb_hsd, ub_hsd, mean_nsd, mean_hsd, mean_relChange, mean_relChange_log2);
summary_mean = table(sum_nsd_mean, sum_hsd_mean, mean_relChange_sum, mean_relChange_log2_sum, sensitivity, abs_del_rxn_int_norm, del_sum_abs_change_norm);


writetable(summaryTable, strcat(subfolder, '\summary_data_v3_avg.xlsx'));
writetable(summary_mean, strcat(subfolder, '\summary_data_v3_avg.xlsx'), 'Sheet', 'summary');

%% Summary extraction - v2
rxns_int = {'MAR04394','MAR04381','MAR04379','MAR04375','MAR04373',...
    'MAR04368','MAR04365','MAR04363','MAR04358','MAR04388','MAR04391'}';

summaryData_v2 = cell(length(rxns_int), 5);
[lb_nsd_v2, ub_nsd_v2, mean_nsd_v2] = deal(zeros(length(rxns_int),1));
[lb_hsd_v2, ub_hsd_v2, mean_hsd_v2] = deal(zeros(length(rxns_int),1));

for j = 1:length(rxns_int)
    for i = 1:length(out_all_fvaBounded)
        idxTmp = find(ismember(out_all_fvaBounded{i,1}.rxns, rxns_int(j)));
        if i == 1
            if ~isempty(idxTmp)
                lb_nsd_v2(j) = out_all_fvaBounded{i,1}.lb(idxTmp);
                ub_nsd_v2(j) = out_all_fvaBounded{i,1}.ub(idxTmp);
                mean_nsd_v2(j) = (lb_nsd_v2(j) + ub_nsd_v2(j)) / 2;
            end
        else
            if ~isempty(idxTmp)
                lb_hsd_v2(j) = out_all_fvaBounded{i,1}.lb(idxTmp);
                ub_hsd_v2(j) = out_all_fvaBounded{i,1}.ub(idxTmp);
                mean_hsd_v2(j) = (lb_hsd_v2(j) + ub_hsd_v2(j)) / 2;
            end
        end
    end
end

rxn_int_idx = find(ismember(reactions, sens_int));
rxn_int_mean_nsd = mean_nsd(rxn_int_idx);
rxn_int_mean_hsd = mean_hsd(rxn_int_idx);
del_rxn_int = rxn_int_mean_hsd - rxn_int_mean_nsd;
abs_del_rxn_int_norm = (del_rxn_int ./ rxn_int_mean_nsd);

del_sum_abs_change_norm = abs(nanmean((mean_hsd_v2 - mean_nsd_v2) ./ mean_nsd_v2));
mean_relChange = mean_hsd_v2 ./ mean_nsd_v2;
mean_relChange_log2 = log2(mean_relChange);

sum_nsd_mean = sum(abs(mean_nsd_v2));
sum_hsd_mean = sum(abs(mean_hsd_v2));
mean_relChange_sum = sum_hsd_mean / sum_nsd_mean;
mean_relChange_log2_sum = log2(mean_relChange_sum);
sensitivity = del_sum_abs_change_norm / abs_del_rxn_int_norm;

% summaryTable = table(reactions, rxn_gene, lb_nsd, ub_nsd_v2, lb_hsd_v2, ub_hsd_v2, mean_nsd_v2, mean_hsd_v2, mean_relChange, mean_relChange_log2);
summary_mean = table(sum_nsd_mean, sum_hsd_mean, mean_relChange_sum, mean_relChange_log2_sum, sensitivity, abs_del_rxn_int_norm, del_sum_abs_change_norm);


% writetable(summaryTable, strcat(subfolder, '\summary_data_v3_v2_avg.xlsx'));
writetable(summary_mean, strcat(subfolder, '\summary_data_v3_v2_avg.xlsx'), 'Sheet', 'summary');

end
