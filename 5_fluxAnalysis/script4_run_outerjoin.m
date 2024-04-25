%% Goal:
%   To combine all the flux results using the run_outerjoin_all function
%   Input:
%       - input_model
%       - data4comp_storage
%       - subsysOut
%           -> These three will be automatically generated if I run the
%           fluxAnalz.
%       - save_dir
%           e.g.,) fluxAnalz_comb
%% Load the model
load(strcat(pwd,'\models\','model_out_cbra2.mat')) % cbra2 is where the model id is updated
load(strcat(pwd,'\3_fluxAnalz\','data4comp_storage.mat')) 
save_dirFlux = '4_outerjoin';
write = 1;

data_store3 = run_outerjoin_f(model_out_cbra2,subsysOut,data4comp_storage,save_dirFlux)

%% Functions
function data_store3 = run_outerjoin_f(model_out_cbra2,subsysOut,data4comp_storage,save_dirFlux)

%% Set the save folder
pathway = pwd;
subfolder = [pathway '\' save_dirFlux];
if ~exist(subfolder, 'dir')
    mkdir(subfolder)
end

modelID_names =[];
for i = 1:length(model_out_cbra2)
   modelID_names{i} = model_out_cbra2{i,1}.modelID; 
end
modelID_names= modelID_names';

%% Combine all the data
data_joined_tmp=[];
for j = 1:length(subsysOut)
    data_joined_tmp = data4comp_storage{1,j};
    for i = 1:length(modelID_names)-1
            data_joined_tmp=outerjoin(data_joined_tmp, data4comp_storage{i+1,j}); %always combine with control; control on the left
            %convert variable name so that I can do outerjoin later 
            data_joined_tmp.Properties.VariableNames(1) = cellstr('rxnIDs');  
%             data_joined_tmp(:,end-1) =[]; % remove the the exp colm
    end
    data_store{1,j} = data_joined_tmp;
end

filename = strcat(subfolder, '\flux_comp_all_manual.xlsx');

%% Fill the empty name with the names
col_name_tmp = 1:2:2*length(modelID_names)-1; %idx of cols that represent rxnID of ctr + all exps
data_store2 = data_store;
for k = 1:length(subsysOut) % k is the subsystem
    k
    subsysOut(k)
    for i = 1:size(data_store2{1,k},1)
        check = ismissing(data_store2{1,k}(i,1));
        if check == 1 % if empty array
           for j = 1:length(col_name_tmp) % now go through the columns (checking whether exp cols have the rxnID
               tf = ismissing(data_store2{1,k}(i,col_name_tmp(j))); 
               if tf ~=1 % if the array is not empty
                   data_store2{1,k}(i,1) = data_store2{1,k}(i,col_name_tmp(j)); % update the first empty col name with this.
                   for jj = 2:length(col_name_tmp) % Now, I need to check other exp col contains this reaction. If so, fill it with the values. Go from 2, because we just fill control (jj=1)
                       idx = find(contains(table2cell(data_store2{1,k}(:,col_name_tmp(jj))),table2cell((data_store2{1,k}(i,col_name_tmp(j))))));
                       check_empty = isempty(idx);
                       if numel(idx) >1  
                           idx = idx(1);
                       end
                % end 
                       if check_empty ~=1 % if index is not empty, meaning there is already this reaction. Then, fill the array with pre-determined values
                          data_store2{1,k}(i,col_name_tmp(jj)+1) = data_store2{1,k}(idx,col_name_tmp(jj)+1); % +1 is the where the value is position. After the rxnID column.
                       end
                   end
               end
           end
        end  
    end
    col_name_tmp_rmov = 3:2:2*length(modelID_names)-1; %idx of cols that represent rxnID of ctr + all exps
    data_store2{1,k}(:,col_name_tmp_rmov) = []; % update the first empty col name with this.
    data_store2{1,k}.Properties.VariableNames = [{'rxnIDs'}; modelID_names];
end
filename = strcat(subfolder, '\flux_comp_all_f_man.xlsx');
disp(filename)

%% Use unique function to report only the the unique
%% Report only unique
for p = 1:length(subsysOut)
    [C,ia,ib] = unique(data_store2{1,p}(:,1));
%     [uniqueA u y] = unique(ib,'first');
%     indexToDupes = find(not(ismember(1:numel(ib),u)))
    data_store3{1,p} = data_store2{1,p}(ia,:);
end
disp('removed redundant')
filename = strcat(subfolder, '\flux_comp_all_ff_man.xlsx');
    for i = 1:length(subsysOut)
            writetable(data_store2{1,i},filename,'Sheet',subsysOut{i});            
    end

    save(strcat(subfolder,'\','data_store3.mat'),'data_store3')
end