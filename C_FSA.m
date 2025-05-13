dir = 'C_sampling';
pathway = pwd;
subfolder = [pathway '/' dir];
if ~exist(subfolder, 'dir')
    mkdir(subfolder)
end
changeCobraSolver('gurobi','all');
solverOK = changeCobraSolver('gurobi','LP');
%% Sampling
for k = 1:length(out_all_fvaBounded)
    subSysModel = out_all_fvaBounded{k,1};
    id = subSysModel.modelID;
    samples = gpSampler(subSysModel,10000);
    disp([num2str(k) '/' num2str(length(out_all_fvaBounded)) ', finished'])
    save_name = [pathway, '/', dir,'/', char(id),'.mat'];
    save(save_name,'samples')
end
