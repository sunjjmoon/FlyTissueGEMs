dir = 'C_sampling';
pathway = '/n/groups/perrimon/lab/sjm/fca_v2_muscle_diabetes/results/v3.2/gapdh_red_zwIncluded/gap_th0d3/OKSl/th0d5';
subfolder = [pathway '/' dir];
if ~exist(subfolder, 'dir')
    mkdir(subfolder)
end
changeCobraSolver('gurobi','all');
solverOK = changeCobraSolver('gurobi','LP');
%% Sampling
for k = 2:length(out_all_fvaBounded)
subSysModel = out_all_fvaBounded{k,1};
   id = subSysModel.modelID;

samples = gpSampler(subSysModel,10000);

   disp([num2str(k) '/' num2str(length(out_all_fvaBounded)) ', finished'])
save_name = [subfolder,'/', char(id),'.mat'];
save(save_name,'samples')
end
