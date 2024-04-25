%% Load the model
load(strcat(pwd,'\models\','model_out_cbra2.mat')) % cbra2 is where the model id is updated
save_dirFlux = '3_fluxAnalz';
write = 0;
[data4comp_storage,subsysOut,FBAsolution_out] = run_fluxAnalz_v5_1(model_out_cbra2,write,save_dirFlux)

function [data4comp_storage,subsysOut,FBAsolution_out] = run_fluxAnalz_v5_1(model_out_cbra2,write,save_dirFlux)
pathway = pwd;
subfolder = [pathway '\' save_dirFlux];
if ~exist(subfolder, 'dir')
    mkdir(subfolder)
end

%% Set-up
model = model_out_cbra2;
%% Run Flux balance analysis
changeCobraSolver('gurobi','LP');

%% subsystem names
subsysName = {'Glycolysis / Gluconeogenesis', 'Pentose phosphate pathway',...
    'Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism',...
    'Nucleotide metabolism', 'Pyruvate metabolism', 'Alanine, aspartate and glutamate metabolism',...
    'Folate metabolism', 'Amino sugar and nucleotide sugar metabolism',...
    'Transport reactions'};

subsysName2 = {'Fatty acid oxidation', 'Beta oxidation of unsaturated fatty acids (n-7) (mitochondrial)',...
    'Beta oxidation of unsaturated fatty acids (n-9) (mitochondrial)',...
    'Beta oxidation of even-chain fatty acids (mitochondrial)',...
    'Beta oxidation of odd-chain fatty acids (mitochondrial)',...
    'Beta oxidation of branched-chain fatty acids (mitochondrial)'...
    'Beta oxidation of di-unsaturated fatty acids (n-6) (mitochondrial)',...
    'Beta oxidation of even-chain fatty acids (peroxisomal)',...
    'Beta oxidation of odd-chain fatty acids (peroxisomal)',...
    'Beta oxidation of phytanic acid (peroxisomal)',...
    'Beta oxidation of unsaturated fatty acids (n-7) (peroxisomal)',...
    'Beta oxidation of unsaturated fatty acids (n-9) (peroxisomal)',...
    };
subsysName3 = {'Fatty acid biosynthesis','Fatty acid biosynthesis (even-chain)',...
                'Cholesterol metabolism','Glycerophospholipid metabolism',...
                'Sphingolipid metabolism','Arginine and proline metabolism', 'CoA synthesis','Pantothenate and CoA biosynthesis'};
subsysName = [subsysName, subsysName2, subsysName3];
subsysSaveName = {'glycolysis','ppp','tca','oxphos','pyr','glu','folate','nucMe','transRxn'};
subsysSaveName2 = {'fao','bouf7','bouf9','boef','boof','bobcf','boduf',...
    'boefp','boofp','bopp','boufp7','boufp9'};
subsysSaveName3 = {'FAS','FAS_even','choles','glycerpho','Sphing','arg','CoASyn','Panto_CoASyn'};
subsysSaveName = [subsysSaveName, subsysSaveName2, subsysSaveName3];

nadpName = {'NADPHc','NADPHm'};
tag = 1; % to make flux label differently so that I can make outerjoin easily
%% Change the reaction bounds (v2)
for p = 1:length(model)
% GLUD - It is expected to be generating NADPH/NADH (Lewis, ARS, 2017)
    list = {'MAR03802','MAR03804'};
    model{p} = changeRxnBounds(model{p},list,0,'u'); 
% IDH - These reactions are reversible.
    list = {'MAR03958', 'MAR00710', 'MAR04111', 'MAR04112'};
    model{p} = changeRxnBounds(model{p},list,-1000,'l'); 
% Trxr-2 - Thioredoxin reductase uses NADPH to reduce oxidized thioredoxins
    list = {'MAR02358'};
    model{p} = changeRxnBounds(model{p},list,-1000,'l'); % 
    model{p} = changeRxnBounds(model{p},list,0,'u');
% Zw - It is considered irreversible reaction (the default directionality is opposite)
    list = {'MAR04306'};
    model{p} = changeRxnBounds(model{p},list,-1000,'l'); % 
    model{p} = changeRxnBounds(model{p},list,0,'u');
% Dhfr - difydrofolate makes THF. Set the reversible reaction to 0.
    list = {'MAR04332', 'MAR04333', 'MAR04335', 'MAR04654', 'MAR04655'};
    model{p} = changeRxnBounds(model{p},list,0,'l'); % 
% CG1236 / GRHPR - This reaction uses NAD(P)H to make glycolate. Set unidirectional.
    list = {'MAR07702'};
    model{p} = changeRxnBounds(model{p},list,0,'l'); % 
% FASN2 - NADPH is used for lipid synthesis in FASN2 reaction. The default directionality is opposite.  
    list = {'MAR02185'};
    model{p} = changeRxnBounds(model{p},list,0,'u'); % 
    model{p} = changeRxnBounds(model{p},list,-1000,'l'); %  
end

%% Obtain fluxes
% For single model, we do not need to compare.
if length(model_out_cbra2) == 1
    model_in = model;
    model_in = changeObjective(model_out_cbra2,'MAR00021'); %set the objective function
    FBAsolution  = optimizeCbModel(model_in,'max'); %FBA analysis
    FBAsolution_out.sol = FBAsolution;
    [P,C,vP,vC] = computeFluxSplits(model_in,{'MAM02555c[c]'},FBAsolution.x);
    T1 = [table(model_in.rxns(find(vP ~= 0)),'VariableNames',{'rxns'}) table(vP((find(vP ~= 0))))];
    filename = strcat(subfolder, '\NADPHc_vP.xlsx');
%     writetable(T1,filename,'Sheet',model_in.modelID)
    T1 = [table(model_in.rxns(find(vC ~= 0)),'VariableNames',{'rxns'}) table(vP((find(vC ~= 0))))];
    filename = strcat(subfolder, '\NADPHc_vC.xlsx');
%     writetable(T1,filename,'Sheet',model_in.modelID)
    
    sum_all = [sum(vP), sum(vC), sum(vP)-sum(vC)];
    filename = strcat(subfolder, '\NADPHc_sum.xlsx');
%     writematrix(sum_all,filename,'Sheet',model_in.modelID)

    
    [P,C,vp,vC] = computeFluxSplits(model_in,{'MAM02555m[m]'},FBAsolution.x);
    T1 = [table(model_in.rxns(find(vP ~= 0)),'VariableNames',{'rxns'}) table(vP((find(vP ~= 0))))];
    filename = strcat(subfolder, '\NADPHm_vP.xlsx');
%     writetable(T1,filename,'Sheet',model_in.modelID)
    T1 = [table(model_in.rxns(find(vC ~= 0)),'VariableNames',{'rxns'}) table(vC((find(vC ~= 0))))];
    filename = strcat(subfolder, '\NADPHm_vC.xlsx');
%     writetable(T1,filename,'Sheet',model_in.modelID)

    sum_all = [sum(vP), sum(vC), sum(vP)-sum(vC)];
    filename = strcat(subfolder, '\NADPHm_sum.xlsx');
%     writematrix(sum_all,filename,'Sheet',model_in.modelID)
    
    filename = strcat(subfolder, '\', model_in.modelID,'.xlsx');

        for j = 1:length(subsysName)
            [data,data4comp] = run_obtainSubsystemRxn(model_in,FBAsolution,subsysName{j},tag);
            writetable(data,filename,'Sheet',subsysSaveName{j})
    %         out{i,j} = data;
            data4comp_storage{1,j}=data4comp;            
            
            tag = tag+1;
        end  
        [data4comp_cyto,data4comp_mito] = run_obtainSubsystemRxnNADPH(model_in,FBAsolution,subfolder,tag);
        data4comp_storage{j+1} = data4comp_cyto; 
        data4comp_storage{j+2} = data4comp_mito; 
%           writetable(data4comp_cyto,filename,'Sheet',nadpName{1})
%          writetable(data4comp_mito,filename,'Sheet',nadpName{2})               
        flux_compare = data4comp_storage;
        subsysOut = [subsysSaveName,nadpName];
else
% Compare multiple models
    for i = 1:length(model)
        model{i} = addReaction(model{i},'ATP_maintenance',...
            'metaboliteList', ...
            {'MAM01371c[c]', 'MAM01371m[m]', 'MAM02040c[c]', 'MAM02040m[m]',... %ATPc,ATPm,H2Oc,H2Om
                'MAM01285c[c]', 'MAM01285m[m]','MAM02751c[c]','MAM02751m[m]',... %ADPc,ADPm,Pic,Pim
                'MAM02039c[c]', 'MAM02039m[m]'},...%Hc,Hm,NADPc,NADPm
            'stoichCoeffList', [-1; -1; -1; -1;...
                                        1; 1; 1; 1;...
                                            1; 1],...
            'reversible',false);

        model{i} = addReaction(model{i},'NAD_demand',...
            'metaboliteList', ...
            {'MAM02552c[c]', 'MAM02039c[c]', 'MAM02552m[m]', 'MAM02039m[m]',... %NADc,Hc,NADm,Hm,
                'MAM02553c[c]', 'MAM02553m[m]'},...%NADHc,NADHm
            'stoichCoeffList', [-1; -1; -1; -1;...
                                            1; 1],...
            'reversible',false);

        model{i}  = addReaction(model{i},'NADPH_demand',...
            'metaboliteList', ...
            {'MAM02555c[c]', 'MAM02555m[m]',... % NADPHc, NADPHm
                'MAM02039c[c]', 'MAM02039m[m]', 'MAM02554c[c]', 'MAM02554m[m]'},...%Hc,Hm,NADPc,NADPm
            'stoichCoeffList', [-1; -1;...
                                        1; 1; 1; 1;],...
            'reversible',false);

%% Set the solver
model_in = changeObjective(model{i},{'ATP_maintenance','NAD_demand','NADPH_demand'},1/3); %ATP synthase       
%     model_in = changeObjective(model{i},'MAR00021'); %set the objective function
    FBAsolution  = optimizeCbModel(model_in,'max'); %FBA analysis
    FBAsolution_out.sol{i} = FBAsolution;
        
    %% Check the flux balance - NADPH
    [P,C,vP,vC] = computeFluxSplits(model_in,{'MAM02555c[c]'},FBAsolution.x);
    T1 = [table(model_in.rxns(find(vP ~= 0)),'VariableNames',{'rxns'}) table(vP((find(vP ~= 0))))];
    filename = strcat(subfolder, '\NADPHc_vP.xlsx');
%     writetable(T1,filename,'Sheet',model_in.modelID);
%     T1.vP{i} = T1;
    T1 = [table(model_in.rxns(find(vC ~= 0)),'VariableNames',{'rxns'}) table(vC((find(vC ~= 0))))];
    filename = strcat(subfolder, '\NADPHc_vC.xlsx');
    
    [P,C,vp,vC] = computeFluxSplits(model_in,{'MAM02555m[m]'},FBAsolution.x);
    T2 = [table(model_in.rxns(find(vP ~= 0)),'VariableNames',{'rxns'}) table(vP((find(vP ~= 0))))];
    filename = strcat(subfolder, '\NADPHm_vP.xlsx');
%     T2.vP{i} = T2;
%     writetable(T1,filename,'Sheet',model_in.modelID);
    T2 = [table(model_in.rxns(find(vC ~= 0)),'VariableNames',{'rxns'}) table(vC((find(vC ~= 0))))];
    filename = strcat(subfolder, '\NADPHm_vC.xlsx');

    filename = strcat(subfolder, '\', model_in.modelID,'.xlsx');

        for j = 1:length(subsysName)
            [data,data4comp] = run_obtainSubsystemRxn(model_in,FBAsolution,subsysName{j},tag);
            if write ==1
                 writetable(data,filename,'Sheet',subsysSaveName{j})
            end
    %         out{i,j} = data;
            data4comp_storage{i,j}=data4comp;            
            tag = tag+1;
        end
        
        [data4comp_cyto,data4comp_mito] = run_obtainSubsystemRxnNADPH(model_in,FBAsolution,subfolder,tag);
        data4comp_storage{i,j+1} = data4comp_cyto; 
        data4comp_storage{i,j+2} = data4comp_mito; 
            if write ==1
                writetable(data4comp_cyto,filename,'Sheet',nadpName{1})
                writetable(data4comp_mito,filename,'Sheet',nadpName{2})
            end
    disp([model_in.modelID])
    end
    
    %% Remove the "MAR_xxxx" if existed. I found this. 
    % e.g., ) data4comp_storage{24, 7}   = 'MAR04442_r_ortho_1001toRest_inHum_u_v10uuu'
    for i = 1:size(data4comp_storage,1)
        for j = 1:size(data4comp_storage,2)
            idx = find(contains(data4comp_storage{i,j}.rxnIDs,'_'));
            data4comp_storage{i,j}.rxnIDs(idx) = extractBefore(data4comp_storage{i,j}.rxnIDs(idx),'_');
        end
    end
    
    save(strcat(subfolder,'\','data4comp_storage.mat'),'data4comp_storage')
    subsysOut = [subsysSaveName,nadpName];
end
    save(strcat(subfolder,'\','subsysOut.mat'),'subsysOut');
end
