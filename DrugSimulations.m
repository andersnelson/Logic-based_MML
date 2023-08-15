% Written by Anirudha Chandrabhatla 
% Updated for Sens Analsysis by Anders Nelson
% Last Updated 11/9/2021
% Version 1.0

% %% If any inputs have changed, run the following section: 
% 
%   warning off;
% % Species information from fib617_references.xlsx, 'species' tab 
% nodeGenesTable = readtable('fib617_references.xlsx', 'Sheet', 'species');
% 
% % Reaction information from fib417_references.xlsx, 'reactions' tab 
% networkReactions = readtable('fib617_references.xlsx', 'Sheet', 'reactions');
% 
% % Drug information from AllPharmActive.xlsx
% approvedTargetsTable = readtable('AllPharmActive.xlsx');
% 
% % Drug ID to drug name matching from drugIDsAndNamesFull516.xlsx
% drugIDToNameConversion = readtable('drugIDsAndNamesFull516.xlsx');
% 
% % Manually curated competitive/non-competitive classifications for all drugs that target a node in the network
% compNonCompClassification = readtable('All Drugs Comp NonComp Classification.xlsx', 'Sheet', 'All Drugs Condensed');
% 
% load('finalDrugOutputNetworkTargets.mat'); % Result from the webscraper
% drugInformation = finalDrugOutputNetworkTargets;
% 
% geneNames = nodeGenesTable.('geneName');  % Gene names from the model
% geneIDs = nodeGenesTable.('ID');          % Gene IDs from the model
% 
% % Drug IDs from the database of approved targets
% drugIDs_fromTargetList = approvedTargetsTable.('DrugIDs'); 
% % Gene names of the drug targets. Retrieved from the database of approved targets.
% drugGeneTargets = approvedTargetsTable.('GeneName'); 
% 
% % Drug names from the database of all drugs
% drugNames = drugIDToNameConversion.('Name');
% % Drug IDs from the database of all drugs 
% drugIDs_fromTargetList_forNameConversion = drugIDToNameConversion.('DrugBankID');
% 
% [drugsToSimulate, formattedReactions] = formatInputs(geneNames, geneIDs, drugGeneTargets, drugIDs_fromTargetList, drugIDs_fromTargetList_forNameConversion, drugNames, networkReactions, compNonCompClassification, drugInformation);

%% Inputs for simulations

drugsToSimulate = load('Anders_drugsToSimulate_srcRas.mat'); %edited drugs list
drugsToSimulate = drugsToSimulate.drugsToSimulate; % Extract from struct

formattedReactions = load('formattedReactions.mat');
formattedReactions = formattedReactions.formattedReactions; % Extract from struct

clusteredRowLabels = load('Anders_clusteredRowLabels.mat');
clusteredRowLabels = clusteredRowLabels.clusteredRowLabels; % Extract from struct

% Set drug dose or doses (as a vector)
COV=0.08485638221813657;
%allocate random basal input values and dose response + increased cytokine
%values
doseResponse = normrnd(0.85,0.85*COV,1);
randn('seed',0);
ensemble=[];
%allocate random basal input values
rand_inputs=normrnd(0.25,0.25*COV,[1,9])
while any(rand_inputs < 0)
    rand_inputs=normrnd(0.1,0.25*COV,[1,9])
end
rand_mech=normrnd(0.725,0.1427)
while rand_mech < 0 | 1 < rand_mech
    rand_mech=normrnd(0.725,0.1427)
end    
%% Part 2: Generate graph of results

inputNodeW = {5,2,[5,2]};
inputLabels = {'Control','IL1','TGFB','TGFB+IL1'};
deltaIn = 0.6;

% Preturbed simulations
sens_col = zeros(height(drugsToSimulate), length(inputLabels));
sens_edafn = zeros(height(drugsToSimulate), length(inputLabels));
sens_prolif = zeros(height(drugsToSimulate), length(inputLabels)); 

sens_mmp1 = zeros(height(drugsToSimulate), length(inputLabels)); 
sens_mmp2 = zeros(height(drugsToSimulate), length(inputLabels)); 
sens_mmp9 = zeros(height(drugsToSimulate), length(inputLabels));
sens_all=[];



drugTargetResponse = [];
drugTargetResponseID = zeros(1, height(drugsToSimulate));

drugsWithNoChangeInCollagen = [];
dataVector = [];
entrestoAgonistResponse = [];
ctr=1;
sensMatrix=[];
for k = 1:length(inputLabels)
    for q = 1:5 %set length for ensemble sim per condition
    % Parameters and initial conditions

    %%% START ENSEMBLE CODE 
    COV=0.08485638221813657;
%allocate random basal input values and dose response + increased cytokine
%values
doseResponse = normrnd(0.85,0.85*COV,1);
randn('seed',0);
ensemble=[];
%allocate random basal input values
rand_inputs=normrnd(0.25,0.25*COV,[1,9])
while any(rand_inputs < 0)
    rand_inputs=normrnd(0.1,0.25*COV,[1,9])
end
rand_mech=normrnd(0.85,0.85*COV,1)
while rand_mech < 0 | 1 < rand_mech
    rand_mech=normrnd(0.85,0.85*COV,1)
end    
%%% END ENSEMBLE CODE
    [params,y0] = fib617_params;
    tspan = [0 500]; 
    options = [];

    % Pull parameters out in order to alter them in the loop without negating the previous alteration
    [rpar,tau,ymax,speciesNames,KI]=params{:}; 
    w = rpar(1,:);
    wNew = rpar(1,:);
    wNew(1:11) = 0.25;
    if k > 1 % First run is a control simulation. 
        wNew(cell2mat(inputNodeW(k-1))) = deltaIn;
    end
    wNew(3) = 0.85; % Setting mechanical input
    wNew(23) = 0.6;
    wNew(24) = 0.6;
    wNew(25) = 0; % Altering the TGFBmRNA --> TGFB reaction weight
    n = rpar(2,:);
    EC50 = rpar(3,:);
    dose = rpar(4,:);
    drugType = rpar(5,:);
    drugBinding = rpar(5,:); 
    drugAgonism = rpar(6,:); 

    rpar = [wNew;n;EC50;dose;drugBinding;drugAgonism];


    % Default simulation creating a vector showing steady state activation of all species
    params = {rpar,tau,ymax,speciesNames,KI};
    [t,y] = ode15s(@fib617_modelODE_drug_static,tspan,y0,options,params);
    ycol0(k) = y(end,101);
    sens_all=y(end,:);
    sensMatrix=[sensMatrix;sens_all]; %add a control simulation for the cytokine context alone

    for u = 1:107
        rpar = [wNew;n;EC50;dose;drugBinding;drugAgonism];
        params = {rpar,tau,ymax,speciesNames,KI};
        ymaxNew=ymax;
        ymaxNew(u)=0; %knockout node
        paramsNew = {rpar,tau,ymaxNew,speciesNames,KI}; % Only need if doing a knockout
        [t2,y2] = ode15s(@fib617_modelODE_drug_static,tspan,y0,options,paramsNew);
        sens_all=y2(end,:);
        %ctr=ctr+1;

        sensMatrix=[sensMatrix;sens_all];
        end
    
    
    dataVector = [];
    % Iterate through 'drugsToSimulate' and simulate the addition of each drug into the network.
    for i = 1:height(drugsToSimulate)
        if i == 32
            disp('here')
        end

        disp(drugsToSimulate{i, 1}); % Prints drug name to the command window  

        for j = 1:length(doseResponse)

            doseNew = dose;
            drugBindingNew = drugBinding;
            drugAgonismNew = drugAgonism;

            % The drug is an agonist
            if strcmp(drugsToSimulate.IsAgonist{i}, 'Yes') == 1
                if isempty(find(drugsToSimulate.AgonistTarget{i} == ';')) % Drug has one agonist target
                    locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, drugsToSimulate.AgonistTargetGeneID{i});
                    if strcmp('Competitive', drugsToSimulate.DrugAction{i}) % Competitive, agonist
                        drugBindingNew(locationOfReactions) = 1;
                        drugAgonismNew(locationOfReactions) = 1;
                        doseNew(locationOfReactions) = -1*doseResponse(j);
                    else % Non-Competitive, agonist
                        drugBindingNew(locationOfReactions) = -1;
                        drugAgonismNew(locationOfReactions) = 1;
                        doseNew(locationOfReactions) = doseResponse(j);
                    end
                else % Drug has multiple targets
                    geneIDsOfTargets = strsplit(drugsToSimulate.AgonistTargetGeneID{i}, ';');
                    for m = 1:length(geneIDsOfTargets)
                        locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, geneIDsOfTargets{m});
                        if strcmp('Competitive', drugsToSimulate.DrugAction{i})
                            drugBindingNew(locationOfReactions) = 1;
                            drugAgonismNew(locationOfReactions) = 1;
                            doseNew(locationOfReactions) = -1*doseResponse(j);
                        else 
                            drugBindingNew(locationOfReactions) = -1;
                            drugAgonismNew(locationOfReactions) = 1;
                            doseNew(locationOfReactions) = doseResponse(j);
                        end
                    end
                end
            end

            % The drug is an antagonist
            if strcmp(drugsToSimulate.IsAntagonist{i}, 'Yes') == 1
                if isempty(find(drugsToSimulate.AntagonistTarget{i} == ';')) % Drug has one antagonist target
                    locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, drugsToSimulate.AntagonistTargetGeneID{i});
                    if strcmp('Competitive', drugsToSimulate.DrugAction{i}) % Competitive, antagonist
                        drugBindingNew(locationOfReactions) = 1;
                        drugAgonismNew(locationOfReactions) = -1;
                        doseNew(locationOfReactions) = doseResponse(j);
                    else
                        drugBindingNew(locationOfReactions) = -1;
                        drugAgonismNew(locationOfReactions) = -1;
                        doseNew(locationOfReactions) = doseResponse(j);
                    end
                else
                    geneIDsOfTargets = strsplit(drugsToSimulate.AntagonistTargetGeneID{i}, ';');
                    for p = 1:length(geneIDsOfTargets)
                        locationOfReactions = strcmp(formattedReactions.ReactionOutputNode, geneIDsOfTargets{p});
                        if strcmp('Competitive', drugsToSimulate.DrugAction{i})
                            drugBindingNew(locationOfReactions) = 1;
                            drugAgonismNew(locationOfReactions) = -1;
                            doseNew(locationOfReactions) = doseResponse(j);
                        else 
                            drugBindingNew(locationOfReactions) = -1;
                            drugAgonismNew(locationOfReactions) = -1;
                            doseNew(locationOfReactions) = doseResponse(j);
                        end
                    end
                end
            end
            disp(num2str(wNew(3)));
            rparNew = [wNew;n;EC50;doseNew;drugBindingNew;drugAgonismNew];
%             paramsNew = {rparNew,tau,ymaxNew,speciesNames,KI}; % Only need if doing a knockout
            paramsNew = {rparNew,tau,ymax,speciesNames,KI};
            [t2,y2] = ode15s(@fib617_modelODE_drug_static,tspan,y0,options,paramsNew);
            sens_all=y2(end,:);
            sensMatrix=[sensMatrix;sens_all]; %add a control simulation for the drug+cytokine with no knockout

           %''''''' SENSITIVITY ANALYSIS OF SEQUENTIAL NODE KD ''''''
            for u = 1:107
            rparNew = [wNew;n;EC50;doseNew;drugBindingNew;drugAgonismNew];
            paramsNew = {rparNew,tau,ymax,speciesNames,KI};
            ymaxNew=ymax;
            ymaxNew(u)=0; %knockout node
            paramsNew = {rparNew,tau,ymaxNew,speciesNames,KI}; % Only need if doing a knockout
            [t2,y2] = ode15s(@fib617_modelODE_drug_static,tspan,y0,options,paramsNew);
            sens_all=y2(end,:);
            %ctr=ctr+1;
           
            sensMatrix=[sensMatrix;sens_all];
            end
           
            
            
            
        end   
    end
end


%% Convert sens analysis to table for exporting to CSV


cytokines={'control','IL1','TGFB','TGFB_IL1'};
drugNames={'Anakinra','Marimistat','Galunisertib',...
    'Fasudil','Pirfenidone','WH4023','HY12289A','Glutathione',...
    'SB203580','Valsartan','Val_BNP','BNP','CWHM12','Salbutamol'};
specNames2={'Control',speciesNames{:}};
colNames={'Treatment',specNames2{:}};
ctr=1;
treatLabels={};
for i=1:length(cytokines)
    treatLabels{ctr}=cytokines{i};
    ctr=ctr+1;
    for j=1:length(drugNames)
        
        if strcmp(cytokines{i},'control') | strcmp(cytokines{i},'TGFB_IL1')
        treatLabels{ctr}=append(drugNames{j},'_2_',cytokines{i})
        else
        treatLabels{ctr}=append(drugNames{j},'_3_',cytokines{i})   
        end
        ctr=ctr+1;
    end
end

treatLabels=repelem(treatLabels,108);
%%
specNames2={'Control',speciesNames{:}};
%treatLabels=repelem(treatLabels,107);
node=repmat(specNames2,1,60);

T=array2table(sensMatrix,'VariableNames',speciesNames);
T=[table(node','VariableNames',{'Node'}), T];%concatenate tables
T=[table(treatLabels','VariableNames',{'Drug'}), T];%concatenate tables

%%
writetable(T,'DrugSimulations_KDsens.csv');


%% Output csv for python regression

cytokines={'control','IL1','TGFB','TGFB_IL1'};
drugNames={'Anakinra','Marimistat','Galunisertib',...
    'Fasudil','Pirfenidone','WH4023','HY12289A','Glutathione',...
    'SB203580','Valsartan','Val_BNP','BNP','CWHM12','Salbutamol'};
colNames={'Treatment',speciesNames{:}};
ctr=1;
treatLables={};
for i=1:length(cytokines)
    treatLabels{ctr}=cytokines{i};
    ctr=ctr+1;
    for j=1:length(drugNames)
        
        if strcmp(cytokines{i},'control') | strcmp(cytokines{i},'TGFB_IL1')
        treatLabels{ctr}=append(drugNames{j},'_2_',cytokines{i})
        else
        treatLabels{ctr}=append(drugNames{j},'_3_',cytokines{i})   
        end
        ctr=ctr+1;
    end
end
%%
T=array2table(sens_all,'VariableNames',speciesNames);
T=[table(treatLabels','VariableNames',{'Treatment'}), T];%concatenate tables
writetable(T,'SimulatedDrugScreen.csv');




