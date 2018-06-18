%% Pannala et al., (2018)
% Script to run timbr predictions for acetaminophen induced liver injury
% Under metabolic flux analysis constraints.
% Load the path to reactions weights generated using the R scripts
path_rcode_directory = '(location of the reaction weights)\Weights_apap/';
path_model_directory = '(current model directory)/';
% Specify output directory for TIMBR predictions
path_timbr_directory = [path_model_directory 'MFA_5h/'];

changeCobraSolver('glpk'); % or glpk
% % Load COBRA models
load iRno_v2.mat;
obj_value_required = 0;
% atp_value_required = 5000;
atp_value_required = 0;
% Set unconstrained reaction bounds to +/- 10^5 because the largest
% consumption rate is fairly close to the default 10^3 reaction bound.
default_bound = 10000; % 10^5
rno_irrev_base = ncomm_blais_model2irrev(rno_cobra);
default_rxns_rno = rno_irrev_base.rxns(rno_irrev_base.ub == 1000);
rno_irrev_default = changeRxnBounds(rno_irrev_base, ...
    default_rxns_rno, default_bound, 'u');
rno_irrev_atp = changeRxnBounds(rno_irrev_default,...
    'RCR11017_f',atp_value_required,'l');
rno_irrev_biomass = changeRxnBounds(rno_irrev_atp,...
    'RCR99999_f',obj_value_required,'l');
% Note: as mentioned above, TIMBR predictions were not performed while 
% requiring biomass production.
rno_irrev1 = rno_irrev_biomass;
%% Bounds for metabolites as input/uptake
ex_u = {'acetoacetate exchange (rev)','alanine exchange (rev)','arginine exchange (rev)', ...
        'asparagine exchange (rev)','glutamine exchange (rev)','glycine exchange (rev)', 'histidine exchange (rev)','L-lactate exchange (rev)',...
        'lysine exchange (rev)','methionine exchange (rev)','ornithine exchange (rev)','phenylalanine exchange (rev)', 'proline exchange (rev)',...
        'serine exchange (rev)','threonine exchange (rev)','tyrosine exchange (rev)','valine exchange (rev)','glycerol exchange (rev)', ...
        'cysteine exchange (rev)','glucose exchange (rev)','isoleucine exchange (rev)','leucine exchange (rev)'};
ex_uv =1.0*[0.0, 111, 30.76,...
        17.1, 17.1, 0.5*46.15, 68.37,0.5*486,... %840
        12.82,6.83, 29.06, 4.27, 1*5.98, ... %21
        23.06, 10.25, 8.54, 6.83, 260,...
        14.2,0, 7.65, 13.21]; % Izamis et al 2010
for i =1:length(ex_u)
    ex_ub(i) = find(strcmp(ex_u(i),rno_irrev1.rxnNames));
end

% % Bounds for metaboolites with no data (fixed based on invitro)
ex_u1 = {'linolenate exchange (rev)','cholesterol exchange (rev)', ...
        'linoleate exchange (rev)','palmitate exchange (rev)', 'oleate exchange (rev)', ...
        'stearate exchange (rev)','palmitolate exchange (rev)','myristic acid exchange (rev)'};
ex_uv1 = 0.5*[6.2,47.84,11.64,11.64,15.13,6.2,7.8,4.66];
for i =1:length(ex_u1)
    ex_ub1(i) = find(strcmp(ex_u1(i),rno_irrev1.rxnNames));
end

% % Select bounds for metabolites which get secreted under gluconeogenic
% conditions (glucose,uea, ketone bodies (acetoacetate and
% B-hydroxybutyrate))
ex_s = {'glucose exchange (fwd)'};
ex_sv = [140];% assumed a random mimimum requirement to follow MFA constraints
for i =1:length(ex_s)
    ex_ub2(i) = find(strcmp(ex_s(i),rno_irrev1.rxnNames));
end
% constrain the model based on above inputs
rno_irrev = changeRxnBounds(rno_irrev1,rno_irrev1.rxns(ex_ub),ex_uv,'b');
rno_irrev = changeRxnBounds(rno_irrev,rno_irrev.rxns(ex_ub1),ex_uv1,'b');
rno_irrev = changeRxnBounds(rno_irrev,rno_irrev.rxns(ex_ub2),ex_sv,'l');

%%      G6Pase    ,  GPI ,      , ALDO        F1,6Bpase  ,  GADPH     , G3P        , PGK      ,     ENO       
gnf = {'RCR11028_f','RCR14185_f','RCR14183_f','RCR14184_f' 'RCR11027_f','RCR11051_f','RCR10230_r','RCR10443_r'...
    'RCR10289_f','RCR11087_r','RCR11072_f','RCR10283_f','RCR11661_r','RCR10349_f','RCR11671_f','RCR11373_f','RCR11674_r'};
% PK        ,    LDH(1337)  , PC        , PEPCK (1319) , CS       , IDH (466)  , OGDH       , PCC        ,  SDH
for k=1:length(gnf)
    gr(k)=find(strcmp(gnf(k),rno_irrev.rxns));
end

% Obtained from the literature (Figure 4b in the main text)
v_lb = 0.9*[378,378,378,378,505.28,252.64,505.28,505.28, ...
        246.58,381.64,628.21,751.85,357.08,357.08,357.08,123.64,480.77];
v_ub = 1.3*[378,378,378,378,505.28,252.64,505.28,505.28, ...
        246.58,381.64,628.21,751.85,357.08,357.08,357.08,123.64,480.77];
% % Flux data for control
rno_irrev.lb(gr)=v_lb;
rno_irrev.ub(gr) = v_ub; 
rno_irrev.ub([308,572,1268,1336,2059,2079,5092,5095])=0;% to get the same result,these reactions need to be blocked as they are non-operational in the reversible network if fluxes are fixed for mfa

%% Determine maxmimum flux through each exchange reaction for iRno
[rno_exchange, rno_demand] = findExcRxns(rno_irrev);
% remove O2 exchange and glycogenin exchanges 
O2f = find(strcmp('O2 exchange (fwd)',rno_irrev.rxnNames));
rno_exchange([O2f])=0;
rno_production = rno_irrev.rxns((cellfun(@length, ...
    regexp(rno_irrev.rxns,'_r')) == 0) & rno_exchange);
rno_fva = zeros(length(rno_production),1);
for exchange_index = 1:length(rno_production)
    rno_consumption = setdiff(intersect(...
        regexprep(rno_production(exchange_index),'_f','_r'),...
        rno_irrev.rxns),...
        rno_production(exchange_index));
    if (~isempty(rno_consumption))
        if (strcmp(rno_consumption, 'RCR30327_r')) % glycerol exchange reversible =RCR30327_r
            rno_production_fba = optimizeCbModel(changeObjective(...
            changeRxnBounds(changeRxnBounds(...
            rno_irrev,rno_production,default_bound,'u'),...
            rno_consumption,0,'l'),...
            rno_production(exchange_index)));
         rno_fva(exchange_index,1) = rno_production_fba.f;
        else
            rno_production_fba = optimizeCbModel(changeObjective(...
            changeRxnBounds(changeRxnBounds(...
            rno_irrev,rno_production,default_bound,'u'),...
            rno_consumption,0,'b'),...
            rno_production(exchange_index)));
        rno_fva(exchange_index,1) = rno_production_fba.f;
        end
    else 
        rno_production_fba = optimizeCbModel(changeObjective(...
            changeRxnBounds(...
            rno_irrev,rno_production,default_bound,'u'),...
            rno_production(exchange_index)));
        rno_fva(exchange_index,1) = rno_production_fba.f;
          
    end
       if rno_production_fba.stat ~=1
            rno_production(exchange_index)
            disp('solution not feasible')
           %pause
       end
    
    disp(rno_production_fba)
end

%% Specify minimimum required flux for each exchange reaction

fba_threshold = 1e-4;
rno_production_ok = rno_fva > fba_threshold;
rno_production_id = rno_production(rno_production_ok,1);
rno_production_requirement = min(0.90 * rno_fva(rno_production_ok,1), 100);
% rno_production_requirement = 0.90*rno_fva(rno_production_ok,1);
rno_production_count = length(rno_production_id);
% % Load TIMBR reaction weights
rno_timbr_weights_load = readtable([path_rcode_directory, ...
    'apap_timbr_weights_rno_liver5h.txt'],'Delimiter','\t');

[rno_model_has_weight, rno_model2weight_index] = ismember(...
    rno_irrev.rxns, rno_timbr_weights_load.rxn_irreversible);
[rno_weight_in_model, rno_weight2model_index] = ismember(...
    rno_timbr_weights_load.rxn_irreversible, rno_irrev.rxns);

if (all(rno_model_has_weight))
    % first 3 columns are organism_id, rxn_id, and rxn_irrerversible
    rno_timbr_weights = rno_timbr_weights_load(rno_model2weight_index,4:end);
    rno_timbr_id = rno_timbr_weights.Properties.VariableNames';
    rno_timbr_count = length(rno_timbr_id);

else
    warning('reaction weights not specified for all reactions in rno')
end

% % Run TIMBR algorithm and save results as .txt files
% Estimate the global network demand of producing each exchange metabolite

for rno_timbr_index = 1:rno_timbr_count
    if rno_timbr_index==1
        % Flux data for control
v_lb = 0.9*[378,378,378,378,505.28,252.64,505.28,505.28, ...
        246.58,381.64,628.21,751.85,357.08,357.08,357.08,123.64,480.77];
v_ub = 1.3*[378,378,378,378,505.28,252.64,505.28,505.28, ...
        246.58,381.64,628.21,751.85,357.08,357.08,357.08,123.64,480.77];
    
rno_irrev.lb(gr)=v_lb;
rno_irrev.ub(gr) = v_ub;

    else    %         % Flux data for treatment with acetaminophen
disp('treatment') 
rno_irrev.lb(gr) = v_lb; 
rno_irrev.ub(gr) = 1.3*v_ub;  % not strictly enforced for the treatment
    end
    rno_timbr_network_demand = zeros(rno_timbr_count, 1);
    disp(rno_timbr_id{rno_timbr_index});
    for rno_production_index = 1:rno_production_count
        disp(rno_production_id(rno_production_index,1));
        rno_timbr_network_demand(rno_production_index,1) = ...
            timbr(rno_irrev, ...
            rno_production_id(rno_production_index,1), ...
            rno_production_requirement(rno_production_index,1), ...
            rno_timbr_weights(:,rno_timbr_index));
    end
    rno_timbr_file_name = [path_timbr_directory ...
        'ncomm_blais_timbr_' rno_timbr_id{rno_timbr_index} '.txt'];
    
    rno_timbr_table =  table(...
        rno_production_id,...
        rno_timbr_network_demand,...
        'VariableNames',{'timbr_id' 'timbr_value'});
    writetable(rno_timbr_table,rno_timbr_file_name,...
        'Delimiter','\t');
end

Cn = readtable('(selected location of the generated files above\MFA_5h\timbr_rno_acetaminophen_t1_d1_ctl.txt');
Tr = readtable('(selected location of the generated files above\MFA_5h\timbr_rno_acetaminophen_t1_d1_trt.txt');

% calculate the raw timbr score using the netwrok demand files generated
% above
Raw_score = -(Tr{:,2}-Cn{:,2})./(Tr{:,2} + Cn{:,2});
Prod_score = (Raw_score-mean(Raw_score))./std(Raw_score);
s1 = Cn{:,1}; % rxn ids of rno_irrev

for i = 1:length(s1)
    rn(i) = find(strcmp(s1(i),rno_irrev.rxns));
end
 
rnames = rno_irrev.rxnNames(rn);
Table_5h = [rnames num2cell(Prod_score)];
