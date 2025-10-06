function [inputs, v] = getPBPKVariant(dpArg, nva)
arguments
    dpArg=[];%                                               ARG 2
    nva.species = {'human'};  %                              ARG 3
    nva.Tags char = ''; %                                    ARG 4
    nva.Distribution char = 'PBPK';%                         ARG 5
    nva.Pulmonary logical = false;%                          ARG 6
    nva.compound char = ''; % compound number                ARG 7
    nva.parent char = ''; % parent , gcn,                    ARG 8
    nva.smiles char = ''; %smiles                            ARG 9
    nva.mw = []; %molecular weight                           ARG 10
    nva.logP = []; % logP can be a value or the source info  ARG 11
    nva.logD = []; % logD can be a value or the cource info  ARG 12
    nva.logD_pH double = 7.4;%                               ARG 13
    nva.sol_pH = [];%                                        ARG 14
    nva.sol_matrix char ='';%                                ARG 15
    nva.FaSSGF = [];%                                        ARG 16
    nva.FaSSIF = []; %                                       ARG 17
    nva.FeSSIF = []; %                                       ARG 18
    nva.FaSSGF_pH = [];%                                     ARG 19
    nva.FaSSIF_pH = [];%                                     ARG 20
    nva.FeSSIF_pH = []; %                                    ARG 21
    nva.FaSSGF_source  char = 'default';%                    ARG 22
    nva.FaSSIF_source char = 'default';%                     ARG 23
    nva.FeSSIF_source char = 'default';%                     ARG 24
    nva.Acidic_pKa = [];%                                    ARG 25
    nva.Basic_pKa = [];%                                     ARG 26
    nva.pkasource char = 'default'; %                        ARG 27
    nva.Kp_method_perf char = 'Lukacova';%                   ARG 28
    nva.Kp_method_perm char = 'Poulin Extracellular';%       ARG 29
    nva.ASF_method char = 'Opt LogD SAV 6.1';%               ARG 30
    nva.RenalFiltration char = 'fup*GFR';%                   ARG 31
    nva.fup_type char = 'adjusted';%                         ARG 32
    nva.Deff_method char = 'default';%                       ARG 33
    nva.Deff double = 0; %                                   ARG 34
    nva.SF double = 0; %                                     ARG 35
    nva.Peff0 = [];%                                         ARG 36
    nva.SR_method char = 'default';%                         ARG 37
    nva.para_abs char = 'default'; %                         ARG 38
    nva.r_he double = 0; %                                   ARG 39
    nva.r_s double = 0; %                                    ARG 40
    nva.fut_type char = 'S+v9.5'; %                          ARG 41
    nva.r double = 25; %                                     ARG 42
    nva.h_max double = 30; %                                 ARG 43
    nva.param_op char = 'arithmetic mean'; %                 ARG 44
    nva.pk_data_op char = 'geometric mean'; %                ARG 45
    nva.iv_data char = 'all'; %                              ARG 46
    nva.PKpredictor smtb.PKpredictor = smtb.PKpredictor.empty; % ARG 47
    nva.fup = []; %                                          ARG 48 % Sp
    nva.B2P = []; %                                          ARG 49
    nva.CLint = []; %                                        ARG 50
    nva.CLint_fuinc = 1; %                                   ARG 51
    nva.dose_vol double = 0; %                               ARG 52
    nva.fu_ent double = 1; %                                 ARG 53
    nva.gut_fpe double = 0; %                                ARG 54
    nva.CL_tot_bl = []; %                                    ARG 55
    nva.CL_tot_pl = []; %                                    ARG 56
    nva.CL_ren_pl = []; %                                    ARG 57
    nva.fub = []; %                                          ARG 58
    nva.sol_mgmL = []; %                                     ARG 59
    nva.acatFlag = true;%                                    ARG 60
    nva.Kpscale = 1;%                                        ARG 61

end

% generate variant containing parameter information based on passed
% arguments
[inputs,v] = parameterizePBPKvariants(dpArg, nva);

end

function [inputs,vContent] = parameterizePBPKvariants(dpArg, nva)



% Specify default settings 
inputs = smtb.pbpk.initializePBPKInputStruct();

% Assing empty Kp Structures
for i = 1:length(inputs.SpeciesSpecific)
    inputs.SpeciesSpecific(i).Kp_perf = ...
        smtb.pbpk.getCompoundKP(inputs.SpeciesSpecific(i).Species);
    inputs.SpeciesSpecific(i).Kp_perm = ...
        smtb.pbpk.getCompoundKP(inputs.SpeciesSpecific(i).Species);
end

% check for the DrugProps if a real compound
switch class(dpArg)  %ARG:2
    case 'smtb.DrugProps'
        inputs.DrugProps = dpArg;
        inputs.compound  = inputs.DrugProps.name;
        inputs.parent = inputs.DrugProps.gcn;
    case 'struct'
        missingFields = setdiff(fields(inputs),fields(dpArg));
        if ~isempty(missingFields)
            error(['Input structure missing the following ' ...
                'fields: ' sprintf(' "%s" ',missingFields{:})]);
        end
        inputs   = dpArg; % TODO check since it allows to set directly
        % the input struct
        %species  = {inputs.SpeciesSpecific.Species};
    case 'char'
        inputs.DrugProps    = smtb.DrugProps(dpArg);
        inputs.compound     = inputs.DrugProps.name;
        inputs.parent       = inputs.DrugProps.gcn;
    otherwise
        % Empty set used for null parameterization, otherwise throw
        % error
        if ~isempty(dpArg)
            error(['Unrecognized drug property/ compund ' ...
                'number input']);
        end
end

species = nva.species; % ARG 3
if ~iscell(species); species = {species}; end

% Parse user inputs and update the default settings

[inputs,species] = parseInputs(inputs,species,nva);


% Build the species variant
vContent = cell(1, length(species));
for i = 1:length(species)
    [vnn, inputs] = buildSpeciesVariants(species{i},inputs,...
        nva.acatFlag);
    vContent{i} = vnn;
end
vContent = [vContent{:}];

% Only return the species specific properties for the requrested species
[~,idx] = intersect(lower({inputs.SpeciesSpecific.Species}), ...
    lower(species));
inputs.SpeciesSpecific = inputs.SpeciesSpecific(idx);

end
