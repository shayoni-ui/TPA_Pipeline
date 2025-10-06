function [inputs] = initializePBPKInputStruct(nva)
arguments
    nva.phys; % physiology structure from Gastroplus
    nva.acat; % acat structure from Gastroplus
    nva.funKbas; % Kbas function
    nva.coeffKbas; % Kbas coefficients
    nva.species; % Species
end
% initializePBPKInputStruct  ...
% Parameters:
% 	 phsy; % physiology strucuture from Gastroplus
%    acat; % acat structure from Gastroplus
%    funKbas; % function for Kbas
%    coeffKbas; % coefficients for Kbas
%    species; {'Human','Minipig','Dog',
%              'Rat','Mouse','Monkey'}
% Outputs
% 	 inputs: PBPK input strucutre
% Date modified: 04-Oct-2022
% File modified by SMTB user: jbp17697
% File created by SMTB user
% Examples:
%

%persistent phys acat kbas_coeffs

if isempty(fieldnames(nva))
    pkgPath = fileparts(which('+smtb/+pbpk/initializePBPKInputStruct.m'));
    load(fullfile(pkgPath,'matFiles/GastroplusPhysiology.mat'), ...
        'pbpk','acat')
    load(fullfile(pkgPath, 'matFiles/k_bas coefficients.mat'), ...
        'fun_kbas','coeffs');
    phys = pbpk;
    kbasCoeffs{1} = fun_kbas;
    kbasCoeffs{2} = coeffs;
    species = {'Human','Minipig','Dog','Rat',...
        'Mouse','Monkey'};
elseif size(fieldnames(nva), 1) == 5
    phys = nva.phys;
    acat = nva.acat;
    kbasCoeffs{1} = nva.funKbas;
    kbasCoeffs{2} = nva.coeffKbas;
    species = nva.species;
else
    error(['Please pass either no arguments or set all ' ...
        'phys, acat, funkbas, coeffkbas and species arguments' ])
end


inputs.compound             = '';
inputs.parent               = '';
inputs.smiles               = '';
inputs.PKModel              = 'PBPK';
inputs.Tags                 = '';
inputs.mw                   = [];
inputs.logP                 = [];
inputs.logP_source          = 'default';
inputs.logD                 = [];
inputs.logD_pH              = 7.4; 
% Initialize logD as a logP value (pH = -1)
inputs.logD_source          = 'default';
inputs.logD_pH74            = [];
inputs.Acidic_pKa           = [];
inputs.Basic_pKa            = [];
inputs.Acidic_pKa_source    = 'default';
inputs.Basic_pKa_source     = 'default';
inputs.Deff                 = 0;
inputs.Deff_source          = 'default';
inputs.Deff_method          = 'default';   
% Adjusts Deff in GI to account for bile salt effect
inputs.sol_mgmL             = [];
inputs.sol_pH               = [];
inputs.sol_source           = 'default';
inputs.sol_matrix           = '';
inputs.SF                   = 0;
inputs.SF_source            = 'default';
inputs.bio_sol              = struct('matrix', ...
    {'FaSSGF';'FaSSIF';'FeSSIF'},'pH',{1.6;6.5;5}, ...
    'value',cell(3,1));
inputs.FaSSGF               = [];
inputs.FaSSIF               = [];
inputs.FeSSIF               = [];
inputs.FaSSGF_pH            = [];
inputs.FaSSIF_pH            = [];
inputs.FeSSIF_pH            = [];
inputs.FaSSGF_source        = 'default';
inputs.FaSSIF_source        = 'default';
inputs.FeSSIF_source        = 'default';
inputs.SR                   = [];
inputs.SR_method            = 'default';  
% method to calculate solubilization ratio for solubility in 
% presence of bile calculation
inputs.Kp_method_perf       = 'Lukacova';
inputs.Kp_method_perm       = 'Poulin Extracellular';
inputs.fut_type             = 'S+v9.5';
inputs.ASF_method           = 'Opt LogD SAV 6.1';
inputs.fup_type             = 'adjusted';
inputs.para_abs             = 'default';
inputs.Peff0                = [];
inputs.Peff_source          = 'default';
inputs.r_he                 = 0;
inputs.r_s                  = 0;
inputs.r                    = 25;
inputs.h_max                = 30;
inputs.SpeciesSpecific      = struct('Species',species,...
    'fup',              1,...
    'adjusted_fup',     1,...
    'fup_type',         'adjusted',...
    'fup_source',       'default',...
    'fup_source_n',     [],...
    'fup_comments',     [],...
    'B2P',              1,...
    'B2P_source',       'default',...
    'B2P_source_n',     [],...
    'B2P_comments',     [],...
    'fub',              [],...
    'fub_source',       'default',...
    'fub_source_n',     [],...
    'fub_comments',     [],...
    'Tissues',          struct('Name',{},'Model',{},'Kp',{}, ...
                               'fut',{},'Kp_method',{}),...
    'Kp_perf',          [],...
    'Kp_perm',          [],...
    'CLint',            [],...
    'CLint_units',      'milliliter/minute/g tissue',...
    'CLint_source',     'default',...
    'CLint_fuinc',      1,... % Default assume fuinc = 1
    'fuinc_source',     'default',...
    'CLint_source_n',   [],...
    'CLint_comments',   [],...
    'RenalFiltration',  'fup*GFR',...
    'CLint_u_Lhr',      [],...
    'CL_hep_pl',        [],...
    'CL_ren_pl',        [],...
    'CL_tot_pl',        [],...
    'CL_tot_bl',        [],...
    'CL_units',         [],...
    'CLh_pctLBF',       [],...
    'CL_source',        'calculated from CLint',...
    'fu_ent',           1,...
    'gut_fpe',          0,...
    'Peff',             0,...
    'dose_vol',         0);
inputs.DrugProps    = smtb.DrugProps.empty;
inputs.PKpredictor  = smtb.PKpredictor.empty;
inputs.Physiologies = phys;
inputs.ACAT         = acat;
inputs.kbas_coeffs  = kbasCoeffs;
inputs.Pulmonary    = false;
inputs.Distribution = 'PBPK';
inputs.calc_opt = struct('iv_data', 'all',...
    'param_op', 'arithmetic mean',...
    'pk_data_op', 'geometric mean');

end



