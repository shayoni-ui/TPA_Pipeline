
function [inputs, species] = parseInputs(inputs,species,nva)
% Helper function to parse inputs

inputs.Tags = nva.Tags; % ARG 4=''
inputs.Distribution = nva.Distribution; % ARG 5='PBPK'
inputs.Pulmonary = nva.Pulmonary; % ARG 6=false
if isempty(inputs.compound)
    inputs.compound = nva.compound; % ARG 7=''
    inputs.parent = nva.parent; % ARG 8=''
end

inputs.smiles = nva.smiles; % ARG 9=''
inputs.mw = nva.mw; % ARG 10- []

if isnumeric(nva.logP) % ARG 11
    if ~isempty(nva.logP)
        inputs.logP = nva.logP;
        inputs.logP_source = "User Defined";
    end
elseif ischar(inputs.logP)
    inputs.logP_source = nva.logP;
else
    error("logP can be of numeric or char type");
end

if isnumeric(nva.logD) % ARG 12 -[]
    if ~isempty(nva.logD)
        inputs.logD = nva.logD;
        inputs.logD_source = "User Defined";
    end
elseif ischar(inputs.logD)
    inputs.logD_source = nva.logD;
else
    error("logD can be of numeric or char type");
end

inputs.logD_pH = nva.logD_pH; % ARG 13 - 7.4

if isnumeric(nva.sol_mgmL) % ARG 59 - []
    if ~isempty(nva.sol_mgmL)
        inputs.sol_mgmL = nva.sol_mgmL;
        inputs.sol_source = 'User Defined';
    end
elseif ~strcmpi(nva.sol_mgmL, 'default')
    inputs.sol_source = nva.sol_mgmL;
end

if isnumeric(nva.sol_pH) % ARG 14=[]
    if ~isempty(nva.sol_pH)
        inputs.sol_pH = nva.sol_pH;
    end
elseif strcmpi(nva.sol_pH, 'default')
    inputs.sol_pH = [];
else
    error('Solubility pH must be input as default or numeric')
end

inputs.sol_matrix = nva.sol_matrix; % ARG 15 ='', char

if isnumeric(nva.FaSSGF) % ARG16 =[], 'QSAR ADMET'
    if ~isempty(nva.FaSSGF)
        inputs.FaSSGF = nva.FaSSGF;
        inputs.FaSSGF_source = 'User Defined';
    end
else
    inputs.FaSSGF_source = nva.FaSSGF;

end

if isnumeric(nva.FaSSIF) % ARG 17 = []
    if ~isempty(nva.FaSSIF)
        inputs.FaSSIF = nva.FaSSIF;
        inputs.FaSSIF_source = 'User Defined';
    end
else
    inputs.FaSSIF_source = nva.FaSSIF;
end

if isnumeric(nva.FeSSIF) % ARG 18 = []
    if ~isempty(nva.FeSSIF)
        inputs.FeSSIF = nva.FeSSIF;
        inputs.FeSSIF_source = 'User Defined';
    end
else
    inputs.FeSSIF_source = nva.FeSSIF_source;
end

if isnumeric(nva.FaSSGF_pH) % ARG 19 = []
    if ~isempty(nva.FaSSGF_pH)
        inputs.FaSSGF_pH = nva.FaSSGF_pH;
    end
elseif ~strcmpi(nva.FaSSGF_pH, 'default')
    error('FaSSGF_pH expecting default or numeric as input');
end

if isnumeric(nva.FaSSIF_pH) % ARG 20 = []
    if ~isempty(nva.FaSSIF_pH)
        inputs.FaSSIF_pH = nva.FaSSIF_pH;
    end
elseif ~strcmpi(nva.FaSSIF_pH, 'default')
    error('FaSSIF_pH expecting default or numeric as input');
end

if isnumeric(nva.FeSSIF_pH) % ARG 21 = []
    if ~isempty(nva.FeSSIF_pH)
        inputs.FeSSIF_pH = nva.FeSSIF_pH;
    end
elseif ~strcmpi(nva.FeSSIF_pH, 'default')
    error('FaSSIF_pH expecting default or numeric as input');
end

%inputs.FaSSGF_source = nva.FaSSGF_source; % ARG 22 = 'default'
%inputs.FaSSIF_source = nva.FaSSIF_source; % ARG 23 = 'default'
%inputs.FeSSIF_source = nva.FeSSIF_source; % ARG 24 = 'default'

if isnumeric(nva.Acidic_pKa) % ARG 25 = [];
    if ~isempty(nva.Acidic_pKa)
        inputs.Acidic_pKa = nva.Acidic_pKa;
        inputs.Acidic_pKa_source = 'User Defined';
    end
else
    inputs.Acidic_pKa_source = nva.Acidic_pKa; % TODO check what
    % possible options
end

if isnumeric(nva.Basic_pKa) % ARG 26 = [];
    if ~isempty(nva.Basic_pKa)
        inputs.Basic_pKa = nva.Basic_pKa;
        inputs.Basic_pKa_source = 'User Defined';
    end
else
    inputs.Basic_pKa_source = nva.Basic_pKa; % TODO check what possible
    % options
end

if ~strcmpi(inputs.Acidic_pKa_source, 'User Defined') % ARG 27 = 'default'
    inputs.Acidic_pKa_source = nva.pkasource;
end
if ~strcmpi(inputs.Basic_pKa_source, 'User Defined')
    inputs.Basic_pKa_source = nva.pkasource;
end

inputs.Kp_method_perf = nva.Kp_method_perf;  % ARG 28 = 'Lukacova'
inputs.Kp_method_perm = nva.Kp_method_perm;% ARG 29 ='Poulin Extracellular'
inputs.ASF_method = nva.ASF_method; % ARG 30 = 'Opt LogD SAV 6.1'

for ii = 1:length(inputs.SpeciesSpecific)       % ARG 31 = 'fup*GFR'
    inputs.SpeciesSpecific(ii).RenalFiltration = nva.RenalFiltration;
end

inputs.fup_type = nva.fup_type;       % ARG 32 = 'adjusted'
inputs.Deff_method = nva.Deff_method; % ARG 33 = 'default'

inputs.Deff = nva.Deff; % ARG 34 = 0
if nva.Deff ~= 0; inputs.Deff_source = 'User Defined'; end

inputs.SF = nva.SF; % ARG 35 =0
if nva.SF ~=0; inputs.SF_source = 'User Defined'; end

if isnumeric(nva.Peff0) % ARG 36 = [];
    if ~isempty(nva.Peff0)
        inputs.Peff0 = nva.Peff0;
        inputs.Peff_source = 'User Defined';
    end
else
    inputs.Peff_source = nva.Peff0; % TODO check what possible options
end

%'Solubilization Ratio Method'
inputs.SR_method = nva.SR_method; % ARG 37 = 'default'

%'Paracellular Absorption Adson/Zhimin/off'
inputs.para_abs = nva.para_abs; % ARG 38 = 'default'

%Hydrodynamic radius of molecular ellipsoid (r_he)
inputs.r_he = nva.r_he; % ARG 39 = 0

% Average projected geometric radius (r_s)
inputs.r_s = nva.r_s; % ARG 40 = 0
inputs.fut_type = nva.fut_type; % ARG 41 = 'S+v9.5'

% Particle size micro meter
inputs.r = nva.r; % ARG 42  =25

% Maximum diffusion layer thickness (micrometer)
inputs.h_max = nva.h_max; % ARG 43 = 30

inputs.calc_opt.param_op = nva.param_op; % ARG 44 = 'arithmetic mean'
inputs.calc_opt.pk_data_op = nva.pk_data_op; % ARG 45 = 'gemoetric mean'

% IV data selection for allometry and clearance definition
inputs.calc_opt.iv_data = nva.iv_data; % ARG 46 = 'all'

inputs.PKpredictor = nva.PKpredictor; % ARG 47 = PKpredictor.empty

% Parse species specific inputs:
% TODO: Change line to find multiple species. strcmpi wont work when
% species is an array

%idx = strcmpi({inputs.SpeciesSpecific.Species}, species);
idx = contains(lower({inputs.SpeciesSpecific.Species}), species);

if isnumeric(nva.fup) % ARG 48 = [];
    if ~isempty(nva.fup)
        [inputs.SpeciesSpecific(idx).fup] = deal(nva.fup);
        [inputs.SpeciesSpecific(idx).fup_source] = deal('User Defined');
    end
else
    [inputs.SpeciesSpecific(idx).fup_source] = deal(nva.fup);
end

if isnumeric(nva.B2P) % ARG 49 = 1;
    if ~isempty(nva.B2P)
        [inputs.SpeciesSpecific(idx).B2P] = deal(nva.B2P);
        [inputs.SpeciesSpecific(idx).B2P_source] = deal('User Defined');
    end
else
    [inputs.SpeciesSpecific(idx).B2P_source] = deal(nva.B2P);
end


if isnumeric(nva.CLint) % ARG 50 = [];
    if ~isempty(nva.CLint)
        [inputs.SpeciesSpecific(idx).CLint] = deal(nva.CLint);
        [inputs.SpeciesSpecific(idx).CLint_source] = deal('User Defined');
    end
else
    [inputs.SpeciesSpecific(idx).CLint_source] = deal(nva.CLint);
end

if isnumeric(nva.CLint_fuinc) % ARG 51 = 1;
    [inputs.SpeciesSpecific(idx).CLint_fuinc] = deal(nva.CLint_fuinc);
    if nva.CLint_fuinc ~=1
        [inputs.SpeciesSpecific(idx).fuinc_source] = deal('User Defined');
    end
else
    [inputs.SpeciesSpecific(idx).fuinc_source] = deal(nva.CLint_fuinc);
end

[inputs.SpeciesSpecific(idx).dose_vol] = deal(nva.dose_vol); % ARG 52 = 0
[inputs.SpeciesSpecific(idx).fu_ent] = deal(nva.fu_ent); % ARG 53 = 1
[inputs.SpeciesSpecific(idx).gut_fpe] = deal(nva.gut_fpe); % ARG 54 = 0


if isnumeric(nva.CL_tot_bl) % ARG 55 = [];
    if ~isempty(nva.CL_tot_bl)
        [inputs.SpeciesSpecific(idx).CL_tot_bl]= deal(nva.CL_tot_bl);
        [inputs.SpeciesSpecific(idx).CL_source] =...
            deal('User Defined (blood)');
    end
else
    [inputs.SpeciesSpecific(idx).CL_source] = deal(nva.CL_tot_bl);
end


if isnumeric(nva.CL_tot_pl) % ARG 56 = [];
    if ~isempty(nva.CL_tot_pl)
        [inputs.SpeciesSpecific(idx).CL_tot_pl]= deal(nva.CL_tot_pl);
        [inputs.SpeciesSpecific(idx).CL_source] =...
            deal('User Defined (plasma)');
    end
else
    [inputs.SpeciesSpecific(idx).CL_source] = deal(nva.CL_tot_pl);
end

if ~isempty(nva.CL_ren_pl) % ARG 57 = [];
    [inputs.SpeciesSpecific(idx).CL_ren_pl] = deal(nva.CL_ren_pl);
    [inputs.SpeciesSpecific(idx).RenalFiltration] = deal('User Defined');
end

if isnumeric(nva.fub) % ARG 58 = [];
    if ~isempty(nva.fub)
        [inputs.SpeciesSpecific(idx).fub]= deal(nva.fub);
        [inputs.SpeciesSpecific(idx).fub_source] = deal('User Defined');
    end
else
    [inputs.SpeciesSpecific(idx).fub_source] = deal(nva.fub);
end

[inputs, species] = parseInputsFromDrugProps(inputs, species, nva);

end


