function [asf,sol,Deff,P_para,ASF_para,SR,k_bas,...
    inputs] = calcGIprops(cmpt,species,state,nva)
arguments
    cmpt;
    species;
    state;
    nva.dp smtb.DrugProps = smtb.DrugProps.empty;
    nva.logP double = [];
    nva.Acidic_pKa double = [];
    nva.Basic_pKa double = [];
    nva.mw = [];
    nva.bio_sol struct = struct.empty; % biorelevant solubilities
    nva.ASF_method char = '';
    nva.Deff_method char = '';
    nva.Deff double = [];
    nva.r_he double = []; % hydrodynamically equivalent sphere radius
    nva.r_s double = []; %projected radius of solute molecule
    nva.SR_method char = ''; % solubilization ratio method
    nva.SR double  = []; % solubilization bile ratio
    nva.RefSol double = []; % intrinsic solubility, reference solubility
    nva.RefSolpH double = [];
    nva.SF double = [];
    nva.ACAT struct = strcut.empty;
    nva.kbas_coeffs cell = {};
    nva.para_model char = ''; %'paracellular model (Zhimin/Adson/off)'
    nva.Peff double = []; %
end
%CALCGIPROPS calculates absorption scale factor and GI bile salt solubilty
%for a given compound
%   [asf,sol,Deff,inputs] = calcGIprops(cmpt,species,state) calculates asf,
%   solubility, and diffusion coefficient in the given compartment,
%   species, and state (not reccomended since this will use meaningless
%   inputs)
%
%   [asf,sol,Deff,inputs] = calcGIprops(...,drugprops), uses a drugprops
%   object to speciefy the needed inputs
%
%   [asf,sol,Deff,inputs] = calcGIprops(...,field,value) allows setting of
%   inputs in field, value pairs
%
%   [asf,sol,Deff,inputs] = calcGIprops(...,drugprops,field,value), uses a
%   drugprops object to set inputs and overrides with values given in
%   field,value pairs
%   Arguments:
%         cmpt;
%         species;
%         state;
%         nva.dp DrugProps = DrugProps.empty;
%         nva.logP double = [];
%         nva.Acidic_pKa double = [];
%         nva.Basic_pKa double = [];
%         nva.mw = [];
%         nva.bio_sol struct = struct.empty; % biorelevant solubilities
%         nva.ASF_method char = ''; Use "Theoretical", "Theoretical SAV",
%            "logD","Opt LogD","Opt LogD SAV", or "Opt LogD SAV 6.1"
%         nva.Deff_method char = '';
%         nva.Deff double = [];
%         nva.r_he double = []; % hydrodynamically equivalent sphere radius
%         nva.r_s double = []; %projected radius of solute molecule
%         nva.SR_method char = ''; % solubilization ratio method
%         nva.SR double  = []; % solubilization bile ratio
%         nva.RefSol double = []; intrinsic solubility,reference solubility
%         nva.RefSolpH double = [];
%         nva.SF double = [];
%         nva.ACAT strcut = strcut.empty;
%         nva.kbas_coeffs cell = {};
%         nva.para_model char = ''; 'paracellular model (Zhimin/Adson/off)'
%         nva.Peff double = [];
%   Outputs:
%


% persistent ACAT




persistent ACAT
if isempty(ACAT)
    pkgPath = fileparts(which('+smtb/+pbpk/initializePBPKInputStruct.m'));
    load(fullfile(pkgPath,'matFiles/GastroplusPhysiology.mat'), ...
        'acat')
    %load('GastroplusPhysiology.mat','acat')
    ACAT = acat; % #ok<NODEF>
else
    acat = ACAT;
end


%kbas_coeffs = [];

asf    = [];
sol    = [];

inputs = struct('logP',         0,...
    'Acidic_pKa',   [],...
    'Basic_pKa',    [],...
    'RefSol',       0,...
    'RefSolpH',    -1,...
    'SF',           0,...
    'mw',           0,...
    'Deff',         0,...
    'Peff',         0,...
    'r_he',         0,...
    'r_s',          0,...
    'ASF_method',   'Opt LogD SA v61',...
    'Deff_method',  'adjusted',...
    'SR_method',    'in vitro',...
    'para_model',   'Zhimin');

if ~isempty(nva.dp)
    inputs.logP         = getPreferredValue(nva.dp,'logP');
    inputs.Acidic_pKa   = getPreferredValue(nva.dp,'pKa_a');
    inputs.Basic_pKa    = getPreferredValue(nva.dp,'pKa_a');
    inputs.RefSol       = nva.dp.sol(strcmpi({nva.dp.sol.matrix}, ...
        'Intrinsic')).value;
    inputs.RefSolpH     = -1;
    inputs.SF           = nva.dp.sf.value;
    inputs.mw           = nva.dp.mw;
    inputs.bio_sol      = nva.dp.sol(strcmpi({nva.dp.sol.matrix}, ...
        'FeSSIF')|...
        strcmpi({nva.dp.sol.matrix},'FaSSIF')|...
        strcmpi({nva.dp.sol.matrix},'FaSSGF')&...
        contains({nva.dp.sol.source},'ADMET','IgnoreCase',true));
    inputs.Deff         = nva.dp.Deff.value;
    inputs.r_he         = nva.dp.descriptors(strcmpi( ...
        {nva.dp.descriptors.type},'r_he')).value;
    inputs.r_s          = nva.dp.descriptors(strcmpi( ...
        {nva.dp.descriptors.type},'r_s')).value;
    inputs.Peff         = nva.dp.perm(strcmpi( ...
        {nva.dp.perm.type},'Peff')&contains( ...
        {nva.dp.perm.source},'QSAR ADMET')).value;
end



if ~isempty(nva.logP); inputs.logP = nva.logP; end
if ~isempty(nva.Acidic_pKa); inputs.Acidic_pKa = nva.Acidic_pKa;end
if ~isempty(nva.Basic_pKa); inputs.Basic_pKa = nva.Basic_pKa; end
if ~isempty(inputs.ASF_method); inputs.ASF_method = nva.ASF_method; end
if ~isempty(inputs.RefSol)
    inputs.RefSol = nva.RefSol;
    inputs.RefSolpH = -1;
end
if ~isempty(nva.RefSolpH); inputs.RefSolpH = nva.RefSolpH; end
if ~isempty(nva.SF); inputs.SF = nva.SF; end
if ~isempty(nva.mw); inputs.mw = nva.mw; end
if ~isempty(nva.bio_sol); inputs.bio_sol = nva.bio_sol; end
if ~isempty(nva.Deff_method); inputs.Deff_method = nva.Deff_method; end
if ~isempty(nva.Deff); inputs.Deff = nva.Deff; end
if ~isempty(nva.r_he); inputs.r_he = nva.r_he; end
if ~isempty(nva.r_s); inputs.r_s = nva.r_s; end
if ~isempty(nva.SR_method); inputs.SR_method = nva.SR_method; end
if ~isempty(nva.SR)
    inputs.SR = nva.SR;
    inputs.SR_method = 'User Defined';
end
if ~isempty(nva.Peff); inputs.Peff = nva.Peff; end
if ~isempty(nva.ACAT); inputs.ACAT = nva.ACAT; end
if ~isempty(nva.kbas_coeffs); inputs.kbas_coeffs = nva.kbas_coeffs; end
if ~isempty(nva.para_model); inputs.para_model = nva.para_model;end



if isempty(nva.kbas_coeffs)
    load('k_bas coefficients.mat','fun_kbas','coeffs');
    kbas_coeffs{1} = fun_kbas;
    kbas_coeffs{2} = coeffs;
else
    kbas_coeffs = nva.kbas_coeffs;
end

if strcmpi(species,'monkey'), species = 'Cynomologous'; end
acat = acat(strcmpi({acat.Species},species) & strcmpi({acat.State},state));
if isempty(acat)
    warning('Unable to find acat physiology for %s %s',state,species);
else
    cmpt = acat.Phys(strcmpi({acat.Phys.Compartment},cmpt));
    if isempty(cmpt), warning('Unable to find GI compartment: %s',cmpt);
    else
        sef_jejunum = acat.Phys(strcmpi({acat.Phys.Compartment}, ...
            'Jejunum1')).SEF;
        asf = calculateASF(cmpt,inputs,sef_jejunum);
        [sol,Deff,SR] = calculateGISolubilityDeff(cmpt,inputs);
        [P_para,ASF_para] = calculateParacellular_ka(cmpt,inputs);
        k_bas = calculate_kbas(cmpt,inputs,species,P_para,asf,kbas_coeffs);
    end
end

end

function asf = calculateASF(cmpt,inputs,sef_jejunum)

switch lower(strrep(strrep(strrep(inputs.ASF_method,' ',''),'_',''), ...
        '.',''))
    case 'theoretical'
        asf = 2/cmpt.Radius;
    case 'theoreticalsav'
        asf = 2/cmpt.Radius * cmpt.SEF / sef_jejunum;
    case 'logd'
        c = [0.19932; 6.26330; 0.27801; 0.23411];
        asf = logDASF(c,cmpt,inputs);
    case 'optlogd'
        c = [0.18148; 1.09440; 0.07340; 0.31836];
        asf = logDASF(c,cmpt,inputs);
    case 'optlogdsav'
        c = [0.88050; 0.01305; 0.09590; 0.70900];
        asf = logDASF(c,cmpt,inputs);
    case 'optlogdsav61'
        c = [0.06944; 0.43028; 0.12147; 0.46632];
        asf = logDASF(c,cmpt,inputs);
    otherwise, error(['Unrecognized ASF method "%s",' ...
            ' please use: "Theoretical", "Theoretical SAV",' ...
            '"logD","Opt LogD","Opt LogD SAV", ' ...
            'or "Opt LogD SAV 6.1"'],inputs.ASF_method)
end
if isempty(asf), asf = 0; end
end

function asf = logDASF(c,cmpt,inputs)

method     = lower(strrep(strrep(strrep(inputs.ASF_method,' ',''), ...
    '_',''),'.',''));
Acidic_pKa = inputs.Acidic_pKa;
Basic_pKa  = inputs.Basic_pKa;
logP       = inputs.logP;
logD       = smtb.pbpk.logDcalc(logP,cmpt.pH,Acidic_pKa,Basic_pKa);
logD_65    = smtb.pbpk.logDcalc(logP,6.5,Acidic_pKa,Basic_pKa);
C          = 6.26; %Gastroplus's fitted constant
% to avoid singular condition

% Calculate SA:V ratio (called A in Gastroplus documentation)
switch method
    case {'logd','optlogd'}
        switch cmpt.Compartment
            case {'Caecum','Colon'}, SAV = 4/cmpt.Radius;
            otherwise,               SAV = 1.2/cmpt.Radius;
        end
    case {'optlogdsav','optlogdsav61'}, SAV = 2/cmpt.Radius * cmpt.SEF;
    otherwise, error('Unrecognized ASF method: %s',inputs.ASF_method);
end

switch strrep(cmpt.Compartment,' ','')
    case {'Duodenum','Jejunum1','Jejunum2','Jejunum3','Jejunum4',...
            'Ileum','Ileum1','Ileum2','Ileum3'}
        switch method
            case 'logd'
                asf = SAV * 10^( c(1) * (logP - logD - c(2)) ...
                    / (logP - logD_65-c(2)));
            otherwise
                asf = SAV * c(2) * 10^(c(1) * (logP - logD - C) ...
                    / (logP - logD_65 - C));
        end
    case 'Stomach'
        asf = 0;
    case {'Caecum','Colon'}
        asf = SAV * c(3) * 10^(c(4)*logD);
end

end

function [sol,Deff,SR] = calculateGISolubilityDeff(cmpt,inputs)

% Solubility calculated using model in Mithani, S.D. et al. Pharm Res.
% 1996, 13(1). pp163-167. Equation 3
% Csx = Cs + SCbs*MW*[bile] where Csx is solubility in presence of bile
% (ug/mL), Cs is solubility in absence of bile (ug/mL), SCbs is
% solubilization capacity, MW is molecular weight, and [bile] is in mM.
% When method is "in vitro," SR is calculated using the measured (or
% predicted from ADMET/QSAR) solubility in the presence of bile (e.g.,
% FaSSIF) using the above equation.  If the method is "theoretical,"
% the theoretical model is used which is eqn 3 in the same reference:
% log(SR) = 2.23 + 0.61*logP

if inputs.SF == 0
    inputs.SF = 20;
end

mw    = inputs.mw;
iFrac_sol = smtb.pbpk.ionfrac(inputs.Acidic_pKa, ...
    inputs.Basic_pKa,inputs.RefSolpH); % Ionization calculations
iFrac_gi = smtb.pbpk.ionfrac(inputs.Acidic_pKa, ...
    inputs.Basic_pKa,cmpt.pH); % Ionization calculations
SF_sol    = inputs.SF .* ones(1,length(iFrac_sol.F_all));
SF_sol(iFrac_sol.NetCharge == 0) = 1;
% Set SF to 0 for neutral species
SF_gi    = inputs.SF .* ones(1,length(iFrac_gi.F_all));
SF_gi(iFrac_gi.NetCharge == 0) = 1;
% Set SF to 0 for neutral species

if inputs.RefSolpH == -1
    iSol = inputs.RefSol;
else
    iSol = inputs.RefSol./min(SF_sol./iFrac_sol.F_all);
end

sol_aq = min(iSol .* SF_gi ./iFrac_gi.F_all);
% Aqueous solubility at GI pH
SC  = (sol_aq*1e-3 / mw) / (1/18);
% Solubilization capacity, mole of drug per mole of water at aqueous pH

switch lower(strrep(strrep(inputs.SR_method,'_',''),' ',''))
    case 'theoretical'
        SR  = 10^(2.23 + 0.61 * inputs.logP);
        % Theoretical bile salt solubilization ratio
        sol = (sol_aq*1e3 + SC * SR * mw * cmpt.Bile)/1e3;
        % Solubility with bile salt effect.  Model output is in ug/mL -
        % the quantity inside () is divided by 1e3 to give mg/mL
    case 'invitro'
        ind = find(strcmpi({inputs.bio_sol.matrix},'FaSSIF')|...
            strcmpi({inputs.bio_sol.matrix},'FeSSIF'));
        SR0 = zeros(1,2);
        sol_bile = zeros(2,1);
        sol_aq_ref = zeros(2,1);
        SR_num = zeros(2,1);
        media = {'FaSSIF','FeSSIF'};
        SR_den = zeros(2,1);
        if isempty(ind)
            warning(['FaSSIF or FeSSIF solubility data not available.' ...
                '  Will calculate SR using theoretical method'])
            SR  = 10^(2.23 + 0.61 * inputs.logP);
            sol = (sol_aq*1e3 + SC * SR * mw * cmpt.Bile)/1e3;
        else
            for i = 1:length(ind)
                mi = media{strcmpi((inputs.bio_sol(ind(i)).matrix),media)};
                switch lower(mi)
                    case 'fassif', bile_conc = 3;
                    case 'fessif', bile_conc = 15;
                end
                iFrac_bile = smtb.pbpk.ionfrac(inputs.Acidic_pKa, ...
                    inputs.Basic_pKa,inputs.bio_sol(ind(i)).pH);
                SF_bile    = inputs.SF .* ones(1,length(iFrac_bile.F_all));
                SF_bile(iFrac_bile.NetCharge == 0) = 1;

                sol_bile(i) = inputs.bio_sol(ind(i)).value;

                sol_aq_ref(i) = min(iSol .* SF_bile ./iFrac_bile.F_all);
                % Aqueous solubility at pH of bile solubility measurement
                SC_bile  = (sol_aq_ref(i)*1e-3 / mw) / (1/18);
                % SC at pH of in vitro solubility measurement
                % Based on feedback from S+, individual
                % SR = (sol_bile - aq_sol)/(SC_bile * mw * bile_conc)
                SR_num(i) = max(sol_bile(i) - sol_aq_ref(i),0);
                SR_den(i) = SC_bile * mw * bile_conc * 1e-3;
                SR0(i) = (SR_num(i))/(SR_den(i));
            end

            idx = strcmpi(media,'fessif');

            if sum(SR0) == 0
                SR = 0;
            elseif any(SR0==0)
                SR = SR0(SR0~=0);
            elseif sol_bile(idx) - sol_bile(~idx) < 0
                SR = mean(SR0);
            elseif (SR_num(idx) - SR_num(~idx))/(SR_den(idx) - ...
                    SR_den(~idx)) < 0
                SR = mean(SR0);
            else
                SR = (SR_num(idx) - SR_num(~idx))/(SR_den(idx) - ...
                    SR_den(~idx));
            end

            sol = (sol_aq*1e3 + SC * SR * mw * cmpt.Bile)/1e3;
            % Solubility with bile salt effect.  Model output is in ug/mL
            % - the quantity inside () is divided by 1e3 to give mg/mL

        end
    case 'userdefined'
        SR = inputs.SR;
        sol = (sol_aq*1e3 + SC * SR * mw * cmpt.Bile)/1e3;
        if isnan(sol)
            sol = sol_aq;
        end
    otherwise
        error(['Unrecognized solubilization ratio method "%s". ' ...
            'Please choose from "in vitro" or "theoretical."'], ...
            inputs.SR_method);
end

if isempty(sol), sol = 0; end

switch lower(strrep(strrep(inputs.Deff_method,'_',''),' ',''))
    case 'adjusted'
        f_mono = min(sol_aq./sol,1);
        if cmpt.Bile>=3
            Dagg = (0.0015*cmpt.Bile) + 0.081;
            % Linear regression coefficients from G+ help files.
            % Dagg output is in units x10^(-5) cm2/s (multiply Dagg
            % by 1e-5 for proper units in subsequent equation).
            Deff = (inputs.Deff * f_mono) + (Dagg*1e-5*(1-f_mono));
        else
            Deff = (inputs.Deff * f_mono) + (0.012e-5*(1-f_mono));
        end
    otherwise
        Deff = inputs.Deff;
end
end

function [P_para,ASF_para] = calculateParacellular_ka(cmpt,inputs)
%
% Refer to GastroPlus manual and the reference below for description of the
% calculations.
% Adson, A., et al. (1994). "Quantitative approaches to delineate
% paracellular diffusion in cultured epithelial cell monolayers." J
% ournal of Pharmaceutical Sciences 83(11): 1529-1536.
%

r_he = inputs.r_he;
r_s = inputs.r_s;
mw = inputs.mw;
poros = cmpt.Porosity;   %Porosity (fraction pores / pore length)
Deff = inputs.Deff;
R = cmpt.PoreRadius;
kB = 1.38e-16;   %Boltzman constant (erg/K)
T = 310;         %Temperature in K
e = 4.8e-10;     %Unit charge of an ion (units esu)
EPG = 3.4e-4;    %electrical potential gradient (EPG),
% 102 mV = 3.4e-4 statvolt
iFrac = smtb.pbpk.ionfrac(inputs.Acidic_pKa,inputs.Basic_pKa,cmpt.pH);
% Ionization calculations

switch lower(inputs.para_model)
    case {'zhimin','off'} % if "off" the calculation is still
        % run but not applied in parameterizePBPKvariants
        F_r = min(max(((1 - (r_s/R)).^2).*(1-(2.104.*(r_he/R))+ ...
            (2.09.*((r_he/R).^3))-(0.95.*((r_he/R).^5))),0),1);
        % Renkin function (equation 4-32 in G+ manual).
        % Limit is between 0 and 1 according to Adson reference.
        % min(max(F_r,0),1) to prevent negatives and values above 1
    case 'adson'
        r_eff = 0.8 + (0.2207*sqrt(mw)); % equation 4-33 in G+ v9.6 manual
        F_r = min(max(((1 - (r_eff/R)).^2).*(1-(2.104.*(r_eff/R))+ ...
            (2.09.*((r_eff/R).^3))-(0.95.*((r_eff/R).^5))),0),1);
        % equation 4-31 in G+ v9.6 manual
end

P_para = zeros(1,length(iFrac.F_all));
for i = 1:length(iFrac.F_all)   % Multiply by fraction of each species
    % to calculate the net amount from each species.
    % Note this should continue to be checked.
    if iFrac.NetCharge(i)==0    % This accounts for neutral and
        % zwitterions.  May need to separate neutral and zwitterions.
        P_para(i) = iFrac.F_all(i) * poros * Deff * F_r;
    else
        kappa = (e.*iFrac.NetCharge(i).*EPG)/(kB.*T);
        P_para(i) = iFrac.F_all(i) * poros * Deff * F_r *...
            (kappa/(1-exp(-1*kappa)));
    end
end

P_para = sum(P_para);

switch cmpt.Compartment
    case 'Stomach'
        ASF_para = 0;
    otherwise
        ASF_para = 2./cmpt.Radius;
end
end

function k_bas = calculate_kbas(cmpt,inputs,sp,p_para,asf,kbas_coeffs)

fun_kbas = kbas_coeffs{1};
coeffs = kbas_coeffs{2};

peff = inputs.Peff * 10^4;
p_para = p_para * 10^4;
if ~any(strcmpi({coeffs.species},sp))
    % remove after adding rest of species for kbas fitting
    sp = 'Human';
end

if (peff - p_para) * cmpt.Volume * asf > 0
    b = coeffs(strcmpi({coeffs.species},sp)).beta_Int;
    k_bas = fun_kbas(b,log10((peff - p_para) * cmpt.Volume * asf));
    % Add lines below back after Stomach coefficients are estimated.
    %     switch cmpt.Compartment
    %         case 'Stomach'
    %             b = coeffs(strcmpi({coeffs.species},sp)).beta_Stomach;
    %             k_bas = fun_kbas(b,log10((peff - p_para) * cmpt.Volume * ...
    %                 asf));
    %         otherwise
    %             b = coeffs(strcmpi({coeffs.species},sp)).beta_Int;
    %             k_bas = fun_kbas(b,log10((peff - p_para) * cmpt.Volume * ...
    %                 asf));
    %     end
else
    k_bas = 3.7;
end
end