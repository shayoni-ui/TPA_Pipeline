function [res,dp] = getCompoundKP(gcn,species,nva)
arguments
    gcn = ''; % GSK compound number , string or DrugProps object
    species char = 'human'; % species type
    nva.Method char = '';
    nva.B2P double = [];
    nva.fup double = [];
    nva.logP double = [];
    nva.pKa_a double = [];
    nva.pKa_b double = [];
    nva.ion_type char = '';
    nva.fup_type char = '';
    nva.fut_type char = '';
    nva.alb_source char = '';
    nva.KaBC double = [];
    nva.Kpscale double = 1;
end
%GETCOMPOUNDKP calculates tissue:plasma ratios for a given compound and
%species
%    [res] = getCompoundKP(gcn,species,field,value) - get the Kp's for the
%    compound by pulling drug properties from the webservice
%
%    [res] = getCompoundKP(drugprops,species,field,value) - get Kp's
%    for the compound by giving the script an existing DrugProperty object.
%    This is faster if the properties have already been pulled from the
%    webservice.
% 
%    Fields:
%       gcn = ''; % GSK compound number , string or DrugProps object
%       species char = 'human'; % species type
%       nva.Method char = '';
%       nva.B2P double = [];
%       nva.fup double = [];
%       nva.logP double = [];
%       nva.pKa_a double = [];
%       nva.pKa_b double = [];
%       nva.ion_type char = '';
%       nva.fup_type char = '';
%       nva.fut_type char = '';
%       nva.alb_source char = '';
%       nva.KaBC double = [];   
% %
%   Valid Kp methods are:
%       'Lukacova' (default)
%       'Rodgers Rowland'
%       'Poulin homogeneous'
%       'Poulin Extracellular'
%       'Berezhkovskiy'
%
%   Valid albumin ratio sources are:
%       'Rowland' (Default for Lukacova and Rodger Rowland methods)
%       'Rothschild'
%       'Constant' (Default for Poulin and Berezkovskiy methods)

tissues = {'Adipose';'Brain';'Heart';'Kidney';'Liver';'Lung';...
    'Muscle';'RedMarrow';'ReproOrg';'RestOfBody';...
    'Skin';'Spleen';'YellowMarrow'}; % Dummy placeholder for tissues, 
% not used in calcualtions, but retuned for empty kp structure 

if nargin == 1 % Return dummy structure for given placeholder
    gcn     = '';
    cmpd    = '';
    dp      = smtb.DrugProps.empty;
else
    if isa(gcn,'smtb.DrugProps') && isempty(gcn)   
        dp = gcn;                                       
        cmpd = '';
        gcn = '';
    elseif isa(gcn,'smtb.DrugProps')              
        dp   = gcn;                                  
        cmpd = dp.name;                                      
        gcn = dp.gcn;
    else                              
        dp   = smtb.DrugProps(gcn);                                        
        cmpd = gcn;                  
        gcn = dp.gcn;
    end
%     gcn = dp.gcn;
end

res = struct('Compound',    cmpd,...
             'gcn',         gcn,...
             'Species',     species,...
             'Kp',          struct('Tissue',tissues,'Kp',1, ...
                'Kpu',1,'fut',1,'fut_ic',1,'pctEW',0,'pctIW',0, ...
                'pctNL',0,'pctAP',0,'pctALB',0),...
             'Vss_L',       [],...
             'Vss_Lkg',     [],...
             'Method',      'Lukacova',...
             'fup_type',    'adjusted',...
             'fut_type',    'S+v9.5',...
             'ion_type',    'default',...
             'alb_source',  'default',...
             'Details',     initializeInputs(dp,species),...
             'DrugProps',   dp);  

if ~isempty(nva.B2P)
    res.Details.B2P = nva.B2P; 
    res.Details.B2P_source = 'User Input';
end
if ~isempty(nva.fup)
    res.Details.fup = nva.fup; 
    res.Details.fup_source = 'User Defined';
end
if ~isempty(nva.logP)
    res.Details.logP = nva.logP;
    res.Details.logP_soure = 'User Defined';
end
if ~isempty(nva.pKa_a)
    res.Details.pKa_a = nva.pKa_a;
    res.Details.pKa_source = 'User Defined';
end
if ~isempty(nva.pKa_b)
    res.Details.pKa_b = nva.pKa_b; 
    res.Details.pKa_source = 'User Defined';
end

if ~isempty(nva.Method); res.Method = nva.Method; end
if ~isempty(nva.ion_type); res.ion_type = nva.ion_type; end
if ~isempty(nva.fup_type); res.fup_type = nva.fup_type; end
if ~isempty(nva.fut_type); res.fut_type = nva.fut_type; end
if ~isempty(nva.alb_source); res.alb_source = nva.alb_source; end

if ~isempty(nva.KaBC)
    res.Details.KaBC = nva.KaBC; 
    res.Details.KaBC_source = 'User Defined';
end




% Parse naming of methods for legacy purposes
 switch strrep(lower(res.Method),' ','')
     case {'rr', 'rr_gas','rodgersrowland','rodgers-rowland'}  
         res.Method = 'Rodgers Rowland';           
         % Rodger-Rowland Method
     case {'lukacova','combrr_gas'}      
         res.Method = 'Lukacova';                   
         % combined Rodgers Equations in Gastroplus
     case {'pou_h','pouh_gas','poulinhomogeneous','pau_h'}
         res.Method = 'Poulin Homogeneous';        
         % Poulin & Theil - Homogeneous
     case {'pou_ex','poulinextracellular','pau_ex'}
         res.Method = 'Poulin Extracellular';       
         % Poulin & Theil - Extracellular
     case {'bere_gas','berezhkovskiy','berezhkovskiy(gastroplus)'}
         res.Method = 'Berezhkovskiy (Gastroplus)'; 
         % Berezhkovskiy method in gastroplus
     case {'bere_sim','berezhkovskiy(simcyp)'}              
         res.Method = 'Berezhkovskiy (SimCYP)';
         % Berezhkovskiy method in simcyp
     otherwise, error('Unrecognized method: %s',res.Method)
 end

%Ensure that pKa's are in row vectors becasue they will
% be concatentated during the calculations
res.Details.pKa_a = res.Details.pKa_a(:)';
res.Details.pKa_b = res.Details.pKa_b(:)';

% Initialize albumin values based on methodology
persistent ALBUMIN
if isempty(ALBUMIN)
    load('GastroplusPhysiology.mat','albumin')
    ALBUMIN = albumin; % #ok<NODEF>
else
    albumin = ALBUMIN;
end

if strcmpi(res.alb_source,'default')
    switch lower(res.Method)
        case {'rodgers rowland','lukacova'}
            res.alb_source = 'Rowland';
        case {'poulin homogeneous','poulin extracellular',...
                'berezhkovskiy (gastroplus)','berezhkovskiy (simcyp)'}
            res.alb_source = 'Constant';
        otherwise, error(['No default albumin source set' ...
                ' for method "%s"'],res.Method)
    end
end
alb    = [albumin.(res.alb_source)];
lip    = [albumin.RowlandLIP];
[~,~,ic] = intersect({res.Details.tis.name},{albumin.Tissue},'stable');
for i = 1:length(res.Details.tis) 
    res.Details.tis(i).albumin = alb(ic(i)); 
    res.Details.tis(i).lipids  = lip(ic(i)); 
end
 
if nargin == 1, return; end

res.Details.E2P = (res.Details.B2P + ...
    res.Details.Hct -1)/res.Details.Hct; %Erythrocyte to plasma ratio

res = smtb.pbpk.KpCalculation(res);

Kpval = num2cell([res.Kp.Kp] * nva.Kpscale);
[res.Kp.Kp]  = deal(Kpval{:}); % Scale the Kp value to match
% in vitro data
res.Vss_L   = Vss_cal(res);
res.Vss_Lkg = res.Vss_L ./ res.Details.BodyWgt;

end

function details = initializeInputs(dp,species)

if isempty(dp)
    logP  = []; logPsource = 'empty';
    pKa_a = []; pKasource  = 'empty';
    pKa_b = [];
    B2P   = []; B2Psource = 'empty';
    fup   = []; fupsource = 'empty';
else
    [logP,logPsource] = getPreferredValue(dp,'logP');
    [pKa_a,~] = getPreferredValue(dp,'pKa_a');
    [pKa_b,pKasource] = getPreferredValue(dp,'pKa_b');
    [B2P,B2Psource]   = getPreferredValue(dp,'B2P',species);
    [fup,fupsource]   = getPreferredValue(dp,'fup',species);
end



details = struct('logP',        logP,...
                 'logP_vow',    [],...
                 'logD_pH74',   [],...
                 'logP_source', logPsource,...
                 'B2P',         B2P,...
                 'B2P_source',  B2Psource,...
                 'E2P',         [],...
                 'KaBC',        [],...
                 'KpuBC',       [],...
                 'KaBC_source', 'default',...
                 'fup',         fup,...
                 'fup_adj',     [],...
                 'fup_source',  fupsource,...
                 'Hct',         [],...
                 'pKa_a',       pKa_a,...
                 'pKa_b',       pKa_b,...
                 'pKa_source',  pKasource,...
                 'mol_type',    '',...
                 'ionization',  [],...
                 'neu',         [],...
                 'msbz',        [],...
                 'tis',         [],...
                 'BodyWgt',     [],...
                 'kpmat',       []);

%Load physiology for given species and software
persistent PHYS
if isempty(PHYS)
    load('GastroplusPhysiology.mat','kp_phys')
    PHYS = kp_phys;
end

if strcmpi(species,'Beagle'), phys_idx = strcmpi({PHYS.species},'Dog');
else,                         phys_idx = strcmpi({PHYS.species},species);
end
if ~sum(phys_idx)
    error(['Unable to find physiology for ' ...
        'species "%s""'],species);
end

details.tis     = PHYS(phys_idx).tissue;
Hct             = PHYS(phys_idx).hematocrit;
details.Hct     = Hct;
details.BodyWgt = PHYS(phys_idx).weight;
% details.E2P     = (B2P + Hct -1)/Hct; %Erythrocyte to plasma ratio.
end

function  Vss = Vss_cal(res)
% VSS_CALCULATION version 2 (15 Jan 2016) 
% Vt: volume of different tissues (cell array)
% Vp: volume of plasma
% Kp: tissue plasma participation coefficient
% bw: body weight
% EP: erythrocyte:plasma ratio, can be calculated from
% blood:plasma ratio and hemtocrit

    nm = {res.Details.tis.name};
    Kp = [res.Kp.Kp];
    Vt = [res.Details.tis(~strcmp(nm,'Plasma') & ...
        ~strcmp(nm,'Blood cells')).volume];
    Vp = res.Details.tis(strcmp(nm,'Plasma')).volume;
    Ep = res.Details.E2P;
    Ve = res.Details.tis(strcmp(nm,'Blood cells')).volume;
    Kp(isnan(Kp)) = 0 ; % replace the NaN elements in Kp with 0
    %Kp(7) = Kp(7);
    Vss = (Ve*Ep + Vp + sum(Vt.*Kp, 'omitnan'))./1000; 
    % method that used by gastroplus and simcyp
    
end
