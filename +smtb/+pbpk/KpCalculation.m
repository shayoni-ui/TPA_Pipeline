function res = KpCalculation(res)
% KPCALCULATIONS Version 2 (14 Dec. 2015)
% An improved, optimized code to in silico predict 
% tissue-plasma partition coefficients (Kp)
% 1) Incorporated the different condition into single 
% equations (not many elseif statement)
% 2) Adipose and yellow marrow are considered as fatty tissue
% 3) Both Gastroplus and Simcyp physiological parameters can be selected
% 4) The logDvo:w calculation used in Gastroplus method consider even
% the ionized species will disstributed in the vegetable phase
% 5) The output can be either saved in .txt file or print on the screen
% To call the function: res = KpCalculation(p);
% where p is a structure array include input parameters 
% (based on the Adamantis format)
% however the parameters that are required to calculate Kp are:
% pKa, logP, Bp, Hct, fup, Method, software

p = res.Details;

% Parse physiology
PL = p.tis(strcmp({p.tis.name}, ...
    'Plasma'));
BC = p.tis(strcmp({p.tis.name}, ...
    'Blood cells'));
TS = p.tis(~strcmp({p.tis.name},'Plasma') & ...
    ~strcmp({p.tis.name},'Blood cells'));

% Calculation ionizatoin states
p.pKa_a = p.pKa_a(p.pKa_a ~= 0 & ~isnan(p.pKa_a));
p.pKa_b = p.pKa_b(p.pKa_b ~= 0 & ~isnan(p.pKa_b));
ionization = smtb.pbpk.ionfrac(p.pKa_a, p.pKa_b,[7.4 7 7.22]); 
% Calculate inonization at pH 7.4 (plasma), 
% pH 7 (intracellular water), and pH 7.22 (blood cells)

Xd_pl  = ionization(1).neutral+ionization(1).zwitterionic;
Xd_ic  = ionization(2).neutral+ionization(2).zwitterionic;
Xd_bc  = ionization(3).neutral+ionization(3).zwitterionic;

fz_pl = ionization(1).zwitterionic;% Fraction zwitterionic in plasma
fc_pl = ionization(1).cationic;    % Fraction cationic in plasma
fa_pl = ionization(1).anionic;     % Fracion anionic in plasma
fn_pl = ionization(1).neutral;     % Fraction neutral in plasma
fc_ic = ionization(2).cationic;    
% Fraction cationic in intracellular water
% fn_ic = ionization(2).neutral;     REMOVED since unused
% Fraction neutral in intracellular water

% Check if the compound is neutral
switch lower(res.ion_type)
    case 'default',neu = isempty([p.pKa_a p.pKa_b]);
    case 'legacy', neu = f_neu_pl > 0.999  & f_neu_ic > 0.999;
    otherwise,     error(['Invlid ionfractype, please choose ' ...
            '"default" or "legacy".'])
end

% Calculate lipophilicity parameters
P           = 10^p.logP*ones(1,length(TS)); 
% Partition coefficient vector: Pow for normal tissue
% and Pvow for fatty tissue
p.logD_pH74 = smtb.pbpk.logDcalc(p.logP,7.4,p.pKa_a,p.pKa_b);
Dow_pH74    = 10^p.logD_pH74;
Pvow        = 10^(1.115*p.logP - 1.35);   
% partition coefficient of vegetable oil: water for adipose 
% and yellow marrow
Dvow        = 10^smtb.pbpk.logDcalc(log10(Pvow),7.4,p.pKa_a,p.pKa_b); 
% partition coefficient of vegetable oil: water for adipose and 
% yellow marrow

% Calculate adjusted fup
p.fup_adj = CalculateAdjustedFup(Dow_pH74,PL.fNL,PL.fNP,PL.fEW,p.fup);
switch lower(res.fup_type)
    case 'experimental', fup = p.fup;
    case 'adjusted',     fup = p.fup_adj;
    otherwise, error(['Invalid type. Please choose' ...
            ' "experimental" or "adjusted"']);
end

% Get values for tissue albumin and lipid levels
% nan(1,length(TS)); LIP = nan(1,length(TS)); REMOVED since unused
[~,~,idx] = intersect({TS.name},{p.tis.name},'stable');
ALB = [p.tis(idx).albumin];
LIP = [p.tis(idx).lipids];

% Adjust partition coefficients in adipose and yellow marrow
switch res.Method
    case {'Rodgers Rowland','Lukacova'}
        P(strcmp({TS.name},'Adipose') | ...
            strcmp({TS.name},'YellowMarrow')) = Pvow; 
    otherwise
        P(strcmp({TS.name},'Adipose') | ...
            strcmp({TS.name},'YellowMarrow')) = Dvow;
        res.fut_type = 'poulin';
        fut = 1./(1+((1-fup)/fup).*ALB);
end

%% Perform partion calculations

switch res.Details.KaBC_source  
    % Assign KaBC from default calculation or user input
    case 'default'
        p.KpuBC = (p.B2P-1+p.Hct)/(p.Hct*fup); 
        % that need to be used to calculate volume of distribution
        p.KaBC  = (p.KpuBC - (1/Xd_bc)/(1/Xd_pl)*BC.fIW - ...
            (10^p.logP*BC.fNL + (0.3*10^p.logP + 0.7)*BC.fNP)/(1/Xd_pl) ...
            ) * ((1/Xd_pl)/(BC.AP*((1/Xd_bc)-1)));
    otherwise
        p.KaBC = res.Details.KaBC;
        p.KpuBC = (p.KaBC./((1/Xd_pl)/(BC.AP*((1/Xd_bc)-1))) ...
            ) + ((1/Xd_bc)/(1/Xd_pl)*BC.fIW) + ((10^p.logP*BC.fNL + ...
            (0.3*10^p.logP + 0.7)*BC.fNP)/(1/Xd_pl));
        p.B2P_fit = (p.KpuBC*p.Hct*fup) - p.Hct + 1;
end
msbz    = ~isempty(p.pKa_b) && max(p.pKa_b) >= 7; 
% logical indicator of moderate-to-strong base/zwitterions 
% (at least 1 basic pKa > 7)

Kpu   = nan(1,length(TS));
kpmat = zeros(5,length(TS));
RAtp  = (neu*LIP + (1-neu)*ALB);

switch res.Method
    case 'Rodgers Rowland'
        kpmat(1,:) = [TS.fEW];                                                                       
        % Drug in extracellular water
        kpmat(2,:) = (1/Xd_ic)/(1/Xd_pl)*[TS.fIW];                                                    
        % Drug in intracellular water
        kpmat(3,:) = (P.*[TS.fNL] +(0.3*P + 0.7).*[TS.fNP])/(1/Xd_pl);                                
        % Drug bound to neutral lipids and phospholipids
        kpmat(4,:) = msbz * max((p.KaBC*[TS.AP]*(1/Xd_ic-1))/(1/Xd_pl),0);                            
        % Drug bound to Acid  Phospholipids, if KaBC is negative, 
        % set that term is 0
        kpmat(5,:) = (1 - msbz) * (1/fup-1-( ...
            P.*PL.fNL + (0.3.*P+0.7)*PL.fNP)/(1/Xd_pl)).* RAtp;      
        % Drug bound to Lipids or Albumin
        
        Kpu = sum(kpmat,1);
        Kp  = Kpu*fup;
        
    case 'Lukacova'
        kpmat(1,:) = [TS.fEW];                                                                       
        % Drug in extracellular water
        kpmat(2,:) = (1/Xd_ic)/(1/Xd_pl)*[TS.fIW];                                                    
        % Drug in intracellular water
        kpmat(3,:) = (P.*[TS.fNL] +(0.3.*P + 0.7).*[TS.fNP])./(1/Xd_pl);                              
        % Drug bound to neutral lipids and phospholipids
        kpmat(4,:) = (fc_pl + fz_pl) * max(( ...
            p.KaBC*[TS.AP]*(1/Xd_ic-1))/(1/Xd_pl),0);                 
        % Drug bound to Acid  Phospholipids, if KaBC is negative,
        % set that term is 0
        kpmat(5,:) = (fa_pl + fn_pl) * (1/fup - 1 - (P*PL.fNL + ...
            (0.3*P+0.7)*PL.fNP)/(1/Xd_pl)).* RAtp;  
        % Drug bound to Lipids or Albumin
        
        Kpu = sum(kpmat,1);
        Kp  = Kpu*fup;
        
    case 'Poulin Homogeneous'
        Kp  = (P.*([TS.fNL]  + 0.3*[TS.fNP]) + ([TS.Vwt] + ...
            0.7*[TS.fNP]))./(P*(PL.fNL + 0.3*PL.fNP) + ...
            (PL.fEW + 0.7*PL.fNP)).*(fup./fut);
    case 'Poulin Extracellular'
        Kp  = fup./fut; %Calculation for permeability limited tissues, 
        % for perfusion limitied tissues, Kp_ec = fEW * fup/fut
    case 'Berezhkovskiy (Gastroplus)'
        Kp  = (P.*([TS.fNL]  + 0.3*[TS.fNP]) + ...
            (0.7*[TS.fNP] + [TS.Vwt]./fut))./...
              (P*(PL.fNL     + 0.3*PL.fNP)   + ...
              (0.7*PL.fNP)  + PL.Vwt/fup);
    case 'Berezhkovskiy (SimCYP)'
        fut = 1/(1+(1-fup)/(2*fup)).*ones(length(TS),1);
        Kp = (P.*([TS.fNL]  + 0.3*[TS.fNP]) + ...
            ([TS.Vwt]./fut + 0.7*[TS.fNP]))./(P*(PL.fNL + 0.3*PL.fNP) ...
            + (PL.Vwt/fup + 0.7*PL.fNP));
        warning('SimCYP Berezhkovskiy not verified')
    otherwise, error('Unrecognized method: %s',p.method)
end

% Calculate compound fut
switch lower(res.fut_type)
    case 'poulin', fut = 1./(1+((1-fup)/fup).*ALB);
    case {'s+9.0','s+v9.0'}, error(['Simulation plus v9.0 fut ' ...
            'calculations not currently implemented'])
        % fut = (fIW*(1+X)/(1+Y) + fEW)./Kpu; 
        % Old method of Fut calculation - does not match GastroPlus fut
    case {'s+','s+9.5','s+v9.5','s-plus'}
        fut = fup ./ Kp .* Xd_pl./Xd_ic;
    otherwise, error(['Invalid fut type, ' ...
            'please choose "poulin" or "s-plus".']);
end
if isinf(p.KaBC)
    KaBC_fut = -1;
else
    KaBC_fut = p.KaBC;
end
fut_ic = (1-[TS.fEW]) ./ ([TS.fIW] + ...
    fc_pl.*max(KaBC_fut,0).*[TS.AP].*fc_ic + Xd_ic .* (P.*[TS.fNL] ...
    +(0.3.*P+0.7).*[TS.fNP])); 
% Intracellulr fut. if KaBC is negative, set that term to 0

%Collect calculation details
p.molType    = [num2str(length(p.pKa_a)),...
    'acid',num2str(length(p.pKa_b)),'base']; 
% basically tell you how many pKa_a and how many pKa_b you have
p.ionization = ionization;
p.kpmat      = kpmat;
p.neu        = neu;
p.msbz       = msbz;
p.logP_vow   = log10(Pvow);

% Build output structure
for i = 1:length(TS)
    res.Kp(i).Tissue = TS(i).name;
    res.Kp(i).Kp     = Kp(i);
    res.Kp(i).Kpu    = Kpu(i);
    res.Kp(i).P      = P(i);
    res.Kp(i).fut    = fut(i);
    res.Kp(i).fut_ic = fut_ic(i);
    res.Kp(i).pctEW  = kpmat(1,i)./sum(kpmat(1:5,i)) * 100;
    res.Kp(i).pctIW  = kpmat(2,i)./sum(kpmat(1:5,i)) * 100;
    res.Kp(i).pctNL  = kpmat(3,i)./sum(kpmat(1:5,i)) * 100;
    res.Kp(i).pctAP  = kpmat(4,i)./sum(kpmat(1:5,i)) * 100;
    res.Kp(i).pctALB = kpmat(5,i)./sum(kpmat(1:5,i)) * 100;
end

res.Details = p;

end