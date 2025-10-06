function fut = calculateFut(tis,fub, Vdss)
arguments
    tis; % structure with information of tissue volumes
    fub; % fraction unbound blood
    Vdss; % Volume of distribution in mL
end
%load('tis.mat');
tisVolumes = [tis.volume];


% Volume of all tissues except blood
Vt = sum(tisVolumes(~(strcmp(string({tis.name}), "Plasma") |...
    strcmp(string({tis.name}), "Blood cells"))));

% Volume of blood
Vb = sum(tisVolumes((strcmp(string({tis.name}), "Plasma") |...
    strcmp(string({tis.name}), "Blood cells"))));


fut = fub * Vt / (Vdss - Vb);

% Alternatively fut = fup/Kpavg
% fut2(ii) = log10(vcmpd(ii).fup ...
%             / calculate_wt_average_kp([vcmpd(ii).pkpk_inputs.Kp], tis));
%  fut(ii) = log10(calculatefut(vcmpd(ii).fub,...
%             vcmpd(ii).Vss_bl_Lkg*1000*70));

end