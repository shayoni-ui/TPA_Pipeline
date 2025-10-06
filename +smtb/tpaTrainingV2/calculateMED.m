function res_med = calculateMED(vcmpd, results, dose_vec, ...
    pkpd_vars, resp_var, thresh, thresh_mod)
arguments
    vcmpd; % virtual compounds structure array;
    results; % results of the runVirtualScreen;
    dose_vec; % vector of dose amounts
    pkpd_vars = {'pd_TE', 'pd_P', 'U',...
        'UnboundAUCinf_nghrmL','UnboundAUCt_nghrmL',...
        'UnboundCss_ngmL','UnboundClast_ngmL',...
        'UnboundCmax_ngmL'}; 
    resp_var = 'pd.P';
    thresh = 0.5;
    thresh_mod = '<';
    % PBPK-PD model variables to calcuate at MED
end
% CALCULATE MED allows to calculate Minimum effective dose resulting in the
% output quantity of interest
% Please look at SMARCA2 and TRPM2 TPA where MED calculationg was first
% used
v_cmpds = vcmpd;
[~,idx] = intersect(double(string({v_cmpds.vc_id})), ...
    double(string({results.vc_id})),'stable');
res_med = v_cmpds(idx);
clearvars v_cmpds


h = waitbar(0,'building MED table');

for i = 1:length(res_med)
    idx = double(string({results.vc_id})) ==...
        double(string({res_med(i).vc_id}));
    respi = results(idx).(resp_var);
    calc_med = 1;
    med = max(dose_vec);
    med_mod = '>';
    switch thresh_mod
        case '>'
            if respi(dose_vec==max(dose_vec)) < thresh
                med = max(dose_vec);
                med_mod = '>';
                calc_med = 0;
            elseif respi(dose_vec==min(dose_vec)) > thresh
                med = min(dose_vec);
                med_mod = '<';
                calc_med = 0;
            end
        case '<'
            if respi(dose_vec==max(dose_vec)) > thresh
                med = max(dose_vec);
                med_mod = '>';
                calc_med = 0;
            elseif respi(dose_vec==min(dose_vec)) < thresh
                med = min(dose_vec);
                med_mod = '<';
                calc_med = 0;
            end
    end
    if calc_med
        med = interp1(respi,dose_vec,thresh,'pchip');
        med_mod = '=';
    end
    res_med(i).MED_mod = med_mod;
    res_med(i).MED_mg = med;

    for j = 1:length(pkpd_vars)
        try
            res_med(i).(pkpd_vars{j}) = interp1(dose_vec, ...
                results(idx).(pkpd_vars{j}), ...
                res_med(i).MED_mg,'pchip');
        catch
            res_med(i).(pkpd_vars{j}) = NaN;
        end
    end
    waitbar(i/length(res_med));
end

close(h);

res_med = struct2table(res_med);


end