function [inputs,species] = parseInputsFromDrugProps(inputs, species, nva)
dp = inputs.DrugProps;
%disp(nva);
% Assing acidic and basic using drop props
inputs = pKapKbFromDrugProps(inputs, dp);

% Assign logP and calculate logD
inputs = logPlogDfromDrugProps(inputs, dp);

% Assingn Deff... If drugProps is empty and Deff value is not passed this
% shoud throw an error
if strcmp(inputs.Deff_source,'default')
    inputs.Deff        = dp.Deff.value;
    inputs.Deff_source = dp.Deff.source;
end

% Assign SF... If drugProps is empty and SF value is not passed this
% shoud throw an error
if strcmpi(inputs.SF_source,'default')
    if ~isempty(dp)
        inputs.SF           = dp.sf.value;
        inputs.SF_source    = dp.sf.source;
    else
        inputs.SF           = 20;
        inputs.SF_source    = 'SF empty in DrugProps, assumed SF=20';
    end
end

% Assign mw... If drugProps is empty and MW value is not passed this
% shoud throw an error
if isempty(inputs.mw)
    if ~isempty(inputs.DrugProps)
        if isempty(inputs.DrugProps.mw)
            inputs.mw = 0;
        else
            inputs.mw = inputs.DrugProps.mw;
        end
    else
        inputs.mw = 0;
    end
end

% Assign Deff correction...
switch inputs.Deff_method
    case {'default','adjusted'}, inputs.Deff_method = 'adjusted';
    case 'unadjusted', inputs.Deff_method = 'unadjusted';
    otherwise, error(['Unrecognized Deff method. ' ...
            ' Please choose from "adjusted" or "unadjusted."'])
end

% Assign solubilization ratio (SR) method
switch lower(strrep(strrep(inputs.SR_method,' ',''),'_',''))
    case {'default','invitro'}, inputs.SR_method = 'in vitro';
    case 'theoretical', inputs.SR_method = 'theoretical';
    otherwise, error(['Unrecognized SR method.  ' ...
            'Please choose from "in vitro" or "theoretical."'])
end

% Assign solubility
inputs = solFromDrugProps(inputs, dp);

% Assign FaSSGF solubility
inputs = faSSGFsolFromDrugProps(inputs, dp);

% Assign FaSSIF solubility
inputs = faSSIFsolFromDrugProps(inputs, dp);
% Assign FeSSIF solubility
inputs = feSSIFsolFromDrugProps(inputs, dp);

% Assign paracellular absorption option
switch lower(inputs.para_abs)
    case {'default','zhimin','on'}, inputs.para_abs = 'Zhimin';
    case 'adson',                   inputs.para_abs = 'Adson';
    case 'off',                     inputs.para_abs = 'off';
    otherwise, error(['Unrecognized paracellular absorption option. ' ...
            ' Please choose from "on" or "off."'])
end

%Compile PKpredictor if it is empty
[inputs, allometry, calc_operation] = compilePKpredictor(inputs, ...
    species, dp);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters that are species-specific
for ii = 1:length(species)
    sp_input = inputs.SpeciesSpecific( ...
        strcmpi({inputs.SpeciesSpecific.Species},species{ii}));
    if ~isempty(sp_input)

        % Assing fup
        sp_input = fupFromDrugProps(sp_input, species, dp, ii, ...
            calc_operation);

        % Assign B2P and fub
        sp_input = b2pFromDrugProps(sp_input, species, dp, ii, ...
            calc_operation);



        % Assign dose_vol (dose volume)
        % Need to check and finish this
        if sp_input.dose_vol==0
            switch lower(species{ii})
                case 'human', sp_input.dose_vol = 250;
                case 'rat', sp_input.dose_vol = 1.25;
                case 'mouse', sp_input.dose_vol = 0.25;
            end
        end

        % Calculate and assign Kp's and adjusted fup
        sp_input.Kp_perf = smtb.pbpk.getCompoundKP(dp,species{ii},...
            'fup',     sp_input.fup,...
            'fup_type',inputs.fup_type,...
            'B2P',     sp_input.B2P,...
            'logP',    inputs.logP,...
            'pKa_a',   inputs.Acidic_pKa,...
            'pKa_b',   inputs.Basic_pKa,...
            'method',  inputs.Kp_method_perf,...
            'fut_type',inputs.fut_type, ...
            'Kpscale', nva.Kpscale);

        inputs.logD_pH74 = sp_input.Kp_perf.Details.logD_pH74;

        sp_input.Kp_perm = smtb.pbpk.getCompoundKP(dp,species{ii},...
            'fup',     sp_input.fup,...
            'fup_type',inputs.fup_type,...
            'B2P',     sp_input.B2P,...
            'logP',    inputs.logP,...
            'pKa_a',   inputs.Acidic_pKa,...
            'pKa_b',   inputs.Basic_pKa,...
            'method',  inputs.Kp_method_perm,...
            'fut_type',inputs.fut_type);

        switch inputs.fup_type
            case 'experimental'
                sp_input.adjusted_fup = sp_input.fup;
            case 'adjusted'
                sp_input.adjusted_fup = sp_input.Kp_perf.Details.fup_adj;
        end
        sp_input.fup_type = inputs.fup_type;

        switch lower(species{ii})
            case 'human' % do nothing
            otherwise % check to ensure no preclinical species
                % only have valid options. Switch to KATE if
                % invalid option was specified
                if ~any(contains({'microsomes','hurel','hepatocytes', ...
                        'kate','calculated from CLint', ...
                        'User Defined (blood)', ...
                        'User Defined (plasma)'}, ...
                        sp_input.CL_source,'IgnoreCase',true))
                    sp_input.CL_source = 'KATE';
                end
        end

        % Estimate CL_tot_pl and CL_tot_bl

        sp_input = estimateCLtotPlCLtotBlFromDrugProps(sp_input, ...
            species, dp, ii,allometry);

        % Calculate CLint. Added abs as negative values have been observed
        % in KATE
        sp_input = calculateCLintFromDrugProps(sp_input, species, dp, ...
            ii,calc_operation);

        %Assign CLint_fuinc
        sp_input = assignCLintFuincFromDrupProps(sp_input, species, dp, ...
            ii); 

        inputs.SpeciesSpecific(strcmpi({inputs.SpeciesSpecific.Species}, ...
            species{ii})) = sp_input;
    end
end

end
