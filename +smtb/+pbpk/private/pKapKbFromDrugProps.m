function inputs = pKapKbFromDrugProps(inputs, dp)
% Allowed pkasource -- "QSAR Chemaxon" "QSAR ADMET" or 
% "KATE Spectroscopic pKa"
if strcmpi(inputs.Acidic_pKa_source,'default')
    [val,source] = getPreferredValue(dp,'pKa_a');
    inputs.Acidic_pKa  = val;
    inputs.Acidic_pKa_source  = source;
elseif strcmpi(inputs.Acidic_pKa_source,'User Defined') % do nothing
else
    try
        [val,source] = getPreferredValue(dp,'pKa_a',...
            inputs.Acidic_pKa_source);
        inputs.Acidic_pKa  = val;
        inputs.Acidic_pKa_source  = source;
    catch
        warning(['%s not a valid option for pKa source - valid ' ...
            'options are "QSAR Chemaxon" "QSAR ADMET" or ' ...
            '"KATE Spectroscopic pKa".  Using default value.'])
        [val,source] = getPreferredValue(dp,'pKa_a');
        inputs.Acidic_pKa  = val;
        inputs.Acidic_pKa_source  = source;
    end
end

% Assign Basic pKa
if strcmpi(inputs.Basic_pKa_source,'default')
    [val,source] = getPreferredValue(dp,'pKa_b');
    inputs.Basic_pKa  = val;
    inputs.Basic_pKa_source = source;
elseif strcmpi(inputs.Basic_pKa_source,'User Defined') % do nothing
else
    try
        [val,source] = getPreferredValue(dp,'pKa_a', ...
            inputs.Basic_pKa_source);
        inputs.Acidic_pKa  = val;
        inputs.Basic_pKa_source  = source;
        [val,source] = getPreferredValue(dp,'pKa_b', ...
            inputs.Basic_pKa_source);
        inputs.Basic_pKa  = val;
        inputs.Basic_pKa_source = source;
    catch
        warning(['%s not a valid option for pKa source - ' ...
            'valid options are "QSAR Chemaxon" "QSAR ADMET" or' ...
            ' "KATE Spectroscopic pKa".  Using default value.'])
        [val,source] = getPreferredValue(dp,'pKa_b');
        inputs.Basic_pKa  = val;
        inputs.Basic_pKa_source = source;
    end
end
end