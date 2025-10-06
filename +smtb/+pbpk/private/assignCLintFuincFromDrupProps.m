function sp_input = assignCLintFuincFromDrupProps(sp_input, species, dp, ii)

switch lower(strrep(strrep(sp_input.fuinc_source,' ',''),'_',''))
    case 'kate'
        fuinc = [];
        sys = '';
        source = '';
        spec = '';
        if contains(sp_input.CLint_source,'microsome','IgnoreCase',true)
            [fuinc,source,spec,~,~,~,sys] = getPreferredValue(dp, ...
                'fuinc_microsomes',species{ii}, ...
                'KATE',calc_operation.param);
        elseif contains(sp_input.CLint_source,'hepatocyte', ...
                'IgnoreCase',true)
            [fuinc,source,spec,~,~,~,sys] = getPreferredValue(dp, ...
                'fuinc_hepatocytes',species{ii},'KATE',calc_operation.param);
        end
        if isempty(fuinc) || isnan(fuinc)
            warning backtrace off
            warning(['Unable to find fuinc in KATE for %s,' ...
                ' CLint_fuinc will be set to 1'],species{ii});
            %Consider changing to a default method if KATE data is unavailable.
            warning backtrace on
            sp_input.CLint_fuinc = 1;
            sp_input.fuinc_source = 'default';
        else
            sp_input.CLint_fuinc = fuinc;
            sp_input.fuinc_source = [upper(sys(1)) lower(sys(2:end)) ...
                ' ' upper(spec(1)) lower(spec(2:end)) ' ' '(' source ')'];
        end
    case 'fup'
        sp_input.CLint_fuinc = sp_input.fup;
        sp_input.fuinc_source = sprintf('fup %s',sp_input.fup_source);
    case 'austinsimple'
        if contains(sp_input.CLint_source,'hepatocytes','IgnoreCase', ...
                true) || contains(sp_input.CLint_source,'hepatocyte', ...
                'IgnoreCase',true)
            warning backtrace off
            warning(['Austin simple method only valid for microsomes. ' ...
                'CLint_fuinc %s calculated using Austin hepatocytes'], ...
                species{ii});
            warning backtrace on
            sp_input.CLint_fuinc = calcFuinc(inputs.logD,inputs.logD_pH, ...
                inputs.Acidic_pKa,inputs.Basic_pKa,'Austin Hepatocytes');
            sp_input.fuinc_source = 'Austin Hepatocytes';
        else
            sp_input.CLint_fuinc = calcFuinc(inputs.logD,inputs.logD_pH, ...
                inputs.Acidic_pKa,inputs.Basic_pKa,sp_input.fuinc_source);
        end
    case 'austinhepatocytes'
        sp_input.CLint_fuinc = calcFuinc(inputs.logD,inputs.logD_pH, ...
            inputs.Acidic_pKa,inputs.Basic_pKa,sp_input.fuinc_source);
        sp_input.fuinc_source = [sp_input.fuinc_source ' ' ...
            sscanf(lower(sp_input.CLint_source),' (KATE)')];
    case 'austin'
        if contains(sp_input.CLint_source,'hepatocytes', ...
                'IgnoreCase',true) || contains(sp_input.CLint_source, ...
                'hepatocyte','IgnoreCase',true)
            sp_input.CLint_fuinc = calcFuinc(inputs.logD,inputs.logD_pH, ...
                inputs.Acidic_pKa,inputs.Basic_pKa,'Austin Hepatocytes');
            sp_input.fuinc_source = 'Austin Hepatocytes';
        elseif contains(sp_input.CLint_source,'microsomes', ...
                'IgnoreCase',true) || contains(sp_input.CLint_source, ...
                'microsome','IgnoreCase',true)
            sp_input.CLint_fuinc = calcFuinc(inputs.logD,inputs.logD_pH, ...
                inputs.Acidic_pKa,inputs.Basic_pKa,'Austin Microsomes');
            sp_input.fuinc_source = 'Austin Microsomes';
        end
    case {'hallifax','hallifaxmicrosomes'}
        if contains(sp_input.CLint_source,'hepatocytes','IgnoreCase', ...
                true) || contains(sp_input.CLint_source,'hepatocyte', ...
                'IgnoreCase',true)
            warning backtrace off
            warning(['Hallifax method only valid for microsomes.  ' ...
                'CLint_fuinc %s will be calculated using Austin method'], ...
                species{ii})  %can change this to set fuinc=1
            warning backtrace off
            sp_input.CLint_fuinc = calcFuinc(inputs.logD,inputs.logD_pH, ...
                inputs.Acidic_pKa,inputs.Basic_pKa,'Austin Hepatocytes');
            sp_input.fuinc_source = 'Austin Hepatocytes';
        else
            sp_input.CLint_fuinc = calcFuinc(inputs.logD,inputs.logD_pH, ...
                inputs.Acidic_pKa,inputs.Basic_pKa,sp_input.fuinc_source);
           sp_input.fuinc_source = [sp_input.fuinc_source  ...
               sscanf(lower(sp_input.CLint_source),'%s (KATE)')];
        end
    case {'default','userdefined'}
        %do nothing
    otherwise
        warning('Invalid fuinc method.  CLint_fuinc will be set to 1.')
        sp_input.fuinc_source = 'default';
end
end