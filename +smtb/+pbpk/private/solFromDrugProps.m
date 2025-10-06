function inputs = solFromDrugProps(inputs, dp)
switch lower(strrep(strrep(inputs.sol_source,' ',''),'_',''))
    case 'default'
        try
            [sol_in,pH_in,sol_source] = getSol(dp,'source', ...
                'KATE CAD','matrix','Aqueous');
            inputs.sol_mgmL = mean(sol_in, 'omitnan');
            inputs.sol_pH = mean(pH_in, 'omitnan');
            inputs.sol_source = sol_source;
            inputs.sol_matrix = 'Aqueous';
            if isnan(inputs.sol_mgmL)
                [sol_in,pH_in,sol_source] = getSol(dp,'source', ...
                    'QSAR ADMET','matrix','Intrinsic');
                inputs.sol_mgmL = sol_in;
                inputs.sol_pH = pH_in;
                inputs.sol_source = sol_source;
                inputs.sol_matrix = 'Intrinsic';
            end
        catch
            [sol_in,pH_in,sol_source] = getSol(dp,'source', ...
                'QSAR ADMET','matrix','Intrinsic');
            inputs.sol_mgmL = sol_in;
            inputs.sol_pH = pH_in;
            inputs.sol_source = sol_source;
            inputs.sol_matrix = 'Intrinsic';
        end
    case 'userdefined'
        inputs.sol_matrix = 'User Defined';
        if isempty(inputs.sol_pH)
            inputs.sol_pH = -1;
        end
    case 'qsaradmetaqueous'
        try
            [sol_in,pH_in,sol_source] = getSol(dp,'source', ...
                inputs.sol_source,'matrix',inputs.sol_matrix);
            inputs.sol_mgmL = sol_in;
            inputs.sol_pH = pH_in;
            inputs.sol_source = sol_source;
        catch
            [sol_in,pH_in,sol_source] = getSol(dp,'source', ...
                'QSAR ADMET','matrix','Intrinsic');
            inputs.sol_mgmL = sol_in;
            inputs.sol_pH = pH_in;
            inputs.sol_source = sol_source;
            inputs.sol_matrix = 'Intrinsic';
        end
    case {'qsaradmet','qsaradmetintrinsic'}
        try
            [sol_in,pH_in,sol_source] = getSol(dp,'source', ...
                inputs.sol_source,'matrix',inputs.sol_matrix);
            inputs.sol_mgmL = sol_in;
            inputs.sol_pH = pH_in;
            inputs.sol_source = sol_source;
        catch
            [sol_in,pH_in,sol_source] = getSol(dp,'source', ...
                'QSAR ADMET','matrix','Intrinsic');
            inputs.sol_mgmL = sol_in;
            inputs.sol_pH = pH_in;
            inputs.sol_source = sol_source;
            inputs.sol_matrix = 'Intrinsic';
        end
    case 'katecad'
        try
            [sol_in,pH_in,sol_source] = getSol(dp,'source', ...
                'KATE CAD','matrix','Aqueous');
            inputs.sol_mgmL = mean(sol_in, 'omitnan');
            inputs.sol_pH = mean(pH_in, 'omitnan');
            inputs.sol_source = sol_source;
            if isnan(inputs.sol_mgmL)
                [sol_in,pH_in,sol_source] = getSol(dp,'source', ...
                    'QSAR ADMET','matrix','Intrinsic');
                inputs.sol_mgmL = sol_in;
                inputs.sol_pH = pH_in;
                inputs.sol_source = sol_source;
                inputs.sol_matrix = 'Intrinsic';
            end
        catch
            [sol_in,pH_in,sol_source] = getSol(dp,'source', ...
                'QSAR ADMET','matrix','Intrinsic');
            inputs.sol_mgmL = sol_in;
            inputs.sol_pH = pH_in;
            inputs.sol_source = sol_source;
            inputs.sol_matrix = 'Intrinsic';
        end
    case 'kateclnd'
        try
            if isempty(dp.sol(strcmpi({dp.sol.source},'KATE CLND')))
                [sol_in,pH_in,sol_source] = getSol(dp,'source', ...
                    'KATE CAD','matrix','Aqueous');
            else
                [sol_in,pH_in,sol_source] = getSol(dp,'source', ...
                    'KATE CLND','matrix','Aqueous');
            end
            inputs.sol_mgmL = mean(sol_in, 'omitnan');
            inputs.sol_pH = mean(pH_in, 'omitnan');
            inputs.sol_source = sol_source;
            inputs.sol_matrix = 'Aqueous';
            if isnan(inputs.sol_mgmL)
                [sol_in,pH_in,sol_source] = getSol(dp,'source', ...
                    'QSAR ADMET','matrix','Intrinsic');
                inputs.sol_mgmL = sol_in;
                inputs.sol_pH = pH_in;
                inputs.sol_source = sol_source;
                inputs.sol_matrix = 'Intrinsic';
            end
        catch
            [sol_in,pH_in,sol_source] = getSol(dp,'source', ...
                'QSAR ADMET','matrix','Intrinsic');
            inputs.sol_mgmL = sol_in;
            inputs.sol_pH = pH_in;
            inputs.sol_source = sol_source;
            inputs.sol_matrix = 'Intrinsic';
        end
    otherwise
        warning(['%s solubility is currently not available. ' ...
            ' Switching solubility to intrinsic from ADMET Predictor'], ...
            inputs.sol_source)
        [sol_in,pH_in,sol_source] = getSol(dp,'source','QSAR ADMET', ...
            'matrix','Intrinsic');
        inputs.sol_mgmL = sol_in;
        inputs.sol_pH = pH_in;
        inputs.sol_source = sol_source;
        inputs.sol_matrix = 'Intrinsic';
end
end