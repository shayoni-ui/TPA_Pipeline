function inputs = logPlogDfromDrugProps(inputs, dp)
% Allowed logP_source options --"KATE Chrom", "KATE CHI", "KATE IAM",
% "Biobyte/Daylight", "Chemaxon", "QSAR Chrom LogP v1 "S+ logP", or
% "Moriguchi"'
%
% Allowed logD_source options -- "KATE Chrom" or "KATE CHI"
%
if strcmp(inputs.logP_source,'default') &&...
        strcmp(inputs.logD_source,'default')
    [val,source] = getPreferredValue(dp,'logP');
    inputs.logP  = val;
    inputs.logP_source = source;
    inputs.logD = smtb.pbpk.logDcalc(inputs.logP,inputs.logD_pH, ...
        inputs.Acidic_pKa,inputs.Basic_pKa);
    inputs.logD_source = 'calculated';
elseif ~strcmpi(inputs.logP_source,'default') &&...
        (strcmpi(inputs.logD_source,'default') ||...
        strcmpi(inputs.logD_source,'calculated'))
    switch lower(strrep(strrep(inputs.logP_source,' ',''),'_',''))
        case 's+logp'
            [inputs.logP,inputs.logP_source] = getPreferredValue(dp, ...
                'logP','S+logP');
        case 'moriguchi'
            [inputs.logP,inputs.logP_source] = getPreferredValue(dp, ...
                'logP','Moriguchi');
        case {'chemaxon','qsarchemaxon'}
            [inputs.logP,inputs.logP_source] = getPreferredValue(dp, ...
                'logP','Chemaxon');
        case {'biobyte/daylight','qsarbiobyte'}
            [inputs.logP,inputs.logP_source] = getPreferredValue(dp, ...
                'logP','QSAR Biobyte/Daylight');
        case {'qsarchromlogp','qsarchromlogpv1','chromlogpv1'}
            [inputs.logP,inputs.logP_source] = getPreferredValue(dp, ...
                'logP','QSAR ChromLogP v1');
        case {'katechrom','katechromlogp','chromlogp'}
            [val,source] = getPreferredValue(dp,'logP','Chrom logP');
            inputs.logP = val; inputs.logP_source = source;
        case {'katechi','katechilogp','chilogp'}
            [val,source] = getPreferredValue(dp,'logP','CHI logP');
            inputs.logP = val; inputs.logP_source = source;
        case {'kateiam','katelogkiam'}
            [val,source] = getPreferredValue(dp,'logP','KATE logK IAM');
            inputs.logP = val; inputs.logP_source = source;
        case {'chinlogp' 'katechinlogp'}
            [val,source] = getPreferredValue(dp,'logP','KATE CHIN logP');
            inputs.logP        = val; inputs.logP_source = source;
        case 'userdefined' %Already specified
        otherwise, error(['Unrecognized logP source "%s", ' ...
                'please choose from "KATE Chrom",' ...
                ' "KATE CHI", "KATE IAM", "Biobyte/Daylight", ' ...
                '"Chemaxon", "QSAR Chrom LogP v1 "S+ logP", or' ...
                ' "Moriguchi"'],inputs.logP_source);
    end

    inputs.logD = smtb.pbpk.logDcalc(inputs.logP,inputs.logD_pH, ...
        inputs.Acidic_pKa, inputs.Basic_pKa);
    inputs.logD_source = 'calculated';
elseif (strcmpi(inputs.logP_source,'default') ||...
        strcmpi(inputs.logP_source,'calculated')) &&...
        ~strcmpi(inputs.logD_source,'default')
    switch lower(strrep(strrep(inputs.logD_source,' ',''),'_',''))
        case {'katechrom','katechromlogd','katechromlogdph7.4'}
            [val,source] = getPreferredValue(dp,'logP', ...
                'KATE Chrom LogD pH 7.4');
            if contains(source,'logp','IgnoreCase',true) ||...
                    contains(source,'QSAR Biobyte/Daylight', ...
                    'IgnoreCase',true)
                inputs.logP        = val;
                inputs.logP_source = source;
                inputs.logD = smtb.pbpk.logDcalc(inputs.logP,inputs.logD_pH, ...
                    inputs.Acidic_pKa,inputs.Basic_pKa);
                inputs.logD_source = 'calculated';
            else
                inputs.logD        = val;
                inputs.logD_source = source;
                inputs.logD_pH     = 7.4;
                inputs.logP = smtb.pbpk.logPcalc(inputs.logD,inputs.logD_pH, ...
                    inputs.Acidic_pKa,inputs.Basic_pKa);
                inputs.logP_source = 'calculated';
            end
        case {'chilogd','katechilogd','katechilogdph7.4'}
            [val,source] = getPreferredValue(dp,'logP', ...
                'KATE CHI LogD pH 7.4');
            if contains(source,'logp','IgnoreCase',true) ||...
                    contains(source,'QSAR Biobyte/Daylight', ...
                    'IgnoreCase',true)
                inputs.logP        = val;
                inputs.logP_source = source;
                inputs.logD = smtb.pbpk.logDcalc(inputs.logP,inputs.logD_pH, ...
                    inputs.Acidic_pKa,inputs.Basic_pKa);
                inputs.logD_source = 'calculated';
            else
                inputs.logD        = val;
                inputs.logD_source = source;
                inputs.logD_pH     = 7.4;
                inputs.logP = smtb.pbpk.logPcalc(inputs.logD,inputs. ...
                    logD_pH,inputs.Acidic_pKa,inputs.Basic_pKa);
                inputs.logP_source = 'calculated';
            end
        case {'qsarchromlogd','qsarchromlogdph7.4'}
            [val,source] = getPreferredValue(dp,'logP', ...
                'QSAR Chrom LogD pH 7.4');
            if contains(source,'logp','IgnoreCase',true) ||...
                    contains(source,'QSAR Biobyte/Daylight', ...
                    'IgnoreCase',true)
                inputs.logP        = val;
                inputs.logP_source = source;
                inputs.logD = smtb.pbpk.logDcalc(inputs.logP,inputs.logD_pH, ...
                    inputs.Acidic_pKa,inputs.Basic_pKa);
                inputs.logD_source = 'calculated';
            else
                inputs.logD        = val;
                inputs.logD_source = source;
                inputs.logD_pH     = 7.4;
                inputs.logP = smtb.pbpk.logPcalc(inputs.logD, ...
                    inputs.logD_pH,inputs.Acidic_pKa,inputs.Basic_pKa);
                inputs.logP_source = 'calculated';
            end
        case 'userdefined' % Already specified
            inputs.logP = smtb.pbpk.logPcalc(inputs.logD,inputs.logD_pH, ...
                inputs.Acidic_pKa,inputs.Basic_pKa);
            inputs.logP_source = 'calculated';
        otherwise, error(['Unrecognized logD source "%s", ' ...
                'please choose from "KATE Chrom" or "KATE CHI"'], ...
                inputs.logD_source);
    end
elseif strcmp(inputs.logD_source,'calculated') ||...
        strcmp(inputs.logP_source,'calculated')
    %Handles the case where an input structure is pre-specified,
    % no action needed
else
    error(['Both logP and logD specified, ' ...
        'fitting of logD vs pH currently not available'])
end
end