function sp_input = fupFromDrugProps(sp_input, species, dp, ii, ...
    calc_operation)
if contains(lower(strrep( ...
        strrep(sp_input.fup_source,' ',''),'_','')),'kate')
    sp_input.fup_source = 'KATE';
elseif contains(lower(strrep(strrep(sp_input.fup_source, ...
        ' ',''),'_','')),'qsaradmet')
    sp_input.fup_source = 'QSAR ADMET';
end
% Assign fup
switch lower(strrep(strrep(sp_input.fup_source,' ',''),'_',''))
    case 'default'
        [val,source,sp,~,comments,n] = getPreferredValue(dp,'fup', ...
            species{ii},sp_input.fup_source,calc_operation.param);
        if isnan(val)||val==0
            switch lower(species{ii})
                case {'rat','human'}
                    fup = getFup(dp,'Species',species{ii},'Source', ...
                        'ADMET');
                    fup_val = abs(exp(mean(log([fup.value]))));
                    sp_input.fup        = fup_val;
                    sp_input.fup_source = fup.source;
                    warning backtrace off
                    warning(['fup for species %s, source %s ' ...
                        'returned a value %s.  ' ...
                        'Using ADMET value instead.'],species{ii}, ...
                        sp_input.fup_source,fup_val)
                    warning backtrace on
                otherwise
                    fup = getFup(dp,'Species','human','Source','ADMET');
                    fup_val = abs(exp(mean(log([fup.value]))));
                    sp_input.fup        = fup_val;
                    sp_input.fup_source = fup.source;
                      warning backtrace off
                    warning(['fup for species %s, ' ...
                        'source %s returned a value %s.  ' ...
                        'Using Human ADMET predicted value instead.'], ...
                        species{ii},sp_input.fup_source,fup_val)
                    warning backtrace on
            end
        else
            sp_input.fup        = val;
            switch lower(source)
                case 'kate'
                    sp_input.fup_source = sprintf('%s %s',source, ...
                        [upper(sp(1)) sp(2:end)]);
                otherwise
                    sp_input.fup_source = source;
            end
            sp_input.fup_source_n = n;
            sp_input.fup_comments = comments;
        end
    case {'admet','kate','gskqsar','ppbhuman','ppbrat',...
            'qsarppbhumanv1.3','qsaradmet','qsarppbhuman',...
            'qsarppbrat','qsarppbratv1.3'}
        [val,source,sp,~,comments,n] = getPreferredValue(dp,'fup', ...
            species{ii},sp_input.fup_source,calc_operation.param);
        if isnan(val)||val==0
            switch lower(species{ii})
                case {'rat','human'}
                    fup = getFup(dp,'Species',species{ii},'Source', ...
                        'ADMET');
                    fup_val = abs(exp(mean(log([fup.value]))));
                    sp_input.fup        = fup_val;
                    sp_input.fup_source = fup.source;
                    sp_input.fup_source_n = length(fup);
                    warning backtrace off
                    warning(['fup for species %s, source %s returned a ' ...
                        'value %s.  Using ADMET value instead.'], ...
                        species{ii},sp_input.fup_source,fup_val)
                    warning backtrace on
                otherwise
                    fup = getFup(dp,'Species','human','Source','ADMET');
                    fup_val = abs(exp(mean(log([fup.value]))));
                    sp_input.fup        = fup_val;
                    sp_input.fup_source = fup.source;
                    sp_input.fup_source_n = length(fup);
                    warning backtrace off
                    warning(['fup for species %s, source %s ' ...
                        'returned a value %s.  Using Human ADMET ' ...
                        'predicted value instead.'],species{ii}, ...
                        sp_input.fup_source,fup_val)
                    warning backtrace on
            end
        else
            sp_input.fup        = val;
            switch lower(source)
                case 'kate'
                    sp_input.fup_source = sprintf('%s %s', ...
                        source,[upper(sp(1)) sp(2:end)]);
                otherwise
                    sp_input.fup_source = source;
            end
            sp_input.fup_source_n = n;
            sp_input.fup_comments = comments;
        end
    case 'userdefined' % do nothing
    otherwise, error(['Unkown fup source, please choose ' ...
            'from QSAR PPB Human, QSAR PPB Rat, QSAR ADMET, or KATE'])
end
end