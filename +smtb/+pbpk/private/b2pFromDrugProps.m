function sp_input= b2pFromDrugProps(sp_input, species, dp, ii, ...
    calc_operation)
if contains(lower(strrep(strrep(sp_input.B2P_source,' ',''),'_','')),...
        'kate')
    sp_input.B2P_source = 'KATE';
elseif contains(lower(strrep(strrep(sp_input.B2P_source,' ',''),'_','')), ...
        'qsaradmet')
    sp_input.B2P_source = 'QSAR ADMET';
end

switch lower(strrep(strrep(sp_input.B2P_source,' ',''),'_',''))
    case 'default'
        [val,source,sp,~,comments,n] = getPreferredValue(dp,'B2P', ...
            species{ii},sp_input.B2P_source,calc_operation.param);
        sp_input.B2P        = val;
        sp_input.B2P_source_n = n;
        sp_input.B2P_comments = comments;
        switch lower(source)
            case 'kate'
                sp_input.B2P_source = sprintf('%s %s',source, ...
                    [upper(sp(1)) sp(2:end)]);
            otherwise
                sp_input.B2P_source = source;
        end
    case {'admet','qsaradmet'}
        [val,source,~,~,comments,n] = getPreferredValue(dp,'B2P', ...
            species{ii},sp_input.B2P_source,calc_operation.param);
        if ~isnan(val)
            sp_input.B2P        = val;
            sp_input.B2P_source_n = n;
            sp_input.B2P_source = source;
            sp_input.B2P_comments = comments;
        else
            warning backtrace off
            warning(['Unable to find B2P for species %s, source %s. ' ...
                ' Using default value instead.'],species{ii}, ...
                sp_input.B2P_source,calc_operation.param)
            warning backtrace on
            [val,source,sp,~,comments,n] = getPreferredValue(dp, ...
                'B2P',species{ii},'default',calc_operation.param);
            sp_input.B2P        = val;
            sp_input.B2P_source_n = n;
            sp_input.B2P_comments = comments;
            switch lower(source)
                case 'kate'
                    sp_input.B2P_source = sprintf('%s %s',source, ...
                        [upper(sp(1)) sp(2:end)]);
                otherwise
                    sp_input.B2P_source = source;
            end
        end
    case 'kate'
        [val,source,sp,~,comments,n] = getPreferredValue(dp,'B2P', ...
            species{ii},sp_input.B2P_source,calc_operation.param);
        if isnan(val)
            fub = getFub(dp,'Species',species{ii},'Source','KATE');
            fub_val = abs(exp(mean(log([fub.fu]))));
            if ~isnan(fub_val)
                sp_input.B2P         = sp_input.fup/fub_val;
                sp_input.B2P_source  = 'fup/fub';
                warning backtrace off
                warning(['Unable to find B2P for species %s, source KATE.' ...
                    '  Using BP = fub/fup with measured fub' ...
                    ' and %s fup.'],species{ii},sp_input.fup_source)
                warning backtrace off
            else
                warning backtrace off
                warning(['Unable to find B2P for species %s, source %s.' ...
                    '  Using default value instead.'], ...
                    species{ii},sp_input.B2P_source)
                warning backtrace on
                [val,source] = getPreferredValue(dp,'B2P',species{ii}, ...
                    'default',calc_operation.param);
                sp_input.B2P        = val;
                sp_input.B2P_source = source;
            end
        else
            sp_input.B2P        = val;
            sp_input.B2P_source_n = n;
            sp_input.B2P_comments = comments;
            switch lower(source)
                case 'kate'
                    sp_input.B2P_source = sprintf('%s %s',source, ...
                        [upper(sp(1)) lower(sp(2:end))]);
                otherwise
                    sp_input.B2P_source = source;
            end
        end
    case 'fup/fub'
        % first check for fub values. if found, define fub here
        if ~isempty(sp_input.fub)
            sp_input.B2P        = sp_input.fup/sp_input.fub;
            sp_input.B2P_source = 'fup/fub';
        else
            [val,source,sp,~,comments,n] = getPreferredValue(dp, ...
                'fub',species{ii}, ...
                sp_input.fub_source,calc_operation.param);
            if ~isnan(val) && ~isempty(val)
                sp_input.fub            = val;
                sp_input.fub_source_n   = n;
                sp_input.fub_comments   = comments;
                sp_input.B2P            = sp_input.fup/sp_input.fub;
                sp_input.B2P_source     = 'fup/fub';
                if contains(source,'QSAR','IgnoreCase',true)
                    sp_input.fub_source = source;
                else
                    sp_input.fub_source = sprintf('%s %s', ...
                        source,[upper(sp(1)) lower(sp(2:end))]);
                end
            else
                warning backtrace off
                warning(['Unable to find fub for species %s' ...
                    ' to derive BP = fup/fub.  Using default value ' ...
                    'instead.'],species{ii})
                warning backtrace on
                [val,source,sp,~,comments,n] = getPreferredValue(dp, ...
                    'B2P',species{ii},'default',calc_operation.param);
                sp_input.B2P        = val;
                sp_input.B2P_source_n = n;
                sp_input.B2P_comments = comments;
                switch lower(source)
                    case 'kate'
                        sp_input.B2P_source = sprintf('%s %s', ...
                            source,[upper(sp(1)) lower(sp(2:end))]);
                    otherwise
                        sp_input.B2P_source = source;
                end
            end
        end
    case 'userdefined' % do nothing
    otherwise, error(['Unkown B2P source, please choose from ADMET, ' ...
            'KATE, or fup/fub'])
end

% Define fub if not defined in B2P step
if isempty(sp_input.fub)
    switch lower(strrep(strrep(sp_input.fub_source,' ',''), ...
            '_',''))
        case {'default','kate'}
            [val,source,sp,~,comments,n] = getPreferredValue(dp, ...
                'fub',species{ii},sp_input.fub_source, ...
                calc_operation.param);
            if ~isnan(val) && ~isempty(val)
                sp_input.fub            = val;
                sp_input.fub_source_n   = n;
                sp_input.fub_comments   = comments;
                if contains(source,'kate','IgnoreCase',true)
                    sp_input.fub_source = sprintf('%s %s',source, ...
                        [upper(sp(1)) lower(sp(2:end))]);
                else
                    sp_input.fub_source = source;
                end
            else
                warning backtrace off
                warning(['Unable to find fub for species %s.' ...
                    ' Switching to fup/B2P.'],species{ii})
                warning backtrace on
                sp_input.fub        = sp_input.fup/sp_input.B2P;
                sp_input.fub_source = 'fup/B2P';
            end
        case {'fup/bp','fup/b2p'}
            sp_input.fub        = sp_input.fup/sp_input.B2P;
            sp_input.fub_source = 'fup/B2P';
    end
end

end