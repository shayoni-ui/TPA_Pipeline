function [inputs, allometry, calc_operation]= compilePKpredictor(inputs, species, dp)
calc_operation.param = inputs.calc_opt.param_op;
calc_operation.pk = inputs.calc_opt.pk_data_op;
calc_operation.iv_data = inputs.calc_opt.iv_data;

if isempty(inputs.PKpredictor)
    % pre-allocate memory for PKpredictor inputs. 5 parameters are compound
    % specicific, 2 are species specific (for now, this will likely change).
    pk_predictor_ops = cell(1, ...
        (length(species(~strcmpi(species,'human'))) * 3 * 2) + (5 * 2));

    pk_predictor_params = {'fup','B2P','fub','logP','logD',...
        'logD_pH','Acidic_pKa','Basic_pKa'};

    ops_ct = 0;
    for i = 1:length(pk_predictor_params)
        switch lower(pk_predictor_params{i})
            case {'fup','b2p','fub'}
                for j = 1:length(species)
                    if ~strcmpi(species{j},'human')
                        idx = find(strcmpi( ...
                            {inputs.SpeciesSpecific.Species},species{j}));
                        if ~isempty(idx)
                            ops_ct = ops_ct + 1;
                            if contains(pk_predictor_params{i},'B2P', ...
                                    'IgnoreCase',true)
                                pk_predictor_ops{ops_ct} = sprintf('%s_%s', ...
                                    'bp',species{j});
                                % currently need this to make this work
                                % with PKpredictor
                            else
                                pk_predictor_ops{ops_ct} = sprintf('%s_%s', ...
                                    pk_predictor_params{i},species{j});
                            end
                            ops_ct = ops_ct + 1;
                            if strcmpi(strrep(lower( ...
                                    inputs.SpeciesSpecific(idx).( ...
                                    sprintf('%s_source', ...
                                    pk_predictor_params{i}))),' ',''), ...
                                    'userdefined')
                                pk_predictor_ops{ops_ct} = ...
                                    inputs.SpeciesSpecific(idx).( ...
                                    pk_predictor_params{i});
                            else
                                pk_predictor_ops{ops_ct} =...
                                    inputs.SpeciesSpecific(idx).( ...
                                    sprintf('%s_source', ...
                                    pk_predictor_params{i}));
                            end
                        end
                    end
                end
            case 'logd_ph'
                ops_ct = ops_ct + 1;
                pk_predictor_ops{ops_ct} = pk_predictor_params{i};
                ops_ct = ops_ct + 1;
                pk_predictor_ops{ops_ct} = inputs.(pk_predictor_params{i});
            otherwise
                ops_ct = ops_ct + 1;
                pk_predictor_ops{ops_ct} = pk_predictor_params{i};
                ops_ct = ops_ct + 1;
                if strcmpi(strrep(lower( ...
                        inputs.(sprintf('%s_source', ...
                        pk_predictor_params{i}))),' ',''),'userdefined')
                    pk_predictor_ops{ops_ct} =...
                        inputs.(pk_predictor_params{i});
                else
                    pk_predictor_ops{ops_ct} = inputs.( ...
                        sprintf('%s_source',pk_predictor_params{i}));
                end
        end
    end
    if ~isempty(dp)
        allometry = smtb.PKpredictor(dp,pk_predictor_ops{:}, ...
            'parameter_operation',calc_operation.param, ...
            'pk_data_operation',calc_operation.pk, ...
            'iv_data',calc_operation.iv_data,...
            'phys',inputs.Physiologies);
        inputs.PKpredictor = allometry;
    else
        allometry = [];
    end
else
    allometry = inputs.PKpredictor;
end
end