function [resStruct, res_sim, failed, me] = runSingleCmpdSim( ...
    exportedModel, vc_id, protocolTable, species, profiles2Save)
%% RUNSINGLECMPDSIM is a copy of the RUNSIMULATION function in the toolbox
% that would allow to run and save time course of a single virtual drug.
% In future iteration the file wont be needed.
% RUNSIMULATION - This function is essentially a wrapper around
% simulatePBPKPDModel().
%
% INPUTS
% 1) exportedModel - (SimBiology.export.Model) The model to simulate.
% 2) vc_id - (double) ID of the virtual compound. Not used for anything,
%               just added to resStruct. 
% 3) protocolTable - (table) A table specifying various simulation setting
%               options, such as methods to calculate certain ADME
%               properties, dose regimen, etc.
% 4) species - (char) Just passed through to simulatePBPKPDModel().
% 5) profiles2Save (cell array) - allows to save the time course of the 
% variables of interest upon simulation of the integrated PBPK-PD model
%
% OUTPUTS
% 1) resStruct - (struct)
% 2) failed - (boolean) Indicates if run failed.
% 3) me - (MException) Matlab exception object generated if the simulation
%           fails.

%%%%%%%%%%%%%%%%%%%%%%%%% START INPUT CHECKING %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exported Model
if ~exist('exportedModel', 'var') || isempty(exportedModel)
    error('exportedModel is required.');
elseif ~isa(exportedModel, 'SimBiology.export.Model')
    error(['exportedModel must be an exported SimBiology model but ' ...
        'received an object of class "%s".'], class(exportedModel));
end

% vc_id
if ~exist('vc_id', 'var') || isempty(vc_id) || ~ischar(vc_id)
    error('vc_id must be a non-empty char.');
end

% Protocol Table
if ~exist('protocolTable', 'var') || isempty(protocolTable)
    error('protocolTable is required.');
elseif ~istable(protocolTable)
    error(['protocolTable must be a table but received an object ' ...
        'of class "%s".'], class(protocolTable));
else
    vnames = protocolTable.Properties.VariableNames;
    if ~isempty(setdiff(vnames, {'param', 'type', 'setting'}))
        error(['protocolTable must have headers "param", "type", ' ...
            'and "setting".']);
    end
    % TO DO - Add additional checks?
end

% Species
% Need to get complete list of acceptable species.
if ~exist('species', 'var') || isempty(species) || ~ischar(species)
    error('Invalid entry for species.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%% END INPUT CHECKING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%% START BUILD DOSE VECTOR %%%%%%%%%%%%%%%%%%%%%%%%%
pTypes = protocolTable.type;
doseIdx = strcmpi(pTypes, 'dose');

pParams = protocolTable.param; 
minIdx = doseIdx & strcmpi(pParams, 'min_amount');
maxIdx = doseIdx & strcmpi(pParams, 'max_amount');
nDoseIdx = doseIdx & strcmpi(pParams, 'n_dose');
routeIdx = doseIdx & strcmpi(pParams, 'route');
unitIdx = doseIdx & strcmpi(pParams, 'amount_units');
intIdx = doseIdx & strcmpi(pParams, 'interval');
totalDoseIdx = doseIdx & strcmpi(pParams, 'total_doses');

minDoseAmount = protocolTable.setting{minIdx};
maxDoseAmount = protocolTable.setting{maxIdx};
nDoseScan = protocolTable.setting{nDoseIdx};
doseVec = logspace( ...
    log10(minDoseAmount), log10(maxDoseAmount), nDoseScan);

doseRoute = protocolTable.setting{routeIdx};
doseUnit = protocolTable.setting{unitIdx};
doseInterval = protocolTable.setting{intIdx};
totalDoses = protocolTable.setting{totalDoseIdx};
%%%%%%%%%%%%%%%%%%%%%%%%%% END BUILD DOSE VECTOR %%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% START GET SIMULATION SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%
response = protocolTable.setting{strcmpi(pParams, 'response')};
effResponse = strsplit(response, ',');
if ischar(effResponse)
    effResponse = {effResponse};
end
% Remove any extraneous commas.
effResponse = effResponse(~strcmpi(effResponse, ','));
% Generate labels appropriate for struct field names.
effResponseLabel = strrep(effResponse, '.', '_');

% We take the sim time entered by the user and convert it to a list of
% output times. I guess we want to ensure we don't store too much
% information and that the simulations all use the same output times so we
% can directly compare derived quantities like AUC.
sim_time = protocolTable.setting{strcmpi(pParams, 'sim_time')};
if length(sim_time) == 1
    if sim_time < 200
        sim_time = [1e-3, 0.5:0.5:sim_time];
    elseif sim_time < 4000
        sim_time = [1e-3, 1:1:sim_time];
    else
        sim_time = [1e-3, 2:2:sim_time];
    end
end

aggregationTime = protocolTable.setting{strcmpi(pParams, 'time_span')};
aggregationType = protocolTable.setting{ ...
    strcmpi(pParams, 'sim_aggregation')};
%%%%%%%%%%%%%%%%%%%%% END GET SIMULATION SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%


resStruct = struct.empty;
failed = false;
me = '';
    
%%% simulate model %%%
try
    [res_sim, ~, errMsgs] = simulatePBPKPDModel( ...
        exportedModel, species, ...
        'route', doseRoute, ...
        'dose amount', doseVec, ...
        'dose amount unit', doseUnit, ...
        'total doses', totalDoses, ...
        'dose_interval', doseInterval, ...
        'sim_time', sim_time, ...
        'Css_calc', 'lastdose', ...
        'profiles', profiles2Save, ...
        'response', effResponse, ...
        'time_span', aggregationTime, ...
        'valtype', aggregationType);
catch me
    failed = true;
end

% If the simulation didn't run because of issues with the input parameters,
% then res_sim will be an empty struct. We want to elevate this to an error
% because from the user's perspective this is the same as a failed run.
if ~failed && isempty(res_sim)
    try
        if isempty(errMsgs)
            error(['Simulation did not run but errMsgs was empty. ' ...
                'This is unexpected, please contact a developer.']);
        else
            error('Could not run simulation. First error message: %s', ...
                errMsgs{1});
        end
    catch me
        failed = true;
    end
end

%%%%%%%%%%%%%%%%%%%% START STORE SIMULATION RESULTS %%%%%%%%%%%%%%%%%%%%%%%
if ~failed
    resStruct = struct('vc_id', vc_id);
    for i=1:length(effResponseLabel)
        resStruct.(effResponseLabel{i}) = [res_sim.(effResponseLabel{i})];
    end
    
    % modelVars is an Nx3 cell array, where
    % N=length(exportedModel.ValueInfo), that is, it's a list of every
    % compartment, species, and parameter in the model. The 1st column is
    % the qualified name, the 2nd column is the initial value, and the 3rd
    % column contains the units.
    modelVars = getModelValues(exportedModel);
    varNames = modelVars(:, 1);
    % This next line is very specific to Bob's PD notation convention. We
    % might want to let the user specify parameters of interest.
    pdIdx = contains(varNames, 'RC50', 'IgnoreCase', true) ...
        | contains(varNames, 'Kd', 'IgnoreCase', true);
    potencyVars = modelVars(pdIdx, :);

    % Calculate the time above Kd and/or RC50.
    KdIdx = find(contains(potencyVars(:,1), 'Kd', 'IgnoreCase', true));
    rc50Idx = find(contains(potencyVars(:,1), 'RC50', 'IgnoreCase', true));

    % Need fup
    fup = modelVars{strcmpi(modelVars(:, 1), 'fup'), 2};

    % Loop over doses
    for i=1:length(res_sim)
        cpSdata = selectbyname(res_sim(i).sim_data, 'plasma.drug');
        % Get data from the last dose interval.
        cpIdx = cpSdata.Time >= (totalDoses - 1)*doseInterval;
        t = cpSdata.Time(cpIdx);
        cpu = cpSdata.Data(cpIdx).*fup;

        % Loop over each instance of a dissociation constant.
        for j=1:length(KdIdx)
            Kdvalue = potencyVars{KdIdx(j), 2};
            Kdname = potencyVars{KdIdx(j), 1};
            if any(cpu > Kdvalue)
                firstTime = t(find(cpu>Kdvalue, 1, 'first'));
                lastTime = t(find(cpu>Kdvalue, 1, 'last'));
                resStruct.(['t_', Kdname])(i) = lastTime - firstTime;
            else
                resStruct.(['t_' Kdname])(i) = 0;
            end
        end

        % Loop over each instance of an RC50
        for j=1:length(rc50Idx)
            rc50value = potencyVars{rc50Idx(j), 2};
            rc50name = potencyVars{rc50Idx(j), 1};
            if any(cpu > rc50value)
                firstTime = t(find(cpu > rc50value, 1, 'first'));
                lastTime = t(find(cpu > rc50value, 1, 'last'));
                resStruct.(['t_', rc50name])(i) = lastTime - firstTime;
            else
                resStruct.(['t_' rc50name])(i) = 0;
            end
        end

        % Get the fields on res_sim(i).nca_res which are not dose_mg or
        % F_profiles. As far as I can tell those aren't fields on nca_res,
        % so maybe this had to do with some old code. I think there are
        % just 5 fields on nca_res: auc_type, num_points, llq, css_calc,
        % and fup.
        ncaFnames = fieldnames( ...
            rmfield(res_sim(i).nca_res, {'dose_mg', 'F_profiles'}));
        for j=1:length(ncaFnames)
            nName = ncaFnames{j};
            if isempty(res_sim(i).nca_res.(nName))
                resStruct.(nName)(i) = NaN;
            else
                resStruct.(nName)(i) = res_sim(i).nca_res.(nName);
            end
        end
    end

end

%%%%%%%%%%%%%%%%%%%%% END STORE SIMULATION RESULTS %%%%%%%%%%%%%%%%%%%%%%%%
end

