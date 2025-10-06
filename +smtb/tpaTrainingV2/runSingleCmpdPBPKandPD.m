function [virtCpdsStructArray, results, simData, doseAmounts, failedRunsStruct] = ...
    runSingleCmpdPBPKandPD(...
    exportedModel, pdVariantContent, virtualCpdsStruct, ...
    protocolTable, saveFilename, profiles2Save)
%% RUNSINGLECMPDPBPKandPD is minorly modified version of RUNVIRTUALSCREEN 
% to allow running a PBPK-PD model with a single virtual compound and
% extract profiles for debugging purposes. In future iteration there should
% be no need of this file
%
%
% INPUTS
% 1) exportedModel - (SimBiology.export.Model) The PBPK/PD model.
% 2) pdVariantContent - (cell) A 1x2 cell array as created by
%                       buildPDVariantContent(). The first entry is just a
%                       variant name and is ignored. The second entry
%                       contains the Nx1 variant content cell array. See
%                       buildPDVariantContent for additional details.
% 3) virtualCpdsStruct - (struct) A 1x1 struct with fields describing the
%                       virtual compounds. Required fields are Ncompounds
%                       and cpdId. Expected fields are solubility,
%                       solubility_pH, SF, logP, Peff, Deff, BP, fup,
%                       CLint_u, mw, pKa_a, and pKa_b. Other fields
%                       may be allowed (I'll update as I come across them).
%                       The values are either single double values (when
%                       the parameter is constant) or an Nx1 double array
%                       listing the value of that parameter for each of the
%                       N compounds.
% 4) protocolTable - (table) An Nx3 table, with variable names 'param',
%                       'type', and 'setting'. param is one of several
%                       values which describe model components, such as
%                       'target_tissue', 'route', 'fup_method', 'sim_time',
%                       etc. I'll describe the full list as I understand
%                       it. type is one of 'pd', 'dose', 'absorption',
%                       'distribution', 'elimination', or 'simulation
%                       setting'. These help to categorize the parameters.
%                       setting is the setting of the parameter, either a
%                       value (like 14 for total_doses) or a name (like
%                       'S+v9.5' for 'fut_method').
% 5) saveFilename - (char) If non-empty, all output variables are saved to
%                       a mat file named saveFilename.
% 
% 6) profiles2Save - (cell array) - cell array of the names of variables to
% save the entire time course from the simulation
% OUTPUTS

%%%%%%%%%%%%%%%%%%%%%%%%% START INPUT CHECKING %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exported Model
if ~exist('exportedModel', 'var') || isempty(exportedModel)
    error('exportedModel is required.');
elseif ~isa(exportedModel, 'SimBiology.export.Model')
    error(['exportedModel must be an exported SimBiology model but ' ...
        'received an object of class "%s".'], class(exportedModel));
end

% PD Variant Content
if ~exist('pdVariantContent', 'var') || isempty(pdVariantContent)
    error('pdVariantContent is required.');
elseif ~iscell(pdVariantContent) || size(pdVariantContent, 1) ~= 1 ...
        || size(pdVariantContent, 2) ~= 2
    error('pdVariantContent must be a 1x2 cell array.');
elseif ~ischar(pdVariantContent{1})
    error('pdVariantContent{1} should be the name of the variant.');
elseif ~iscell(pdVariantContent{2}) || ...
        (~isempty(pdVariantContent{2}) && size(pdVariantContent{2}, 2) ~= 1)
    error('pdVariantContent{2} should be an Nx1 variant cell array.');
end
for i=1:length(pdVariantContent{2})
    vcRow = pdVariantContent{2}{i};
    if ~iscell(vcRow) || size(vcRow, 1) ~= 1 || size(vcRow, 2) ~= 4
        error('Variant content element %i is not a 1x4 cell array.');
    end
end

% Virtual Compound Struct
EXPECTED_VIRTUAL_CPDS_FIELDS = { ...
    'solubility', 'solubility_pH', 'SF', 'Peff', 'logP', 'Deff', ...
    'pKa_a', 'pKa_b', 'BP', 'fup', 'CLint_u', 'mw', ...
    'Ncompounds', 'cpdId'}';
if ~exist('virtualCpdsStruct', 'var') || isempty(virtualCpdsStruct)
    error('virtualCpdsStruct is required.');
elseif ~isstruct(virtualCpdsStruct)
    error(['virtualCpdsStruct must be a struct but received an ' ...
        'object of class "%s".'], class(virtualCpdsStruct));
else
    fnames = fieldnames(virtualCpdsStruct);
    if ~isempty(setdiff(EXPECTED_VIRTUAL_CPDS_FIELDS, fnames))
        error(['virtualCpdsStruct is missing at least one of the 14 ' ...
            'required fields. See function documentation for details.']);
    end
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
    % Check that total_doses*interval <= sim_time. If this is not true then
    % the code which gets the plasma concentration over the last dose does
    % not work.
    protocolParams = protocolTable.param;
    protocolSettings = protocolTable.setting;
    totalDoses = protocolSettings{strcmp('total_doses', protocolParams)};
    doseInterval = protocolSettings{strcmp('interval', protocolParams)};
    simTime = protocolSettings{strcmp('sim_time', protocolParams)};
    if totalDoses*doseInterval > simTime
        error(['Must have total_doses * interval <= sim_time, but ' ...
            'total_doses*interval = %0.2f and sim_time = %0.2f.'], ...
            totalDoses*doseInterval, simTime);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%% END INPUT CHECKING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Need to define exportedModel as a different variable so it works in
% parfor with applyVariant. (Is that really required?)
m0 = exportedModel;
if ~isAccelerated(m0)
    accelerate(m0);
end
virtualCpdsStruct.Ncompounds = 1; % added line for single compound simulation
Ncpds = virtualCpdsStruct.Ncompounds;

% Only working with human at the moment
species = 'Human';

% Original code loaded all GastroPlus physiologies and passed that data on
% to runVirtualSimulation and simulatePBPKPDModel, but in each case the
% only information being used was the species bodyweight. No point in
% loading all that data just for one number, so I just wrote a sub-function
% which gives that data.
% load('GastroplusPhysiology.mat', 'pbpk');
% phys = pbpk(strcmpi({pbpk.Species},species));
% bw_kg = phys.Properties.Weight;
bw_kg = speciesWeight(species);

% Get the variant content.
pdVariant = pdVariantContent{2};

% Modifying virtualCpdsStruct. Create a new struct with just the "extra"
% fields, those which are not one of the expected fields. This step does
% not affect the original struct.
virtCpdsExtraFieldsStruct = rmfield( ...
    virtualCpdsStruct, EXPECTED_VIRTUAL_CPDS_FIELDS);
extraFieldnames = fieldnames(virtCpdsExtraFieldsStruct);
% Now, for some reason, we take the extraFieldsStruct and create a cell
% array containing the exact same data. I'll be looking to remove this step
% once I understand why we do this.
virtCpdsExtraCell = cell(length(extraFieldnames), 2);
for i=1:length(extraFieldnames)
    efname = extraFieldnames{i};
    virtCpdsExtraCell{i, 1} = efname;
    virtCpdsExtraCell{i, 2} = virtCpdsExtraFieldsStruct.(efname);
end

%%%%%%%%%%%%%%% START UNPACK VIRTUAL COMPOUND STRUCT %%%%%%%%%%%%%%%%%%%%%%
% Now we completely unpack the virtual compound structure. I think we have
% to do this so that parfor can efficiently pass this data to the workers.
cpdIds = virtualCpdsStruct.cpdId;
if length(virtualCpdsStruct.solubility) == 1
    solubility = repmat(virtualCpdsStruct.solubility, Ncpds, 1);
else
    solubility = virtualCpdsStruct.solubility;
end

if length(virtualCpdsStruct.solubility_pH) == 1
    solubility_pH = repmat(virtualCpdsStruct.solubility_pH, Ncpds, 1);
else
    solubility_pH = virtualCpdsStruct.solubility_pH;
end

if length(virtualCpdsStruct.SF) == 1
    sol_factor = repmat(virtualCpdsStruct.SF, Ncpds, 1);
else
    sol_factor = virtualCpdsStruct.SF;
end

if length(virtualCpdsStruct.Peff) == 1
    peff = repmat(virtualCpdsStruct.Peff, Ncpds, 1);
else
    peff = virtualCpdsStruct.Peff;
end

if length(virtualCpdsStruct.logP) == 1
    logP = repmat(virtualCpdsStruct.logP, Ncpds, 1);
else
    logP = virtualCpdsStruct.logP;
end

if length(virtualCpdsStruct.Deff) == 1
    deff = repmat(virtualCpdsStruct.Deff, Ncpds, 1);
else
    deff = virtualCpdsStruct.Deff; 
end

if length(virtualCpdsStruct.pKa_a) == 1
    acidic_pka = repmat(virtualCpdsStruct.pKa_a, Ncpds, 1);
else
    acidic_pka = virtualCpdsStruct.pKa_a;
end

if length(virtualCpdsStruct.pKa_b) == 1
    basic_pka = repmat(virtualCpdsStruct.pKa_b, Ncpds, 1);
else
    basic_pka = virtualCpdsStruct.pKa_b;
end

if length(virtualCpdsStruct.BP) == 1
    bp = repmat(virtualCpdsStruct.BP, Ncpds, 1);
else
    bp = virtualCpdsStruct.BP;
end

if length(virtualCpdsStruct.fup) == 1
    fup = repmat(virtualCpdsStruct.fup, Ncpds, 1);
else
    fup = virtualCpdsStruct.fup;
end

if length(virtualCpdsStruct.CLint_u) == 1
    clint = repmat(virtualCpdsStruct.CLint_u, Ncpds, 1);
else
    clint = virtualCpdsStruct.CLint_u;
end

if length(virtualCpdsStruct.mw) == 1
    mw = repmat(virtualCpdsStruct.mw,Ncpds,1);
else
    mw = virtualCpdsStruct.mw;
end
%%%%%%%%%%%%%%%% END UNPACK VIRTUAL COMPOUND STRUCT %%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% START UNPACK PROTOCOL TABLE %%%%%%%%%%%%%%%%%%%%%%%%
% First get protocols related to ADME
pTypes = protocolTable.type;
admeIdx = strcmpi(pTypes, 'absorption') ...
    | strcmpi(pTypes, 'distribution') ...
    | strcmpi(pTypes, 'elimination');
admeTable = protocolTable(admeIdx, :);
admeStruct = table2struct(admeTable);
% Need to change field name 'setting' to 'val'. Remarkably, the most
% straightforward way to do this is to copy the struct to a cell and then
% copy back.
admeStruct = cell2struct( ...
    struct2cell(admeStruct), {'param', 'type', 'val'});

admeParams = {admeStruct.param};
Kp_method = admeStruct(strcmpi('Kp_method', admeParams)).val;
fup_method = admeStruct(strcmpi('fup_method', admeParams)).val;
fut_method = admeStruct(strcmpi('fut_method', admeParams)).val;

% Next get target tissue.
pParams = protocolTable.param;
targetTissue = protocolTable.setting{strcmpi('target_tissue', pParams)};

% Build up a dose vector. The dose vector isn't used in this function, and
% it is recreated in runVirtualSimulation. It is only here so it can be
% returned to the user and saved with the output. (I guess it's hard to
% pull it from runVirtualSimulation because that function is used in the
% parfor loop. It's not impossible, but I guess this was most
% straightforward.)
minIdx = strcmpi(pTypes, 'dose') & strcmpi(pParams, 'min_amount');
maxIdx = strcmpi(pTypes, 'dose') & strcmpi(pParams, 'max_amount');
nDoseIdx = strcmpi(pTypes, 'dose') & strcmpi(pParams, 'n_dose');
minDoseAmount = protocolTable.setting{minIdx};
maxDoseAmount = protocolTable.setting{maxIdx};
totalDoses = protocolTable.setting{nDoseIdx};
doseAmounts = logspace( ...
    log10(minDoseAmount), log10(maxDoseAmount), totalDoses);
%%%%%%%%%%%%%%%%%%%%%%% END UNPACK PROTOCOL TABLE %%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%% START SIMULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
virtCpdsStructArray = cell(Ncpds, 1);
results = cell(Ncpds, 1);
failedCpdIds = cell(Ncpds, 1);
failedErrorObjs = cell(Ncpds, 1);
idxSimSuccess = true(Ncpds,1);

% Set up the callback to keep track of simulation progress.
% dataQueue = parallel.pool.DataQueue;
% afterEach(dataQueue, @nUpdateWaitbar);

simStartTime = tic;
% wbar = waitbar(0, sprintf(['Completed 0 out of %i runs. Time elapsed: ' ...
%     '0 seconds'], Ncpds), ...
%     'Name', 'Virtual Screen Progress');
% waitbarCounter = 0;

%for i=1:Ncpds
    
    % ----------------- START BUILD INPUT STRUCTURES -------------------- %
    % Bob's function buildPBPKInputStruct requires two param structs, a
    % compound-specific struct and a species-specific struct. Build those
    % now. (Originally Bob had this written as a struct array with fields
    % 'param' and 'val', but then he needed a bunch of strcmp statements to
    % pick out the desired parameter when needed. Much simpler to keep this
    % a 1x1 struct.)

    % First build the compound specific struct.
    cmpdStruct = struct( ...
        'logP', logP, ...
        'solubility', solubility, ...
        'sol_pH', solubility_pH, ...
        'SF', sol_factor, ...
        'Peff', peff, ...
        'Deff', deff, ...
        'FaSSIF', 0, ...
        'mw', mw, ...
        'Acidic_pKa', acidic_pka, ...
        'Basic_pKa', basic_pka);

    % calculate logD from logP
    pKa_a = acidic_pka;
    if isnan(pKa_a)
        pKa_a = [];
    end
    pKa_b = basic_pka;
    if isnan(pKa_b)
        pKa_b = [];
    end
    cmpdStruct.logD = logDcalc(logP, 7.4, pKa_a, pKa_b);

    
    % Now build the species-specific input struct. Again, we'll simplify
    % from an Nx1 struct array to a 1x1 struct array.
    spStruct = struct( ...
        'species', species, ...
        'CLint', clint, ...
        'BP', bp, ...
        'fup', fup, ...
        'fub', fup./bp);  % Assume fub = fup/BP

    % Now build the PBPK input struct. (What is returned is actually a 1x3
    % cell array. The input struct is the 2nd element of that array.) Any
    % error messages are contained in the third entry - WHY DO WE IGNORE
    % THESE MESSAGES?
    pbpkInputCellArray = buildPBPKInputStruct( ...
        cmpdStruct, spStruct, admeStruct, '', '');
    pbpkInputStruct = pbpkInputCellArray{2};
    
    % We use the toolbox function parameterizePBPK, but we don't capture
    % the returned model because parameterizePBPK does not modify exported
    % models (for non-exported models it adds the variants directly to the
    % model.
    % The following function generates both a fasted and fed physiology for
    % each species listed in the input pbpkInputStruct. Thus,
    % pbpkVariantContent is a 2x2 cell array, where the first column are
    % just the names '<species> Fasted' and '<species> Fed', and the second
    % column contains the variant content cell arrays.
    [~, pbpkInputStruct, pbpkVariantContent] = parameterizePBPK( ...
        m0, pbpkInputStruct, 'fuinc_human', 1);
    kpPerfStruct = pbpkInputStruct.SpeciesSpecific.Kp_perf;
    kpPermStruct = pbpkInputStruct.SpeciesSpecific.Kp_perm;
    % ------------------ END BUILD INPUT STRUCTURES --------------------- %
    
    % -------------- START MERGE PBPK AND PD VARIANTS ------------------- %
    pbpkpdVariantContent = cell(1, 2);
    pbpkpdVariantContent{1} = 'pbpk_pd_variant';
    % For some reason we're only interested in the fasted physiologies.
    fastedIdx = contains(pbpkVariantContent(:, 1), 'Fasted', ...
        'IgnoreCase', true);
    pbpkpdVariantContent{2} = [ ...
        pbpkVariantContent{fastedIdx, 2}; pdVariant];
    
    % The virtual compounds might be defined by some PD-related parameters.
    % Before starting this loop we extracted those parameters into
    % virtCpdsExtraCell (storing the 'extra' params). Now we see which of
    % those parameters correspond to parameters defined on the PD variant.
    pdVariantTable = variantContentTable(pdVariant);
    pdExtraParams = intersect(virtCpdsExtraCell(:, 1), ...
        pdVariantTable.name, 'stable'); %#ok<PFBNS> 
    
    for j=1:length(pdExtraParams)
        % Uses a try block just to throw an error with a better error
        % message. The function setVariantContent is so simple we should
        % probably just write a version which throws more useful errors,
        % then we could avoid this try block.
        try
            idxPV = strcmpi(virtCpdsExtraCell(:, 1), pdExtraParams{j});
            pbpkpdVariantContent{2} = setVariantContent( ...
                pbpkpdVariantContent{2}, ...
                pdExtraParams{j}, virtCpdsExtraCell{idxPV, 2}(1));
        catch
            error(['Error setting potency value(s). Check parameter ' ...
                'names to ensure they are consistent with the model.']);
        end
    end
    % --------------- END MERGE PBPK AND PD VARIANTS -------------------- %

    % ------------------- START SET fut, fut_ic ------------------------- %
    % TO DO: FIRST LOOK TO SEE IF WE EVEN NEED TO CALCULATE fut! 
    % List all tissues with Kp values.
    kpTissues = {kpPerfStruct.Kp.Tissue};
    
    tisName = targetTissue;
    fu = NaN;
    fuIC = NaN;
    fuEC = NaN;
    if ~any(strcmpi(kpTissues, targetTissue))
        % The target tissue is not one of the tissues in the PBPK 
        % model. In this case we need to set the tissue PK parameters.
        % Bob assumes that the tissue name will include the type of
        % tissue (fast or slow) in parentheses after the tissue name,
        % something like "MyTissue (fast)".
        tis = strsplit(targetTissue, '(');
        if length(tis) == 1
            % Do nothing in order to preserve legacy behavior.
        elseif length(tis) ~= 2
            error(['Target tissue was "%s". The code expects the name ' ...
                'of the target tissue be followed by the "type" ' ...
                'of tissue in parentheses, where "type" is one of ' ...
                'fast, slow, or all. Thus, it should look something ' ...
                'like "NewTissue (slow)".'], targetTissue);
        else
            
            kpStruct = getCompoundKP( ...
                DrugProps.empty, species, ...
                'logP', cmpdStruct.logP, ...
                'acidic_pka', pKa_a, ...
                'basic_pka', pKa_b, ...
                'fup', spStruct.fup, ...
                'b2p', spStruct.BP, ...
                'method', Kp_method,...
                'fup_type', fup_method,...
                'fut_type', fut_method);

            tisName = strtrim(strrep(tis{1}, ' ', ''));
            tisType = strtrim(strrep(strrep(tis{2}, ')', ''), ' ', ''));
            [Kp, fu, fuIC] = weightedAverageKp(kpStruct, tisType);
            pbpkpdVariantContent{2} = setVariantContent( ...
                pbpkpdVariantContent{2}, ['Kp_' tisName], Kp);
            
            fuEC = NaN;
        end

    elseif strcmpi(targetTissue, 'blood')
        % Not yet supported
        error(['Target tissue cannot be "blood". This case will be ' ...
            'implemented when the need arises.']);

    elseif strcmpi(targetTissue, 'plasma')
        % Not yet supported
        error(['Target tissue cannot be "plasma". This case will be ' ...
            'implemented when the need arises.']);

    else
        idxPerf = strcmpi(targetTissue, {kpPerfStruct.Kp.Tissue});
        idxPerm = strcmpi(targetTissue, {kpPermStruct.Kp.Tissue});
        fu = kpPerfStruct.Kp(idxPerf).fut;
        fuEC = kpPermStruct.Kp(idxPerm).fut;
        fuIC = kpPermStruct.Kp(idxPerm).fut_ic;
    end
    
    % set fu values if required
    fuVarNames = {['fu_' tisName], ['fu_ic_' tisName], ['fu_ec_' tisName]};
    fuVals = [fu fuIC fuEC];
    pbpkpdVariantTable = variantContentTable(pbpkpdVariantContent);
    for j=1:length(fuVarNames)
        fuName = fuVarNames{j};
        if ismember(fuName, pbpkpdVariantTable.name)
            pbpkpdVariantContent{2} = setVariantContent( ...
                pbpkpdVariantContent{2}, fuName, fuVals(j));
        end
    end
    % -------------------- END SET fut, fut_ic -------------------------- %
    
    % ------------------- START RUN SIMULATION -------------------------- %
    m = applyVariant(m0, pbpkpdVariantContent{2});
    [resi, simData, failed, me] = runSingleCmpdSim( ...
        m, cpdIds{1}, protocolTable, species, profiles2Save);
    
    if failed
        failedCpdIds = cpdIds{1};
        failedErrorObjs = me;
        idxSimSuccess = 0;
    else
        results = resi;
    end
    % -------------------- END RUN SIMULATION --------------------------- %

    % --------------- START STORE VIRTUAL CPD PARAMS -------------------- %
    piS = pbpkInputStruct;
    spS = piS.SpeciesSpecific;
    vcdi = struct( ...
        'vc_id', cpdIds, ...
        'mw', piS.mw, ...
        'solubility', piS.sol_mgmL, ...
        'SF', piS.SF, ...
        'SR', piS.SR, ...
        'logP', piS.logP, ...
        'logD_pH74', piS.logD_pH74, ...
        'Peff', piS.Peff0, ...
        'Deff', piS.Deff, ...
        'BP', spS.B2P, ...
        'fup', spS.fup, ...
        'fub', spS.fub, ...
        'CLint_u', spS.CLint./spS.CLint_fuinc, ...
        'CL_bl_mLminkg', spS.CL_tot_bl.*1000./(60*bw_kg), ...
        'CL_pl_mLminkg', spS.CL_tot_pl.*1000./(60*bw_kg), ...
        'CLren_pl_mLminkg', spS.CL_ren_pl.*1000./(60*bw_kg), ...
        'CL_pctLBF', spS.CLh_pctLBF, ...
        'Vss_pl_Lkg', spS.Kp_perf.Vss_Lkg, ...
        'Vss_bl_Lkg', spS.Kp_perf.Vss_Lkg/ spS.B2P, ...
        'Thalf_clv_hr', (log(2)*spS.Kp_perf.Vss_L)/ spS.CL_tot_pl);

    vcdi.acidic_pKa = NaN;
    if ~isempty(piS.Acidic_pKa)
        vcdi.acidic_pKa = piS.Acidic_pKa;
    end

    vcdi.basic_pKa = NaN;
    if ~isempty(piS.Basic_pKa)
        vcdi.basic_pKa = piS.Basic_pKa;
    end

    for j=1:length(virtCpdsExtraCell)
        vcdi.(virtCpdsExtraCell{j, 1}) = virtCpdsExtraCell{j, 2}(1);
    end
    
    try
        vcdi.(['Kp_' tisName]) = getVariantContent( ...
            pbpkpdVariantContent, ['Kp_' tisName]);
    catch
    end

    try 
        vcdi.(['fu_' tisName]) = getVariantContent( ...
            pbpkpdVariantContent, ['fu_' tisName]);
    catch
    end

    try
        vcdi.(['fu_ic_' tisName]) = getVariantContent( ...
            pbpkpdVariantContent, ['fu_ic_' tisName]);
    catch
    end

    try 
        vcdi.(['fu_ec_' tisName]) = getVariantContent( ...
            pbpkpdVariantContent, ['fu_ec_',tisName]);
    catch
    end
    
    % added lines to get the weighted average Kp and fut
    try
    vcdi.Kpavg = calculateWeightedAverageKp([kpPerfStruct.Kp.Kp], ...
        kpPermStruct.Details.tis);
    vcdi.KpTable = struct2table(kpPerfStruct.Kp);
    catch
    end
    % added line to get the fut
    try
        vcdi.fut_method1 = calculateFut(kpPermStruct.Details.tis, ...
            vcdi.fub,...
            vcdi.Vss_bl_Lkg*1000*70);
        vcdi.fut_method2 = vcdi.fup / vcdi.Kpavg;
    catch
    end


    virtCpdsStructArray = vcdi;
    % ---------------- END STORE VIRTUAL CPD PARAMS --------------------- %

    %send(dataQueue, i);
%end

%close(wbar);

% This next line converts the cell array of structs into a struct array.
%virtCpdsStructArray = [virtCpdsStructArray{:}];

results = results(logical(idxSimSuccess));
%results = [results{:}];

failedCpdIds = failedCpdIds(~logical(idxSimSuccess));
failedErrorObjs = failedErrorObjs(~logical(idxSimSuccess));
failedRunsStruct = struct( ...
    'cpdId', failedCpdIds, 'ME', failedErrorObjs);

if ~isempty(saveFilename)
    save(saveFilename, ...
        'virtCpdsStructArray', 'results', ...
        'doseAmounts', 'failedRunsStruct');
end


figure('DefaultAxesFontSize', 8);
tiledlayout('flow');
for ii = 1:length(profiles2Save)
    nexttile;
    plot(simData(end).sim_data.Time, simData(end).sim_data.Data(:, ii), ...
        'LineWidth', 2);
    xlabel('Time, hr');
    title(sprintf('%s-%s',profiles2Save{ii}, ...
        simData(end).sim_data.DataInfo{ii}.Units));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% END SIMULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%% START INNER HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%
% function nUpdateWaitbar(~)
%     waitbarCounter = waitbarCounter + 1;
%     waitbar(waitbarCounter/Ncpds, wbar, ...
%         sprintf(['Completed %i out of %i compounds.\nTime elapsed: ' ...
%             '%0.2f seconds'], waitbarCounter, Ncpds, toc(simStartTime)));
% end
%%%%%%%%%%%%%%%%%%%%% END INNER HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%% START OUTER HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%
function bodyWeight = speciesWeight(species)
    switch lower(species)
        case 'mouse'
            bodyWeight = 0.025;
        case 'rat'
            bodyWeight = 0.25;
        case 'dog'
            bodyWeight = 10;
        case 'minipig'
            bodyWeight = 14.2;
        case 'monkey'
            bodyWeight = 4;
        case 'human'
            bodyWeight = 70;
        otherwise
            error('Unexpected species: %s', species);
    end
end
%%%%%%%%%%%%%%%%%%%%% END OUTER HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%

