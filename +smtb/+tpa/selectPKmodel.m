function [m,mexp] = selectPKmodel(modelType, nva)
arguments
    modelType {mustBeMember(modelType, {'PBPK',...
        'Compartmental', 'LargeMolecule'})};
    nva.newCompartment {validateNewCompartment(nva.newCompartment)};
    nva.exportFlag = true;
end
% selectPKmodel allows to generate a PBPK, compartmental or
% large molecule PK model
%
% Parameters:
%   modelType: {mustBeMember(modelType,
%          {'PBPK','Compartmental', 'LargeMolecule'})}
%   nva.newCompartment: Structure to add a new compartment
%          Required fields: name, capacity, Q
%   nva.exportFalg; flag to check whether to export the model or not
%
% Outputs:
%   m: SimBiology Model
%   mexp: SimBiology exported model
%
% Date modified: 19-Oct-2022
% File modified by SMTB user: jbp17697
% File created by SMTB user: Jaimit Parikh
%
% Examples:
%   s.name = 'targetTissue'; s.capacity = 1; s.Q = 100;
%   [m, mexp] = selectPKmodel('PBPK', 'newCompartment', s)
%   [m, mexp] = selectPKmodel('PBPK', 'newCompartment', s,...
%                              'exportFlag', false)
switch modelType
    case 'PBPK'
        m = smtb.pbpk.buildPBPKmodelTemplate();
        if ~isempty(nva.newCompartment)
            addNewCompartment(m, nva.newCompartment);
        end
        cs = m.configset;
        cs.TimeUnits = 'hour';
        cs.SolverOptions.AbsoluteTolerance = 1e-8;
        cs.SolverOptions.RelativeTolerance = 1e-6;
        cs.CompileOptions.UnitConversion = true;
        cs.CompileOptions.DimensionalAnalysis = true;
        verify(m,cs);
    case 'Compartmental'
        m = smtb.pbpk.buildPBPKmodelTemplate("pkmodel",...
            "compartmental");
    case 'LargeMolecule'
        m = Simbiology.Model.empty;
end

if nva.exportFlag
    mexp = export(m); accelerate(mexp)
else
    mexp = SimBiology.export.Model.empty;
end

end

function addNewCompartment(m, compartmentInfo)
tt = addcompartment(m, compartmentInfo.name, ...
    'Capacity',compartmentInfo.capacity,...
    'CapacityUnits','milliliter');%brain blood
addspecies(tt, 'drug', 'InitialAmount',0, ...
    'InitialAmountUnits','nanogram/milliliter');
addspecies(tt, 'drugU', 'InitialAmount',0, ...
    'InitialAmountUnits','nanogram/milliliter');

addparameter(m,sprintf('Q_%s', compartmentInfo.name), ...
    'Value',compartmentInfo.Q, 'ValueUnits', 'milliliter/hour');

addreaction(m, sprintf('null -> %s.drug',compartmentInfo.name), ...
    'ReactionRate', sprintf('Q_%s*blood.drug',compartmentInfo.name),...
    'Name', sprintf('R_drug_syn_%s',compartmentInfo.name));
addreaction(m, sprintf('%s.drug -> null',compartmentInfo.name), ...
    'ReactionRate', sprintf('Q_%s * B2P/Kp_Spleen * %s.drug',...
    compartmentInfo.name, compartmentInfo.name), ...
    'Name', sprintf('R_drug_deg_%s',compartmentInfo.name));

addrule(m, sprintf('%s.drugU = %s.drug*fup/Kp_Spleen',...
    compartmentInfo.name, compartmentInfo.name), ...
    'RuleType', 'RepeatedAssignment');

end

function validateNewCompartment(x)
if isa(x, "struct")
    requiredFields = {'name', 'capacity', 'Q'};
    checkFields = isfield(x,requiredFields);
    if all(checkFields)
    else
        error('Missing required fields, \n%s', ...
            strjoin(requiredFields(~isfield(x, requiredFields)),...
            '\n'));
    end

else
    error('The newCompartment argument needs to be structure');
end
end



