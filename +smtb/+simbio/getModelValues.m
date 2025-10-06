function infoTable = getModelValues(m,sd)
arguments
    m {checkValidationM(m)};
    sd = [];
end
% getModelValues provides the initial values assigned to the parameters,
% compartments and species of the simbiology exported or non exported 
% model
%
% Parameters:
%   m ; % simbiology exported or non-exported model
%   sd; % simbiology simulation data object
% 
% Outputs:
%   infoTable; % output information table about the model values 
%
%

if isa(m, 'SimBiology.Model')
    infoTable  = getModelValuesM(m);
else
    infoTable = getModelValuesMexp(m);
end

disp(sd);

end


function checkValidationM(m)
if ~(isa(m, 'SimBiology.Model') || isa(m, 'SimBiology.export.Model'))
    error(['The model should be ' ...
        'simbiology exported or non exported model'])
end
end

function infoTable = getModelValuesM(m)

% Save info about whether parameter and compartment are constant or not
parameterConstantInfo = [m.Parameters.Constant];
compartmentsConstantInfo = [m.Compartments.Constant];


[m.Parameters.Constant] = deal(false);
[m.Compartments.Constant] = deal(false);




p = m.Parameters.get;
c = m.Compartments.get;
s = m.Species.get;
names = {s.Name, p.Name, c.Name}';
%[names, idxs] = sort(names);
values = [p.Value, c.Value, s.Value]';
%values = values(idxs);
units = {s.Units, p.Units, c.Units}';
%units =  units(idxs);
infoTable = table();
infoTable.Names = names;
infoTable.Values = values;


cs = m.getconfigset();
originalStopTime = cs.StopTime;
cs.StopTime = 0;
sd = sbiosimulate(m);

%[~,idxs] = sort(sd1.DataNames);
updatedValues = sd.Data';
infoTable.UpdatedValues = updatedValues;
%infoTable.DataNames = names;

infoTable.Units = units;
infoTable.Type = cellfun(@(x)x.Type, sd.DataInfo, ...
    'UniformOutput', false);

% Put back the model in the original state
cs.StopTime = originalStopTime;
pCons = num2cell(parameterConstantInfo);
[m.Parameters.Constant] = pCons{:} ;
cCons = num2cell(compartmentsConstantInfo);
[m.Compartments.Constant] = cCons{:};

end


function infoTable = getModelValuesMexp(mexp)

values = [mexp.ValueInfo.InitialValue];
names = {mexp.ValueInfo.Name};
units = {mexp.ValueInfo.Units};
type = {mexp.ValueInfo.Type};

infoTable = table();
infoTable.Names = names';
infoTable.Values = values';


originalStopTime = mexp.SimulationOptions.StopTime;

mexp.SimulationOptions.StopTime = 0;
sd = simulate(mexp);
updatedValues = nan(length(values), 1);
idx = contains(names, regexpPattern(cellfun(@(x)(sprintf('^%s$',x)), ...
    sd.DataNames', 'UniformOutput', false)));
updatedValues(idx) = sd.Data';


infoTable.UpdatedValues = updatedValues;

infoTable.Units = units';
infoTable.Type = type';

% Put back the model in the original state
mexp.SimulationOptions.StopTime = originalStopTime;




end