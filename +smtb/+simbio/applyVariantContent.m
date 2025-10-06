function applyVariantContent(m, v)
arguments
    m {checkValidationOfM(m)};
    v {checkValidationOfV(v)};
end
% applyVariantContent applies the variant content to the simbiology model m
% (m can either be a exported or nonexported model)
% Currently works for the exported model. Need to include updates for the
% non exported model
% Parameters:
% 	 m: exported or nonexported simbiology model
% 	 v: simbiology variant object or a cell array of variant information
% Date modified: 04-Oct-2022
% File modified by SMTB user: jbp17697
% File created by SMTB user
% Examples:
%
vTb = smtb.simbio.variantContentTable(v); %convert variant content to table
paramNames = string(vTb.Name);
paramValues = vTb.Value;

if isa(m, 'SimBiology.export.Model')
    updateValuesExportedModel(m, paramNames, paramValues);
else
    updateValuesNonExportedModel(m, paramNames, paramValues);
end

end

function updateValuesExportedModel(m, paramNames, paramValues)
allParamModel = string({m.ValueInfo.Name});
if ~all(ismember(paramNames, allParamModel))
    errID = 'myComponent:inputError';
    msgtext = ['input parameters structure contains'...
        ' parameter name not present in the model \n.'];
    ME = MException(errID,msgtext);
    throw(ME)
end

for ii = 1:length(paramNames)
    m.InitialValues(m.getIndex(paramNames(ii))) = paramValues(ii);
end

end

function checkValidationOfV(v)
if isa(v, "SimBiology.Variant")
elseif iscell(v)
    if length([v{:}]) ~= length(v) * 4
        error(['v must be a Simbiology Variant object or a '...
            'cell array with each element as 1x4 cell']);
    end

else
    error(['v must be a Simbiology Variant object or a '...
        'cell array with each element as 1x4 cell'])
end
end

function checkValidationOfM(m)
if ~(isa(m,'SimBiology.export.Model') || isa(m,'SimBiology.Model'))
    error (['Model m must be a simbiology model or a simbiology'...
        'exported model']);
end
end