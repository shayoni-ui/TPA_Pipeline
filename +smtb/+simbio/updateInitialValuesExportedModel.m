function  updateInitialValuesExportedModel(mexp, pars)
arguments
    mexp,
    pars,
end
% updateInitialValuesExportedModel  ...
% Parameters:
% 	 mexp:
% 	 pars:
% Date modified: 27-Sep-2022
% File modified by SMTB user: jpb17697
% File created by SMTB user: Jaimit Parikh
% Example:
%

% check if the name of the parameters to be updated are present in
% the model
parameterNames = string(fields(pars))';
allParametersModel = string({mexp.ValueInfo.Name});
if ~all(ismember(parameterNames, allParametersModel))
    errID = 'myComponent:inputError';
    msgtext = ['input parameters structure contains parameter name ' ...
        'not present in the model \n.'];
    ME = MException(errID,msgtext);
    throw(ME)
end

for ii = parameterNames
    mexp.InitialValues(mexp.getIndex(ii)) = pars.(ii);
end

end
