function setParamDist(obj, paramName, paramDistType, ...
    distParamNamesValues, bounds)
% setParamDistallows to set the desired distribution of each parameters
% Any distribution to MATLAB makedist allowed here
%
% See also makedist

arguments
    obj;
end
arguments(Repeating)
    paramName;
    paramDistType;
    distParamNamesValues;
    bounds;
end

obj.ParamNames = paramName;
obj.NP = length(paramName);

% create marginal distributions for each parameters as defined
% in the arguments using makedist

desiredParamDists = cellfun(@(x, y)makedist(x, y{:}),paramDistType, ...
    distParamNamesValues, 'UniformOutput',false);

% created truncated distributions from marginals if the bounds are provided
obj.DesiredParamDists = cellfun(@truncate2, desiredParamDists, bounds,...
    'UniformOutput',false);

end

function tDist = truncate2(dist, bounds)
    if ~isempty(bounds)
        tDist = truncate(dist,bounds{:});
    else
        tDist = dist;
    end

end