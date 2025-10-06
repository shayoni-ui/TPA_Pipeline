function sampleobj = getSamples(nva)
arguments
    nva.NS = 8;
    nva.methodType = 'lhs';
    nva.paramDist = {...
        'a', 'normal', {'mu', 0, 'sigma', 0.1},...
        'b', 'uniform', {'Lower', 0, 'Upper', 5} ...
        };
    nva.seed = '';

end
%% getSamples allows to generate samples of arbitrary 
% distributions
% getSamples
% Parameters:
%   NS - number of samples
%   methodType = "lhs", "montecarlg", "structured", "sobol"
%   seed = seed for reproducibility
% Outputs:
%   sampleObj - Resultsing samples object
% 
% Examples:
%   getSamples(4, 'strucutred', {}
sampleobj = smtb.gsa.sampleSet(nva.methodType);

if isstruct(nva.paramDist)
    paramDistArr = struct2cell(nva.paramDist')';
    paramDistArr = cellfun(@expandStruct, paramDistArr, ...
        'UniformOutput',false);
    paramDistArr = reshape(paramDistArr', 1, []);
    sampleobj.setParamDist(paramDistArr{:})
elseif iscell(nva.paramDist)
    paramDistArr = nva.paramDist;
    if contains(class(nva.paramDist{1}), 'prob')
       sampleobj.DesiredParamDists = paramDistArr;
    else
       sampleobj.setParamDist(paramDistArr{:});
    end
else
    error("Unexpected type for the paramDistArr")
end



sampleobj.generateSamples(nva.NS, nva.seed);


end

function xo = expandStruct(x)
if isstruct(x)
    xo = namedargs2cell(x);
else
    xo = x;
end
end