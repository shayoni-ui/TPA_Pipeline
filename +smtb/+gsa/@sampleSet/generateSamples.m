function samples = generateSamples(obj, ns, seed)
% generateSamples(obj, ns) is a method to generate and store the samples
% Final samples are stored in the samples property of the obj
arguments
    obj;
    ns;
    seed ='';
end

if ~isempty(seed)
    obj.Seed = seed;
end
obj.TotalSamples = obj.TotalSamples + ns;
obj.NS = ns;
if obj.SampleLevel < 1
switch obj.MethodType
    case {"sobol", "structured"}
        samples = getSamplesSobol(obj);

    case "lhs"
        samples = getSamplesLHS(obj);

    case "montecarlo"
        samples = getSamplesMonteCarlo(obj);

    otherwise
        error("The method type can be either sobolset or lhsdesign")
end
else
    if ~isempty(obj.Sobolgenerator)
        samples = obj.Sobolgenerator(obj.TotalSamples+1:...
            obj.TotalSamples+ns,:);
    else
        error(['Generating different set of samples only supported ',...
            'for "Sobol", "Structured", options type']);
    end
end

if ~isempty(obj.DesiredParamDists)
    % apply inverseTransformSampling to each columns of the samples
    % to convert the samples to the desired distributions given by
    % paramDists
    originalParamDist = cell(1, obj.NP);
    [originalParamDist{:}] = deal(makedist('Uniform'));
    samples = splitapply(...
        @(a,b,c)inverseTransformSampling(a, b{:}, c{:}), ...
        samples, originalParamDist ,obj.DesiredParamDists,...
        1:obj.NP);

end
obj.SampleLevel = obj.SampleLevel + 1;

obj.Samples{obj.SampleLevel} = samples;

obj.DateGenerated{obj.SampleLevel} = datetime('today');
if ispc
    user = getenv('username');
else
    user = getenv('USER');
end
obj.User{obj.SampleLevel} = user ;

end


function transformedSamples = inverseTransformSampling(samples, ...
    sampleDist, desiredDist)
    arguments
        samples;
        sampleDist;
        desiredDist;
    end
    %% transformedSamples = inverseTransformSampling(samples, sampleDist, 
    % desiredDist) returns the transformedSamples with the desired
    % distribution given orginal samples sampled from sampleDist
    
    transformedSamples = desiredDist.icdf(sampleDist.cdf(samples));
    

end