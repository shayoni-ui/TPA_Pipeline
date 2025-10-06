classdef sampleSet < handle
    % sampleSet class provides objects to generate samples of arbitrary
    % dimension and distributions. Currently we assume independence between
    % the parameters. We can use any distribution supported by makedist
    %
    % sampleSet Properties:
    %   NP - number of parameters to generate the samples
    %   NS - number of samples along each dimension to be sampled. 
    %       For sobol the number of final samples are (NP + 2) * NS
    %   MethodType - string for the method type to be used for sampling.
    %       Current options are "lhs", "montecarlo", "structured", "sobol"
    %   ParamNames - Name of parameters
    %   DesiredParamDists - Desired distribution of the parameters
    %   Sobolgenerator - sobolgenerator
    %   Samples - final samples generated
    %   DateGenerated - Date the distribution was generated
    %   User - Username who generated the distribution
    %   SampleLevel - The number of times the samples were generated
    %   TotalSamples - The total number of samples generated
    %   Seed - Random seed for reproducibility (int or 'default')
    %
    % sampleSet Methods:
    %   setParamDist - method to set the marginal distribution for each 
    %           parameter
    %   generateSamples - method to generate the samples using desired
    %       method and distributions
    %
    % Date modified: 03-Oct-2022 
    % File modified by SMTB user: jbp17697 
    % File created by SMTB user: Jaimit Parikh
    %
    % See also MAKEDIST
    %
    % TODO: include possibility to sample with some correlation structure
    %     : option to build distribution from given data for a parameter  
    %
    % Examples:
    %
    %

    properties
        NP; % number of parameters
        NS; % Number of base samples along each dimension
        MethodType; % supported methods include structured, sobol, lhs, mc
        ParamNames; % 
        DesiredParamDists;
        Samples;
        DateGenerated;
        User;
        SampleLevel;
        TotalSamples;
        Seed; % Seed for reproducibility
    end
    properties(Access = 'private')
        Sobolgenerator;
    end

    methods

        function obj = sampleSet(methodType)
            % constructor
            arguments
                methodType {mustBeMember(methodType, ["sobol",...
                    "lhs", "montecarlo",...
                    "structured"])} = "lhsdesign";
            end
            obj.MethodType = methodType;
            obj.DesiredParamDists = [];
            obj.SampleLevel = 0;
            obj.TotalSamples = 0;
        end

        setParamDist(obj, paramName, paramDistType, distParamNamesValues)

        samples = generateSamples(obj, ns, seed);

        f = visualizeDist(obj, lvl);

    end

    methods(Access = private)

        samples = getSamplesSobol(obj);

        samples = getSamplesLHS(obj);


    end

end