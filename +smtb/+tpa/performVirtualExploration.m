function sd = performVirtualExploration(mexp, d, nva)
arguments
    mexp SimBiology.export.Model;
    d;
    nva.variantArr cell = {}; % variant cell array
    nva.TeDuration double = []; % Target engagement duration array 
    %   for mode 1
    nva.selectNames  = ["plasma.drug", "TargetTissue.drugU", "pd.CGRP"];
    nva.shrinkData double = 0; % shrink data to save last shrinData doses
    % If 0 would save the entire profile.
end
% performVirtualExploration allows to run in parallel the exported PKPD
% for mTPA analysis. 
% Mode1 Analysis: Pass exported model, dose object and TeDuration array
% Mode2 Analysis: Pass the exported model, dose object and variant array
%
% Parameters:
%   mexp SimBiology.export.Model;
%   d;
%   nva.variantArr cell = {}; % cell array of simbiology variant
%   nva.TeDuration double = []; % Target engagement duration array 
%                                 for mode 1
%   nva.selectNames  = ["plasma.drug", "TargetTissue.drugU", "pd.CGRP"];
%   nva.shrinkData double; % shrinks data to save only
%                               last shrinkData doses for long sim
% Outputs:
%   sd - Simbiology simulation data object or similar structure
%
% Date modified: 10-Oct-2022 
% File modified by SMTB user: jbp17697 
% File created by SMTB user: Jaimit Parikh
%
% Examples: 
%   sd = performVirtualExploration(mexp,d, 'variantArr', variantArr)
%   sd = performVirtualExploration(mexp,d, 'TeDuration', Teduration);
%
variantArr = nva.variantArr; TeDuration = nva.TeDuration;
selectNames = nva.selectNames; shrinkData = nva.shrinkData;

if isempty(variantArr)
    variantArr = cell(1, length(TeDuration));
else
    TeDuration = double.empty(length(variantArr), 0);
end

% Set up callback for displaying progress of the simulations
D = parallel.pool.DataQueue;
N = length(variantArr);
p = 1;
afterEach(D, @nUpdateWaitbar);

sd = cell(1, N);

parfor ii = 1:N
    sd{ii} = performSingleSim(mexp,d, variantArr{ii}, ...
        TeDuration(ii,:), selectNames, shrinkData);
    send(D, ii);
end

% Inner helper function
    function nUpdateWaitbar(~)
        %waitbar(p/N, h);
        tt = 100 * p/N;
        if mod(tt, 10) == 0
            fprintf('Done ... %2.0f %% \n', tt);
        end
        p = p + 1;
    end
end

function sd = performSingleSim(mexp,d, variant, TeDuration,...
    selectNames, shrinkData)
arguments
    mexp;
    d;
    variant;
    TeDuration;
    selectNames =["plasma.drug"];
    shrinkData = 10; % shrinks data to save only last 'shrinkData' doses
end

if ~isempty(variant)
    % Get initial values from the variant
    x0 = smtb.simbio.initialValuesFromVariant(mexp, variant);
end

if ~isempty(TeDuration)
    clTE = mexp.ValueInfo(string({mexp.ValueInfo.Name}) == ...
        'clTE').InitialValue;
    d.Amount = clTE*TeDuration(1)*TeDuration(2);
    d.Rate = clTE*TeDuration(1);
end

if ~isempty(variant)
    sd = simulate(mexp,x0, d);
else
    sd = simulate(mexp, d);
end
sd = selectbyname(sd, selectNames);


% Resample to save last x doses provided in shrinkData
if shrinkData
    t = sd.Time(sd.Time > d.Interval * (d.RepeatCount - shrinkData));
    sd = resample(sd, t); 
end

end