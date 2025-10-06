
function [variant, qois] = generateVCs(pKparamName, ...
    pKparamValues, pKconstPars, pdParamNames,pdParamValues, pdConstPars)
arguments
    pKparamName cell; 
    pKparamValues double;
    pKconstPars struct;
    pdParamNames cell = {};
    pdParamValues double = [];
    pdConstPars struct = struct.empty;
end
% generateVCs allows to generate set of virtual compounds given
% the name and value of PK and the PD parameters
%
% Parameters:
%    pKparamName cell; % name of the PK parameters that are varied
%    pKparamValues double; values of the PK parameter that are varied
%    pKconstPars struct; values of PK parameters that are constant but
%       different than default
%    pdParamNames cell; % name of the PD parameters that are varied
%    pdParamValues double; % values of the PD parameters that are varied
% Outputs:
%    variants; % returns a simbiology variant info for
%       generated virtual compound
%    qois; addition quanities of interest added here
%
% Date modified: 04-Oct-2022
% File modified by SMTB user: jbp17697
% File created by SMTB user
% Examples:

% preallocating variant array
variant = cell(1, length(pKparamValues));
qois = cell(1, length(pKparamValues));

if isempty(pdParamValues)
    pdParamValues = double.empty(length(pKparamValues), 0);
end


parfor ii = 1:length(pKparamValues) % PARFOR
    [qoi, variantContent] = generateVariant(pKparamName, ...
        pKparamValues(ii, :), pKconstPars,...
        pdParamNames, pdParamValues(ii,:), pdConstPars);
    variant{ii} = variantContent;
    qois{ii} = qoi;

end
end

function [qois, variantContent] = generateVariant(pKparamName, ...
    pKparamValues, pKconstPars, pdParamNames, pdParamValues, pdConstPars)

pNameValue = cat(1, pKparamName, num2cell(pKparamValues));
if ~isempty(pKconstPars)
pKconstParsCell = namedargs2cell(pKconstPars);
pbpkArgs = [pNameValue(:)', pKconstParsCell(:)'];
else
    pbpkArgs = pNameValue(:)';
end
[inputs,v] = smtb.pbpk.getPBPKVariant(pbpkArgs{:});
vHumanFasted = v(strcmpi({v.Name},'Human Fasted'));
% TODO: add other quantities of interest (for now avg Kp and VssL added)
kpAvg = calculateWtAvgKp([inputs.SpeciesSpecific.Kp_perf.Kp.Kp],...
    inputs.SpeciesSpecific.Kp_perf.Details.tis);
VssL = inputs.SpeciesSpecific.Kp_perf.Vss_L;
qois.kpAvg = kpAvg; qois.VssL = VssL;

if ~isempty(pdParamNames)
for jj = 1:length(pdParamNames)
    cc = {'parameter', pdParamNames{jj}, 'Value',...
        pdParamValues(jj)};
    vHumanFasted.addcontent(cc);
end
end

if ~isempty(pdConstPars)
    for jj = fieldnames(pdConstPars)'
        cc = {'parameter', jj{1} , 'Value', pdConstPars.(jj{1})};
        vHumanFasted.addcontent(cc);
    end
end

variantContent = vHumanFasted;
end

function kpAvg = calculateWtAvgKp(kp, tis)
tis = [[tis(1:13).volume]' [tis(1:13).density]'];
wgts = tis(:,1).*tis(:,2)./(sum(tis(:,1).*tis(:,2)));
kpAvg = sum(kp'.*wgts);
end
