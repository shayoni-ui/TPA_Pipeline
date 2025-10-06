function nca  = ncaQuickCalculation(t,c,nva)
arguments
    t;
    c;
    nva.auctype {mustBeMember(nva.auctype,{'linuplogdown', ...
        'linlog'})} = 'linuplogdown';
    nva.llq = 0;
    nva.numPoints = 5;
    nva.doseLength = 0;
end
% ncaQuickCalculation is a condensed version of ncaCalculation to
% improve calculation speed for large n PBPK simulations  ...
%
% THIS IS ONLY RECOMMENDED FOR USE WITH SIMULATED PROFILES
%
% Parameters:
% 	 t: time for a single dose
% 	 c: conc for a single dose
% 	 nva.auctype {mustBeMember(nva.auctype,{'linuplogdown', ...
%         'linlog'})} = 'linuplogdown';
%    nva.llq = 0; Lower limit of quantization
%    nva.numPoints = 5; number of points to calculate lambdaz
%    nva.doseLength = 0;
% Outputs
% 	 nca: Output is a class NCA with AUC, AUCinf, AUMC, AUMCinf, t1/2,
%         Cmax, Tmax and others
% Date modified: 09-Oct-2022
% File modified by SMTB user: jbp17697
% File created by SMTB user
% Examples:
%

if length(t) ~= length(c), error(['Length of time and' ...
        ' concentration must be the same']), end

t = t(:); t = t(~isnan(c)); if isempty(t), t = 0; end
c = c(:); c = c(~isnan(c)); if isempty(c), c = 0; end
if length(t) ~= length(c), error(['Time and concentration ' ...
        'vectors must be the same length']); end

[t,idx] = sort(t);
c       = c(idx);

nca = smtb.pbpk.NCA;
nca.LLQ = nva.llq;
nca.DoseLength = nva.doseLength;

t = t(c > nca.LLQ);
c = c(c > nca.LLQ);

[nca.Cmax,idx] = max(c); % add logic to make sure cmax is 
% after end of infusion
nca.Tmax       = t(idx);
nca.Clast      = c(find(~isnan(c),1,'last'));
nca.Tlast      = t(find(~isnan(c),1,'last'));

[nca.AUC,nca.AUMC] = trapzlog(t,c,nva.auctype);

nca.MRT = nca.AUMC / nca.AUC - nca.DoseLength/2;

if length(t) >= nva.numPoints
    idx = length(t) - nva.numPoints + 1:1:length(t);

    lm   = fitlm(t(idx),log(c(idx)));
    nca.lambdaZ      = -lm.Coefficients.Estimate(2);
    nca.lambdaZ_pts  = lm.NumObservations;
    nca.lambdaZ_rng  = [min(lm.Variables.x1) max(lm.Variables.x1)];
    nca.Rsq          = lm.Rsquared.Ordinary;
    nca.Rsq_adjusted = lm.Rsquared.Adjusted;

    nca.AUCinf  = nca.AUC + c(end)./nca.lambdaZ;
    nca.AUMCinf = nca.AUMC + c(end)*t(end)/nca.lambdaZ + ...
        c(end)/nca.lambdaZ^2;
    nca.MRTinf  = nca.AUMCinf / nca.AUCinf - nca.DoseLength/2;
    nca.t12     = log(2)/nca.lambdaZ;
    nca.t12_CI  = [log(2)./-(lm.Coefficients.Estimate(2) - ...
        1.96 * lm.Coefficients.SE(2))...
        log(2)./-(lm.Coefficients.Estimate(2) + ...
        1.96 * lm.Coefficients.SE(2))];
    nca.t12_CI(nca.t12_CI < 0) = Inf;
end
end
