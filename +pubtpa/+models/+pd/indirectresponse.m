function mPD = indirectresponse()
m1 = sbiomodel('indirectResponse');
cmpProps = {'Name', 'Value', 'Units', 'Constant'};
cmps = {
    {'pd', 1, 'milliliter', 1},...
    };
addcompartments(m1, cmps, cmpProps);

spProps = {'CompName', 'Name', 'Value', 'Units'};
species = {
    {'pd', 'P', 1, 'nanomole/liter'},...
    {'pd', 'TE', 0, 'molecule'}
    };
addspecies(m1, species, spProps)

% 1. P <-> null, null -> P,
rnSet1 = {
    {'P -> null', 'kdeg_P * P','R_degradation_P'},...
    {'TE -> P + TE', 'ksyn_P * drug_effect','R_synthesis_P'},...
    };
prSet1 = {
    {'kdeg_P',  1, "1 / hour", "deg", 1}, ...
    {'ksyn_P', 1, "nanomole / (liter * hour)", "syn", 1},...
    {'drug_effect', 1, 'dimensionless', "", 0},...
    {'oneTE', 1, 'molecule', "", 1},...
    {'Imax', 1, "dimensionless", "", 1},...
    };

% Reactions
rnProps = {'Reaction', 'ReactionRate', 'Name'};
reactions = rnSet1;
addreactions(m1, reactions, rnProps);

% Parameters
paramProps = {'Name', 'Value', 'Units', 'Tag', 'Constant'};
params = prSet1;
addparameters(m1, params, paramProps);

addrule(m1,'drug_effect = 1 - (Imax * (pd.TE/oneTE))', ...
    'RuleType','RepeatedAssignment');


m1.getequations
mPD = m1;
clearvars -except mPD;

cs = getconfigset(mPD);
cs.CompileOptions.DimensionalAnalysis = true;
cs.CompileOptions.UnitConversion = true;
cs.TimeUnits = 'hour';

addparameter(mPD,'clTE','Value',0,'ValueUnits','1/hour');
addreaction(mPD,'TE -> null','Name',...
    'R_clTE','ReactionRate','TE * clTE');
end

