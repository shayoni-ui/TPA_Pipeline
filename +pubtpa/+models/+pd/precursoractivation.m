function mPD  = precursoractivation()
m1 = sbiomodel('precursoractivation');
cmpProps = {'Name', 'Value', 'Units', 'Constant'};
cmps = {
    {'pd', 1, 'milliliter', 1},...
    };
addcompartments(m1, cmps, cmpProps);

spProps = {'CompName', 'Name', 'Value', 'Units'};
species = {
    {'pd', 'P', 100, 'nanomole/milliliter'},...
    {'pd', 'A', 1, 'nanomole/milliliter'},...
    {'pd', 'TE', 0, 'molecule'}
    };
addspecies(m1, species, spProps)

% 1. P -> A, null -> P,

kout = 0.2; % 1/hr
gamma = 0.2;
alpha = 100;
kact = kout / (alpha * (1 + gamma));
koth = gamma * kact;
%kin = kout;


rnSet1 = {
    {'P  -> A', 'koth * P','R_degradation_P'},...
    {'P + TE  -> TE + A', 'kact * P * drug_effect','R_activation_P'},...
    {'null -> P', 'kin','R_synthesis_P'},...
    };
prSet1 = {
    {'koth',  koth, "1 / hour", "deg", 1}, ...
    {'kin',  kout, "nanomole / (milliliter * hour)", "syn", 1}, ...
    {'kact', kact, "1 / hour", "deg", 1},...
    {'drug_effect', 1, 'dimensionless', "", 0},...
    {'oneTE', 1, 'molecule', "", 1},...
    {'Imax', 1, "dimensionless", "", 1},...
    };

% 2. P -> A A->null

rnSet2 = {
    {'A  -> null', 'kout * A','R_degradation_A'},...
    };
prSet2 = {
    {'kout',  kout, "1 / hour", "deg", 1}, ...
    };



% Reactions
rnProps = {'Reaction', 'ReactionRate', 'Name'};
reactions = [rnSet1, rnSet2];
addreactions(m1, reactions, rnProps);

% Parameters
paramProps = {'Name', 'Value', 'Units', 'Tag', 'Constant'};
params = [prSet1, prSet2];
addparameters(m1, params, paramProps);

addrule(m1,'drug_effect = 1 + (Imax * (pd.TE/oneTE))', ...
    'RuleType','RepeatedAssignment');

% adddose(m1,sbiodose('te_dose','Target','pd.TE','Amount', 0, ...
%     'AmountUnits','molecule', 'Interval', 24, ...
%     'TimeUnits','hour','RateUnits','molecule/hour', ...
%     'RepeatCount',0));

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



