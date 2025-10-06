function [m, mexp] = createPDmodel(nva)
arguments
    nva.P0 = 0.075;
    nva.t12P = 0.33;
    nva.Imax = 1;
    nva.exportFlag = false;
end
% A simple indirect response model for protein P (CGRP)
m = sbiomodel('TRPM3model');
pd = m.addcompartment('pd', 'Value', 1, 'Units', 'milliliter', ...
    'constant', 1);
m.addparameter('P0', 'Value', nva.P0, 'Units', ...
    'nanogram/milliliter', 'Constant', 1);

pd.addspecies('P', 'Value', 0, 'Units', 'nanogram/milliliter');
m.addrule('P = P0', 'initialassignment');

m.addreaction('P -> null', 'ReactionRate', 'koutP * P',...
    'Name', 'R_degradation_P');
m.addreaction('null -> P', 'ReactionRate',...
    'kinP', 'Name', 'R_synthesis_P');

m.addparameter('t12P', 'Value', nva.t12P, 'Units', 'hour',...
    'constant', 1);
m.addparameter('koutP', 'Value', 1, 'Units', '1/hour', ...
    'constant', 1);
m.addparameter('log2', 'value', log(2), 'Units', 'dimensionless');
m.addrule('koutP = log2/t12P', 'initialassignment');

m.addparameter('kinP', 'value', 1, 'Units',...
    'nanogram / (milliliter * hour)');
m.addrule('kinP = P0 * koutP', 'initialassignment' );

cs = m.getconfigset();
cs.TimeUnits = 'hour';

verify(m);

if nva.exportFlag
    mexp = export(m); accelerate(mexp);
else
    mexp = SimBiology.export.Model.empty;
end
end