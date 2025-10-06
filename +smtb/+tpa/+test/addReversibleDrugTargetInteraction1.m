function addReversibleDrugTargetInteraction(mPD, rnName)
arguments
    mPD SimBiology.Model;
    rnName char;
end

%rnName = 'R_synthesis_P';
cmp = mPD.addcompartment('drugCmp', 'capacity', 1, 'Units', 'milliliter');
cmp.addspecies('drug', 'Value',0, 'Units', 'milligram/milliliter');

mPD.addparameter('TE', 'Value', 0, 'Constant',false,...
    'Units','dimensionless');
mPD.addparameter('drugEffect', 'Value', 1, 'Constant',false, ...
    'Units', 'dimensionless');
mPD.addparameter('Imax', 'Value', 1, 'Units','dimensionless');
mPD.addparameter('XC50', 'Value', 1, 'Units', 'milligram/milliliter')

mPD.addrule('drugEffect = 1 - Imax * TE', ...
    'RuleType', 'RepeatedAssignment');
mPD.addrule('TE = drugCmp.drug / (drugCmp.drug + XC50)', ...
    'RuleType', 'RepeatedAssignment');

mPD.addparameter('kdegDrug', 'Value', 1e3, 'Units', '1/hour');
mPD.addreaction('drugCmp.drug -> null', 'ReactionRate',...
    'kdegDrug * drugCmp.drug');

duration = 20; TEamplitude = 0;
doseAmplitude = TEamplitude / (1-TEamplitude);
doseAmount = doseAmplitude * duration;
doseTE = sbiodose('doseTE', 'repeat', 'Amount', doseAmount*1e3,...
    'TargetName', 'drugCmp.drug', 'RepeatCount',20,...
    'Rate', doseAmplitude*1e3, 'RateUnits', 'milligram/hour',...
    'AmountUnits','milligram', 'TimeUnits', 'hour', 'Interval', 24);

mPD.adddose(doseTE);

rn = sbioselect(mPD,'Type', 'Reaction', 'name', rnName);
rn.ReactionRate = sprintf('(%s) * drugEffect', rn.ReactionRate);
end


% TEsamples = linspace(0.01, 0.99, 10000);
% doseSamples = TEsamples./(TEsamples + 1);
% aa = logspace(-2, 2, 10000); bb = aa./(aa + 1);
% 
% plot(aa,bb , 'ko'); set(gca, 'Xscale', 'log');
% figure; histogram(bb, 100);