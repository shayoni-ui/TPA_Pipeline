function [m, mexp] = createPBPK(exportFlag)
arguments
exportFlag = true;
end

m = smtb.pbpk.buildPBPKmodelTemplate();
tt = addcompartment(m, 'TargetTissue', 'Capacity',1,...
    'CapacityUnits','milliliter');%brain blood
addspecies(tt, 'drug', 'InitialAmount',0, ...
    'InitialAmountUnits','nanogram/milliliter');
addspecies(tt, 'drugU', 'InitialAmount',0, ...
    'InitialAmountUnits','nanogram/milliliter');

addparameter(m,'Q_TargetTissue','Value',100,'ValueUnits', ...
    'milliliter/hour');

addreaction(m, 'null -> TargetTissue.drug', 'ReactionRate', ...
    'Q_TargetTissue*blood.drug');
addreaction(m, 'TargetTissue.drug -> null', 'ReactionRate', ...
    'Q_TargetTissue * B2P/Kp_Spleen * TargetTissue.drug');

addrule(m, 'TargetTissue.drugU = TargetTissue.drug*fup/Kp_Spleen', ...
    'RuleType', 'RepeatedAssignment');


cs = getconfigset(m);
cs.CompileOptions.DimensionalAnalysis = true;
cs.CompileOptions.UnitConversion = true;
cs.TimeUnits = 'hour';
verify(m,cs);

if exportFlag
mexp = export(m);
%accelerate(mexp);
else
    mexp = '';
end

end