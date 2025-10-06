function m = createPBPKmodel()
% generate the PBPK model structure from the template file in the toolbox
m = BuildPBPKStructure();


% You can use the unbound plasma or unbound concentration from the target
% tissue of interes from the PBPK to drive the model
% Optionally you can expand the PBPK model with a novel compartment that
% dosenot already exist. Check the TRPM2/TRPM3 TPA code to see example of 
% addition of 4 compartment Brain model to expand our standard PBPK model 

targetCmp = addcompartment(m, 'TargetTissue', 'Capacity',1,...
    'CapacityUnits','milliliter'); % example
addspecies(targetCmp, 'drug', 'InitialAmount',0, ...
    'InitialAmountUnits','nanogram/milliliter'); % total drug concentration
addspecies(targetCmp, 'drugU', 'InitialAmount',0, ...
    'InitialAmountUnits','nanogram/milliliter'); % unbound drug conc

addparameter(m,'Q_TargetTissue','Value',100,'ValueUnits', ...
    'milliliter/hour'); % Blood flow to target tissue. Really high here

addreaction(m, 'null -> TargetTissue.drug', 'ReactionRate', ...
    'Q_TargetTissue*blood.drug');
addreaction(m, 'TargetTissue.drug -> null', 'ReactionRate', ...
    'Q_TargetTissue * B2P/Kp_Spleen * TargetTissue.drug');
% B2P is blood to plasma ratio and using Kp of spleen as a surrogate for
% the target tissue. Create a new Kp parameter if you know the value of 
% Kp for the target tissue

addrule(m, 'TargetTissue.drugU = TargetTissue.drug*fup/Kp_Spleen', ...
    'RuleType', 'RepeatedAssignment');


cs = getconfigset(m);
cs.CompileOptions.DimensionalAnalysis = true;
cs.CompileOptions.UnitConversion = true;
cs.TimeUnits = 'hour';
verify(m,cs);

end