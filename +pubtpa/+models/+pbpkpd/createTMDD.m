function m = createTMDD()
m = sbiomodel('tmdd');
tsCmp = m.addcompartment('t', 'Constant',true, 'Value', 3.1, ...
    'Units', 'liter', 'Notes', 'Peripheral tissue compartment');
plasmaCmp = m.addcompartment('p', 'Constant',true, 'Value', 3.06,...
    'Units', 'liter', 'Notes', 'Plasma compartment');
tsCmp.addspecies('D', 'Value', 0, 'Units', 'nanomole/liter', ...
    'Notes', 'Drug concentration in the peripheral tissue');
plasmaCmp.addspecies('D', 'Value', 0, 'Units', 'nanomole/liter', ...
    'Notes', 'Drug concentration in the plasma compartment');
plasmaCmp.addspecies('T', 'Value', 1, 'Units', 'nanomole/liter', ...
    'Notes', 'Target concentration in the plasma compartment');
plasmaCmp.addspecies('DT', 'Value',0, 'Units', 'nanomole/liter', ...
    'Notes', 'Drug-Target complex in the plasma compartment');

m.addparameter('Kpt', 'Value', 0.186, 'Units', '1/day');
m.addparameter('Ktp', 'Value', 0.184, 'Units', '1/day');
m.addparameter('KelT', 'Value', 1, 'Units', '1/day');
m.addparameter('KelD', 'Value', 0.085, 'Units', '1/day');
m.addparameter('KelDT', 'Value', 0.085/10, 'Units', '1/day');
m.addparameter('KsynT', 'Value', 1, 'Units', 'nanomole/(liter * day)');
m.addparameter('Kon', 'Value', 86.4, 'Units', 'liter / (nanomole * day)');
m.addparameter('Koff', 'Value', 1, 'Units', '1/day');
m.addparameter('T0', 'Value', 1, 'Units','nanomole/liter');
m.addparameter('t12T', 'Value', 1, 'Units', 'day');
m.addparameter('KD', 'Value', 100, 'Units', 'nanomole/liter');

m.addrule('KelT = 0.693/t12T', 'initialAssignment');
m.addrule('KsynT = KelT  * T0', 'initialAssignment');
m.addrule('Koff = KD * Kon', 'initialAssignment');


m.addreaction('null -> p.T', 'ReactionRate', 'KsynT', 'Name', 'RsnyT');
m.addreaction('p.T -> null', 'ReactionRate', 'KelT * p.T', 'Name', ...
    'RelT');
m.addreaction('p.D -> null', 'ReactionRate', 'KelD * p.D', 'Name', ...
    'RelDp');
m.addreaction('p.D <-> t.D', 'ReactionRate', 'Kpt * p.D - Kpt * t.D', ...
    'Name', 'RmigDpt');
m.addreaction('p.D + p.T <-> p.DT', 'ReactionRate', ...
    'Kon * p.D * p.T - Koff * p.DT', 'Name', 'RbinDT');
m.addreaction('p.DT -> null', 'ReactionRate', 'KelDT * DT', ...
    'Name', 'RelDT');

end