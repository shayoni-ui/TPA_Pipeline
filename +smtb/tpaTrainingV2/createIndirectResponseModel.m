function m = createIndirectResponseModel()
% CREATEINDIRECTRESPONSE model file is an example file of a PD model
% that we may use for TPA. Here it creates simple indirect response model
% to be run for TPA analysis
% Please see the lines 13 - 31 below which includes few lines of code that 
% will be added to the your own PD model to be used for TPA analysis

m = sbiomodel('indirectResponse');

% add compartment
pdCmp = m.addcompartment('pd', 'Value', 1, 'Units', 'liter');

% TPA script requirement
% ---------------------------------------------------------------               
pdCmp.addspecies('TE', 'Value', 0, 'Units', 'molecule'); 
m.addparameter('one_TE', 'Value', 1, 'Units', 'molecule', 'Notes', ...
    'for non dimensionalizing TE');
m.addparameter('drugEffect', 'Value', 1, 'Units', 'dimensionless', ...
    'Constant', 0);
m.addparameter('Imax', 'Value', 1, 'Units', 'dimensionless');
m.addparameter('cl_TE', 'Value', 0, 'Units', '1/hour');
m.addreaction('pd.TE -> null', 'Name', 'TE clearance for TE Simulations', ...
    'ReactionRate', 'pd.TE * cl_TE');
m.addrule('drugEffect = 1 - Imax * pd.TE / one_TE ', 'RuleType', ...
    'RepeatedAssignment'); % Drug effect terms feedsback to modulate 
% one of the processes in the PD .. Please see line 49 below where in the 
% current  example it is used to modulate the Kin parameter in the model 
adddose(m,sbiodose('te_dose','Target','pd.TE','Amount', 0, ...
    'AmountUnits','molecule', 'Interval', 24, ...
    'TimeUnits','hour','RateUnits','molecule/hour','RepeatCount',0));
% ------------------------------------------------------------------              

% add species
pdCmp.addspecies('P', 'Value', 1, 'Units', 'nanomole/liter');

% add parameters
m.addparameter('kin', 'Value', 1, 'Units', 'nanomole/(liter * hour)', ...
    'Notes', 'synthesis rate P');
m.addparameter('kout', 'Value', 1, 'Units', '1/hour', 'Notes', ...
    'degradation rate of P');
m.addparameter('P0', 'Value', 1, 'Units', 'nanomole/liter', ...
    'Notes', 'Initial value of the protein P');
m.addparameter('Pinhb', 'Value', 0, 'Units', 'dimensionless', ...
    'Constant', 0, 'Notes', 'percent inhibition of protein P to plot');

% add reactions
m.addreaction('pd.P -> null', 'ReactionRate', 'kout * pd.P', ...
    'Name', 'R_degradation_P');
m.addreaction('null -> pd.P', 'ReactionRate', 'kin * drugEffect', ...
    'Name', 'R_synthesis_P');

% calculating percent inhibition of protein for plotting
m.addrule('Pinhb = 100 * (P0 - P) / P0', 'RuleType', ...
    'RepeatedAssignment', 'Notes', ['% inhibition of protein P compared' ...
    'to its initial value']);

% to make sure P and P0 values start at steady state 
% if you change kin and kout
% For complex models you will have to check if the PD model parameters
% are at steady state. Potentially add function
m.addrule('P0 = kin/kout', 'RuleType', 'initialAssignment'); 
m.addrule('P = P0', 'RuleType', 'initialAssignment');


cs = getconfigset(m);
cs.CompileOptions.DimensionalAnalysis = true;
cs.CompileOptions.UnitConversion = true;
cs.TimeUnits ='hour';

end
