function mPKPD = integratePDandPBPK(mPD, mPK)
    arguments
        mPD SimBiology.Model;
        mPK SimBiology.Model;
    end
% integratePDandPBPK allows to merge the PD and PBPK model for TPA analysis
% Inputs:
%   mPD - PD model developed by the user;
%   mPK - PK model developed using the template in the toolbox and
%   potentially expanded to include new compartments
% Outputs:
%   mPKPD - integrade PBPK-PD model
% Note: the srcipt currently supports only small molecule PBPK model 
%   but plans to support different modes of action in future


% Remove the cl_TE equation introducted in the PD model for mode 1 analysis
removeClTEequation(mPD);

% merge PBPK and PD model
mPKPD = mrgMod(mPK, mPD);

% Define a rule describing the species of the PBPK model driving the PD 

mPKPD.addparameter('RC50', 'Value', 1, 'Units', 'nanogram/milliliter', ...
    'Notes', 'Potency parameter');
mPKPD.addrule(['pd.TE = (TargetTissue.drugU ' ...
    '/ (RC50 + TargetTissue.drugU)) * one_TE'], ...
    'RuleType', 'RepeatedAssignment');


cs = mPKPD.configset;
cs.TimeUnits = 'hour';
cs.SolverOptions.AbsoluteTolerance = 1e-8;
cs.SolverOptions.RelativeTolerance = 1e-6;
cs.CompileOptions.UnitConversion = true;
cs.CompileOptions.DimensionalAnalysis = true;

verify(mPKPD);

end


function removeClTEequation(mPD)
    rxn = mPD.Reactions;
    params = mPD.Parameters;
    idxP = strcmpi({params.Name},'cl_TE');
    idxR = strcmpi({rxn.Name},'TE clearance for TE Simulations');
    rxnTE = rxn(idxR);
    delete(rxnTE);
    paramTE = params(idxP);
    delete(paramTE);

end