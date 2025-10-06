function [m, mexp] =createMergedPBPKPD(pd, exportFlag)
arguments
    pd = 'indirectresponse';
    exportFlag = true;
end
switch lower(replace(pd, ' ' , ''))
    case 'indirectresponse'
        mPD = pubtpa.models.pd.indirectresponse;
    case 'precursoractivation'
        mPD = pubtpa.models.pd.precursoractivation;
    case 'tumor'
        mPD = pubtpa.models.pd.tumor;
    case 'tolerance'
        mPD = pubtpa.models.pd.tolerance;
    case 'transduction'
        mPD = pubtpa.models.pd.transduction;
    otherwise
        error('The entered PD model is not supported')
end

removeClTEequation(mPD);

[mPBPK, ~] = pubtpa.models.pbpk.createPBPK(false);
m = smtb.simbio.mrgMod(mPD, mPBPK);

m.addparameter('RC50', 'Value',1,'Units','nanogram/milliliter', ...
    'Notes', 'Potency', 'Constant', 1);
m.addrule(['pd.TE = (TargetTissue.drugU / ' ...
    '(RC50 + TargetTissue.drugU))*oneTE'],'RuleType', ...
    'RepeatedAssignment');

cs = m.configset;
cs.TimeUnits = 'hour';
cs.SolverOptions.AbsoluteTolerance = 1e-8;
cs.SolverOptions.RelativeTolerance = 1e-6;
cs.CompileOptions.UnitConversion = true;
cs.CompileOptions.DimensionalAnalysis = true;

verify(m);

if exportFlag
    mexp = export(m); accelerate(mexp);
else
    mexp = '';
    %sbioaccelerate(m);
end

end

function removeClTEequation(m)
    rxn = m.Reactions;
    params = m.Parameters;
    idxP = strcmpi({params.Name},'clTE');
    idxR = strcmpi({rxn.Name},'R_clTE');
    rxnTE = rxn(idxR);
    delete(rxnTE);
    paramTE = params(idxP);
    delete(paramTE);

end