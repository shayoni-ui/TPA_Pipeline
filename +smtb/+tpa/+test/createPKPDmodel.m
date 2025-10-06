function [m, mexp] = createPKPDmodel(pkdriver, exportFlag)
arguments
    pkdriver char; % name of the PK species that drives the PD model
    exportFlag = true;
end
% Started creation of the function. Implementation still left
% The idea is to replace the drug species added to the model with
% createPDmode that allows integration of drug with 
[mPD, ~] = createPDmodelwithMoa();
[mPK, ~] = createPKmodel(false);
m = smtb.simbio.mrgMod(mPK, mPD);

replaceDrugWithPKdriver(m, pkdriver);

cs = m.configset;
cs.TimeUnits = 'hour';
cs.SolverOptions.AbsoluteTolerance = 1e-8;
cs.SolverOptions.RelativeTolerance = 1e-6;
cs.CompileOptions.UnitConversion = true;
cs.CompileOptions.DimensionalAnalysis = true;

verify(m);

if exportFlag
    mexp = export(m); accelerate(mexp);
end

end

function replaceDrugWithPKdriver(m, pkdriver)
    disp(m);
    disp(pkdriver);
end