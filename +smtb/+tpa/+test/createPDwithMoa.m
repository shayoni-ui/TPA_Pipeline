function [m, mexp] = createPDwithMoa()
[m, ~] = smtb.tpa.test.createPDmodel();
smtb.tpa.includeMoaPD(m, 'Reversible',...
    'targetProcesses', 'R_synthesis_P', 'functionalForms','Reversible');
cs = getconfigset(m);
cs.TimeUnits = 'hour'; cs.StopTime = 24*20;
cs.CompileOptions.UnitConversion = true;
cs.SolverOptions.AbsoluteTolerance = 1e-8;
cs.SolverOptions.RelativeTolerance = 1e-7;

mexp = export(m); accelerate(mexp);

end

