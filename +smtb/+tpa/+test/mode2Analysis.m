disp('Starting parameterization to get dervied parameters...')
tic
N = 500;
sampleSetPK = smtb.gsa.getSamples('methodtype', 'lhs', 'N', N, ...
    'paramDist', ...
        {'CLint', 'LogUniform', {'Lower', 0.1, 'Upper', 100}, {}...
         'logP', 'uniform', {'Lower', 0, 'Upper', 5},         {}} ...
    );
%
% if debug
%     sampleSetPK.visualizeDist;
% end


[variantArr, qoiArr] = smtb.tpa.generateVCs(sampleSetPK);

for ii = 1:length(variantArr)
    variantArr{ii} = smtb.simbio.updateVariantContent(variantArr{ii},...
        'fup', fup);
end
toc
disp('Finished parameterization')


%% Perform the simulation using the model and the generated virtual
% compound variants
disp('Starting simulation of virtual compounds')
tic
mexp = trpm3Model.mexp;

d = mexp.getdose('PO');
d.Amount = 100; d.Interval = 24; d.RepeatCount = 20;
mexp.SimulationOptions.StopTime = 480;

sd=smtb.tpa.performVirtualExploration(mexp, d,'variantArr',variantArr, ...
    'selectNames',["plasma.drug", "TargetTissue.drugU", ...
    "pd.CGRP"], "shrinkData",true);
toc
disp('Finished simulation of virtual compounds');