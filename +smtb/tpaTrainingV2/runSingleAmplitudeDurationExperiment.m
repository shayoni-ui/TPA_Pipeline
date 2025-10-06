function sdata = runSingleAmplitudeDurationExperiment(mexp, teAmount, ...
    duration, stopTimeDays, pdVariantFileName, plotVariables, ylimit)
arguments
    mexp SimBiology.export.Model;
    teAmount double;
    duration double;
    stopTimeDays double;
    pdVariantFileName char;
    plotVariables cell = {}; 
    ylimit cell = {};

end
exportedModel = mexp;
kdeg_te = 1e6;
variantContent = buildPDVariantContent(pdVariantFileName);
variantContent = addVariantContent(variantContent, ...
    {'parameter', 'cl_TE', 'Value', kdeg_te});
% Remove unused species, parameters, etc. from variant content.
vcTable = variantContentTable(variantContent{2});
mNames = {exportedModel.ValueInfo.Name};
namesToRemove = setdiff(vcTable.name, mNames);
for i=1:length(namesToRemove)
    variantContent{2} = removeVariantContent( ...
        variantContent{2}, namesToRemove{i});
end

exportedModel = applyVariant(exportedModel, variantContent{2});

doseObj = exportedModel.getdose('te_dose');
doseObj.Rate = kdeg_te*teAmount;
doseObj.Amount = doseObj.Rate * duration;
doseObj.RepeatCount = stopTimeDays - 1;
stopTimeHours = 24*stopTimeDays;
exportedModel.SimulationOptions.StopTime = stopTimeHours;

sdata = simulate(exportedModel, doseObj);


if ~isempty(plotVariables)
    figure('DefaultAxesFontSize', 16);
    tiledlayout('flow');
    plotData = selectbyname(sdata, plotVariables);
    for ii = 1:length(plotVariables)
        nexttile;
        plot(plotData.Time, plotData.Data(:, ii), 'LineWidth', 2, ...
            'Color','k');
        xlabel('Time, hr');
        ylabel(plotData.DataNames(ii));
        xlim([0, stopTimeHours]);
        if ~isempty(ylimit)
        ylim(ylimit{ii});
        end
        
    end
end


end