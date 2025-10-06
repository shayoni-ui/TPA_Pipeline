%% Plotting mode 2 analysis. Look at the function mode2Plots below for 
% more details

% 1. Indirect response fast turover rate

indirectfastRes = load('./matFiles/indirectFast01.mat');
mode2Plots("indirectresponse",indirectfastRes.results, ...
    'CLintKpIndirectFast01.png',...
    'CminCavgIndirectFast01.png', ...
     'CminCavgScatterIndirectFast01.png', ...
    'CmintRC50ScatterIndirectFast01.png', ...
    'featuresIndirectFast01.csv', "P", "max",...
    [0.15, 10], [0.15, 10], [0, 100], [0.1, 10], [0.1, 10], [0, 100], 50, ...
    'Indirect response (fast turnover), fub = 0.01');
%% 


indirectfastRes = load('./matFiles/indirectFast05.mat');
mode2Plots("indirectresponse",indirectfastRes.results, ...
    'CLintKpIndirectFast05.png',...
    'CminCavgIndirectFast05.png', ...
    'CminCavgScatterIndirectFast05.png', ...
    'CmintRC50ScatterIndirectFast05.png', ...
    'featuresIndirectFast05.csv', "P", "max",...
    [0.15, 10], [0.15, 10], [0, 100], [0.1, 10], [0.1, 10], [0, 100], 50, ...
    'Indirect response (fast turnover), fub = 0.05');


indirectfastRes = load('./matFiles/indirectFast1.mat');
mode2Plots("indirectresponse",indirectfastRes.results, ...
    'CLintKpIndirectFast1.png',...
    'CminCavgIndirectFast1.png', ...
    'CminCavgScatterIndirectFast1.png', ...
    'CmintRC50ScatterIndirectFast1.png', ...
    'featuresIndirectFast1.csv', "P", "max",...
    [0.15, 10], [0.15, 10], [0, 100], [0.1, 10], [0.1, 10], [0, 100], 50, ...
    'Indirect response (fast turnover), fub = 0.1');

indirectfastRes = load('./matFiles/indirectFast2.mat');
mode2Plots("indirectresponse",indirectfastRes.results, ...
    'CLintKpIndirectFast2.png',...
    'CminCavgIndirectFast2.png', ...
    'CminCavgScatterIndirectFast2.png', ...
    'CmintRC50ScatterIndirectFast2.png', ...
    'featuresIndirectFast2.csv', "P", "max",...
    [0.15, 10], [0.15, 10], [0, 100], [0.1, 10], [0.1, 10], [0, 100], 50, ...
    'Indirect response (fast turnover), fub = 0.2');

%%
% 2. Indirect response slow turnover rate


indirectslowRes = load('./matFiles/indirectSlow01.mat');
mode2Plots("indirectresponse",indirectslowRes.results, ...
    'CLintKpIndirectSlow01.png',...
    'CminCavgIndirectSlow01.png', ...
    'CminCavgScatterIndirectSlow01.png', ...
    'CmintRC50ScatterIndirectSlow01.png', ......
    'featuresIndirectSlow01.csv', "P", "max",...
    [0.15, 10], [0.15, 10], [0, 100], [0.1, 10], [0.1, 10], [0, 100], 50, ...
    'Indirect response (slow turnover), fub = 0.01');



indirectslowRes = load('./matFiles/indirectSlow05.mat');
mode2Plots("indirectresponse",indirectslowRes.results, ...
    'CLintKpIndirectSlow05.png',...
    'CminCavgIndirectSlow05.png', ...
    'CminCavgScatterIndirectSlow05.png', ...
    'CmintRC50ScatterIndirectSlow05.png', ......
    'featuresIndirectSlow05.csv', "P", "max",...
    [0.15, 10], [0.15, 10], [0, 100], [0.1, 10], [0.1, 10], [0, 100], 50, ...
    'Indirect response (slow turnover), fub = 0.05');

indirectslowRes = load('./matFiles/indirectSlow1.mat');
mode2Plots("indirectresponse",indirectslowRes.results, ...
    'CLintKpIndirectSlow1.png',...
    'CminCavgIndirectSlow1.png', ...
    'CminCavgScatterIndirectSlow1.png', ...
    'CmintRC50ScatterIndirectSlow1.png', ......
    'featuresIndirectSlow1.csv', "P", "max",...
    [0.15, 10], [0.15, 10], [0, 100], [0.1, 10], [0.1, 10], [0, 100], 50, ...
    'Indirect response (slow turnover), fub = 0.1');

indirectslowRes = load('./matFiles/indirectSlow2.mat');
mode2Plots("indirectresponse",indirectslowRes.results, ...
    'CLintKpIndirectSlow2.png',...
    'CminCavgIndirectSlow2.png', ...
    'CminCavgScatterIndirectSlow2.png', ...
    'CmintRC50ScatterIndirectSlow2.png', ...
    'featuresIndirectSlow2.csv', "P", "max",...
    [0.15, 10], [0.15, 10], [0, 100], [0.1, 10], [0.1, 10], [0, 100], 50, ...
    'Indirect response (slow turnover), fub = 0.2');
%%
% 3. Precursor activation model
precursorRes = load('./matFiles/precursor.mat');
mode2Plots("precursoractivation",precursorRes.results, ...
    'CLintKpPrecursor.png',...
    'CminCavgPrecursor.png','CminCavgScatterPrecursor.png', ...
    'CmintRC50ScatterPrecursor.png', ...
    'featuresPrecursor.csv', "A", 'max',...
    [0.15, 10], [0.15, 10], [1, 1.4], [0.1, 10], [0.1, 10], [1, 1.4], 1.3, ...
    'Precursor activation');

%%
% 4. Transduction model
transductionRes = load('./matFiles/transduction.mat');
mode2Plots("transduction",transductionRes.results, ...
    'CLintKpTransduction.png',...
    'CminCavgTransduction.png','CminCavgScatterTransduction.png', ...
    'CmintRC50ScatterTransduction.png', ...
    'featuresTransduction.csv', "M1", "max",...
    [0.15, 10], [0.15, 10], [1, 1.8], [0.1, 10], [0.1, 10], [1, 1.8], 1.45, ...
    'Transduction');

%%
% 5. Tolerance model
toleranceRes = load('./matFiles/tolerance.mat');
mode2Plots("tolerance",toleranceRes.results, ...
    'CLintKpTolerance.png',...
    'CminCavgTolerance.png','CminCavgScatterTolerance.png', ...
    'CmintRC50ScatterTolerance.png', 'featuresTolerance.csv', ...
    "P1", "max", ...
    [0.15, 10], [0.15, 10], [1, 1.7], [0.1, 10], [0.1, 10], [1, 1.7], 1.5, ...
    'Tolerance');

%%
% 6. Tumor growth model
tumorRes = load('./matFiles/tumor.mat');
mode2Plots("tumor",tumorRes.results, ...
    'CLintKpTumor.png',...
    'CminCavgTumor.png', 'CminCavgScatterTumor.png', ...
    'CmintRC50ScatterTumor.png', ...
    'featuresTumor.csv', "TGI", "max",...
    [0.15, 10], [0.15, 10], [0, 50], [0.1, 10], [0.1, 10], [0, 50], 30, ...
    'Tumor growth');



%% Function to plot cmin vs cavg  and XC50 x CLint vs Kp

function mode2Plots(pd,results, fname1, fname2, fname3,fname4, fname5, ...
    responsename, responsetype, xlim1, ylim1, crange1, ...
    xlim2, ylim2, crange2, lvl, titletext)

%% Extract fup from the first virtual compound variant. Assuming
% fup to be constant 
vTb= smtb.simbio.variantContentTable(results.variantArr{1});
fup = vTb(string(vTb.Name) == "fup", 'Value').Variables;

out = getResponseVar(pd,results.sd, results.dose, results.qoiArr, 1,...
    responsename, responsetype, fup);
out.CLint = results.sampleSetPK.Samples{1}(:,1);
%%
contourXY(out.CLint, out.kpAvg', out.response', ...
    'CLint x XC50', 'Kp', xlim1, ylim1, crange1, ...
    fname1, lvl,titletext);

contourXY(arrayfun(@(x)x.cmin, out.ncaO)',...
    arrayfun(@(x)x.cavg, out.ncaO)', ...
    out.response', ...
    'Cmin / XC50', 'Cavg / XC50', xlim2, ylim2, crange2, ...
    fname2, lvl, titletext);

contourXY(out.CLint, fup'./out.kpAvg', out.response', ...
    'CLint x XC50', 'fut', xlim1, flip(fup./ylim1), crange1, ...
    strcat('fut_',fname1), lvl,titletext);


f = figure('DefaultAxesFontSize', 12); 
scatter((arrayfun(@(x)x.cmin, out.ncaO)'), ...
    (arrayfun(@(x)x.cavg, out.ncaO)'), 15,out.response', 'filled','o');
colormap("jet"); 
xlim(xlim2); ylim(ylim2);
set(gca, 'Xscale', 'log'); set(gca, 'Yscale', 'log')
xlabel('Cmin/RC50'); ylabel('Cavg/RC50'); colorbar();
fname3 = strcat('./pngFiles/', fname3);
exportgraphics(f, fname3)

f = figure('DefaultAxesFontSize', 12); 
scatter((arrayfun(@(x)x.cmin, out.ncaO)'), ...
    (arrayfun(@(x)x.tRC50, out.ncaO)'), 15,out.response', 'filled','o');
colormap("jet"); 
xlim(xlim2); ylim([0, 24]);
set(gca, 'Xscale', 'log'); %set('gca', 'Yscale', 'log')
xlabel('Cmin/RC50'); ylabel('t > RC50'); colorbar();
fname4 = strcat('./pngFiles/', fname4);
exportgraphics(f, fname4)


features=struct2table(out.ncaO);
features.response = out.response';
fname5 = strcat('./csvFiles/', fname5);
writetable(features,fname5);
end

% Helper functions

function contourXY(x,y,z, xlab, ylab, xl, yl, crange, fname, lvl,...
    titletext)
arguments
    x,
    y,
    z,
    xlab = 'CLint x XC50';
    ylab = 'Kp';
    xl = [];
    yl = [];
    crange = []
    fname = 'temp.png';
    lvl = 50;
    titletext = '';
end
    n_pts = 50;
    xpoints = linspace(log10(min(x)),log10(max(x)),n_pts);
    ypoints = linspace(log10(min(y)),log10(max(y)),n_pts);

    [xgrid,ygrid] = meshgrid(xpoints,ypoints,n_pts);
    zgrid = griddata(log10(x),log10(y),z,xgrid,ygrid,'cubic');


    f = figure('DefaultAxesFontSize', 12);
    contourf(10.^xgrid, 10.^ygrid, zgrid,15, ...
        'LineStyle', 'none'); 
    if ~isempty(crange); clim(crange); end
    colorbar();
    colormap(jet);
    %colormap(flipud(magma));
    hold on;
    contour(10.^xgrid, 10.^ygrid, zgrid,[lvl lvl], ...
        'color','k', 'LineWidth',2);

    xlabel(xlab)
    ylabel(ylab)%10^-7 papp
    title(titletext)

    if ~isempty(xl); xlim(xl); end
    if ~isempty(yl); ylim(yl); end
    set(gca, 'Xscale', 'log');
    set(gca, 'Yscale', 'log');

    fname = strcat('./pngFiles/', fname);
    exportgraphics(f, fname, 'resolution', 300);
end

function out = getResponseVar(pd, sd, dose, qoiArr, RC50, responsename, ...
    responsetype, fup)
    %tlast = dose.Interval * (dose.RepeatCount - 1);
    
    response = zeros(1, length(sd));
    %ncaO= cell(1, length(sd));
    kpAvg = zeros(1, length(sd));
    idxR = find(strcmpi(sd{1}.DataNames, responsename));
    for ii = 1:length(sd)
        switch responsetype
            case 'last'
                response(ii) = sd{ii}.Data(end,idxR);
            case 'max'
                response(ii) = max(sd{ii}.Data(:, idxR));
            case 'min'
                response(ii) = min(sd{ii}.Data(:, idxR));
            case 'avg'
                response(ii) = mean(sd{ii}.Data(:, idxR));
        end
        %% extract plasma concentration to determine
        % cmin, cmax, cavg, tRC50, tmax, auc etc.
        drugConc = selectbyname(sd{ii}, 'plasma.drug');
        time = drugConc.Time;
        ncaO(ii) = getDrugSummary(time, fup*drugConc.Data,dose, RC50);
        kpAvg(ii)=qoiArr{ii}.kpAvg;
    end

    switch pd
        case "indirectresponse"
            response = (1-response)*100;
        otherwise
            fprintf('response left as it is \n');
    end
    
    out.response = response;
    out.ncaO = ncaO;
    out.kpAvg = kpAvg;
end

function ncaO = getDrugSummary(t,c,d, RC50)
c = c(t > d.Interval * (d.RepeatCount - 1));
t = t(t > d.Interval * (d.RepeatCount - 1));

time = t;
conc = c;

[maxV, maxI] = max(conc);
ncaO.RC50 = RC50;
ncaO.cmax = maxV; % maximum conc
ncaO.cavg = trapz(time, conc) / (time(end) - time(1)); % avg conc
ncaO.tmax = time(maxI) - time(1); % time of max conc
ncaO.cmin =  conc(find(~isnan(conc), 1, 'last')); % trough conc
ncaO.tmin = time(find(~isnan(conc), 1, 'last')) - time(1); % time at trough
[ncaO.auc, ncaO.aucm] = trapzlog(time, conc); % AUC

ncaO.cmax = ncaO.cmax / RC50;
ncaO.cmin = ncaO.cmin / RC50;
ncaO.cavg = ncaO.cavg / RC50;

% Estimate time the conc remains above RC50
if ~isempty(RC50)
    idx2 = [];
    idx = find(conc > RC50);
    if ~isempty(idx)
        idx1 = idx(1);
        if length(idx) > 1
            idx2 = idx(end);
        else
            ncaO.tRC50 = 0;
        end
        
        if ~isempty(idx1) && ~isempty(idx2)
            ncaO.tRC50 = time(idx2) - time(idx1);
        end
    else
        ncaO.tRC50 = 0;
    end

end

end