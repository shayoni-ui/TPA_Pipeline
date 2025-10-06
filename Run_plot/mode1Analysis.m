% Perform mode 1 Analysis
% Look at mode1AnalysisRun function below for help

%% 1. Indirect Response model with fast turnover rate
v = sbiovariant('a');
v.addcontent({'Parameter', 'ksyn_P', 'Value', 1});
v.addcontent({'Parameter', 'kdeg_P', 'Value', 1});
mode1AnalysisRun("indirectresponse", "indirectfastTEduration.png", ...
    50, v, 'Indirect response (fast turnover)');



%% 2. Indirect response model with slow turnover rate
v = sbiovariant('a');
v.addcontent({'Parameter', 'ksyn_P', 'Value', 0.01});
v.addcontent({'Parameter', 'kdeg_P', 'Value', 0.01});
mode1AnalysisRun("indirectresponse", "indirectslowTEduration.png", 50, v, ...
    'Indirect response (slow turnover)');

%% 
% 3. Precursor activation model
mode1AnalysisRun("precursoractivation", "precursorTEduration.png", 1.1, {}, ...
    'Precursor activation');

%% 4. Transduction model
mode1AnalysisRun("transduction", "transductionTEduration.png", 1.3, {}, ...
    'Transduction');

%% 5. Tolerance model
mode1AnalysisRun("tolerance", "toleranceTEduration.png", 1.4, {}, ...
    'Tolerance');

%% 6. Tumor growth model
mode1AnalysisRun("tumor", "tumorTEduration.png", 30, {}, ...
    'Tumor growth');


%%
function mode1AnalysisRun(pd, fname, lvl, v, titletext)
switch pd % select the PD model
    case 'indirectresponse'
        m = pubtpa.models.pd.indirectresponse;
        responsename = "P"; scale = 1; 
    case 'precursoractivation'
        m = pubtpa.models.pd.precursoractivation;
        responsename = "A"; scale = 1;
    case 'transduction'
        m = pubtpa.models.pd.transduction;
        responsename = "M1"; scale = 1;
    case 'tolerance'
        m = pubtpa.models.pd.tolerance;
        responsename = "P1"; scale = 1;
    case 'tumor'
        m = pubtpa.models.pd.tumor;
        responsename = "TGI"; scale = 1;
    otherwise
        error('The entered model is not supported or defined')
end

% Add value for the TE degradation parameter and dose for 
% to simulate step functions
clTE= sbioselect(m, 'Name', 'clTE'); clTE.value = 1e3;
m.adddose(sbiodose('TEdose','Target','TE','Amount', 0, ...
    'AmountUnits','molecule', ...
    'TimeUnits','hour','RateUnits','molecule/hour', ...
    'RepeatCount',0));
cs = m.getconfigset; cs.TimeUnits = 'hour';
cs.CompileOptions.UnitConversion = true;
mexp = export(m); %accelerate(mexp);
d = mexp.getdose;
d.Interval = 24; % hours
d.RepeatCount = 20; % 
mexp.SimulationOptions.StopTime = d.Interval * d.RepeatCount;

%% Generate TE vs Duration parameters /virtual compounds
N = 20;
TeDuration = combvec(linspace(0, 1, N), linspace(1, 24, N));

switch pd
    case 'indirectresponse'
        smtb.simbio.applyVariantContent(mexp, v);
end
%% Perform simulation of the PD model at different TE duration and 
% amplitude. Save the results of only last dose using shrinkData = 1
% The performVirtualExploration function is part of the SMTB package
% written by Jaimit and not our SMTBToolbox used for the tutorial
tic
sd = smtb.tpa.performVirtualExploration(mexp, d,'TeDuration',...
    TeDuration, 'selectNames', ["TE", responsename],...
    'shrinkData',1);
toc

%% Plotting
%plotProfiles(sd(end))

%% Extracting biomarkers
response = getResponseVar(sd, responsename, scale);

switch pd
    case 'indirectresponse'
        response = (1-response)*100;
    otherwise
        fprintf('response left as is \n');
end
%% Contour plot
contourdurationTE(TeDuration(:,1), TeDuration(:,2), response, fname, lvl, ...
    titletext, N)

end

%% Helper functions
function contourdurationTE(x,y,z, fname, lvl, titletext, N)
xgrid = reshape(x*100, N, N); ygrid = reshape(y, N, N);
zgrid = reshape(z, N, N);

f = figure('DefaultAxesFontSize', 12);
contourf(xgrid, ygrid, zgrid,15, ...
    'LineStyle', 'none');
clim([min(zgrid(:)) max(zgrid(:))]);
colormap(jet); colorbar();
hold on;
contour(xgrid, ygrid, zgrid,[lvl, lvl], ...
    'color','k', 'LineWidth',2)
xlabel('% TE')
ylabel('Duration , hour')%10^-7 papp
title(titletext)
fname = strcat('./pngFiles/', fname);
exportgraphics(f, fname, 'resolution', 300);
end

function response = getResponseVar(sd, name, scale)
response = zeros(1, length(sd));
for ii = 1:length(sd)
    rr = selectbyname(sd{ii}, name);
    response(ii) = scale*max(rr.Data);
end
end

% function plotProfiles(sd, fname)
% arguments
%     sd;
%     fname=false;
% end
% f = figure('DefaultAxesFontSize', 12);
% tiledlayout('flow')
% for ii = 1:length(sd{1}.DataNames)
%     nexttile;
%     for jj = 1:length(sd)
%         if strcmpi(sd{jj}.DataNames{ii}, 'CGRP')
%             scale = 1000;
%             units = 'pg/mL';
%         else
%             scale = 1;
%             units = sd{jj}.DataInfo{ii}.Units;
%         end
%         plot(sd{jj}.Time, sd{jj}.Data(:, ii)*scale, 'LineWidth',2); hold on;
%     end
%     xlabel(sprintf('Time, %s', sd{jj}.TimeUnits));
%     ylabel(sprintf('%s, %s',sd{jj}.DataInfo{ii}.Name, units));
%     xlim([0, 480]);
% 
% end
% 
% if fname
%     exporgraphics(f, fname, 'resolution', 300);
% end
% end
