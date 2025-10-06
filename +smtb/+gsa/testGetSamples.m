clear; close all;
% parameterName, distType, distParameter, truncatingBounds
s = {...
    'logP', 'normal', {'mu', 2, 'sigma', 1}, {0,10}...
    'CLint', 'uniform', {'Lower', 0, 'Upper', 100},{} ...
    };


sampleobj = smtb.gsa.getSamples(methodType="structured",...
    N=2^3, paramDist=s);
sampleobj.generateSamples(40, 0);
sampleobj.generateSamples(100, 0);
sampleobj.generateSamples(300, 0);
disp(sampleobj);

f = figure('DefaultAxesFontSize', 14);
t = tiledlayout(f, 1, 1);
scatter(sampleobj.Samples{1}(:,1), ...
    sampleobj.Samples{1}(:,2),70, 'b', 'filled',...
    'DisplayName', 'Lvl1');
hold on;
scatter(sampleobj.Samples{2}(:,1), ...
    sampleobj.Samples{2}(:,2), 50, 'r', 'filled',...
    'DisplayName', 'Lvl2');
scatter(sampleobj.Samples{3}(:,1), ...
    sampleobj.Samples{3}(:,2), 30, 'g', 'filled',...
    'DisplayName', 'Lvl3');
scatter(sampleobj.Samples{4}(:,1),...
    sampleobj.Samples{4}(:,2),10,'k', 'filled',...
    'DisplayName', 'Lvl4');
xlabel(sampleobj.ParamNames{1}); ylabel(sampleobj.ParamNames{2}); legend()


%
% s(1).name = 'a';
% s(1).dist = 'uniform';
% s(1).parameters.lower = 0;
% s(1).parameters.upper = 6;
% 
% s(2).name = 'b';
% s(2).dist = 'uniform';
% s(2).parameters.lower = 0;
% s(2).parameters.upper = 2;
