%% Creating PD model
[m, ~] = smtb.tpa.test.createPDmodel();

%% perform mode1Analysis
npulses = 20; shrinkData = 0; 
[sd, TeDuration] = smtb.tpa.performMode1Analysis(m, 'R_synthesis_P', ...
    ["TE", "P"], shrinkData, npulses);

%% Plotting 
plotProfiles(sd(3))
response = getResponseVar(sd);
contourPdurationTE(TeDuration(:,1), TeDuration(:,2), response)



%% Helper functions

function contourPdurationTE(x,y,z)
    xgrid = reshape(x*100, sqrt(length(x)), sqrt(length(x))); 
    ygrid = reshape(y, sqrt(length(x)), sqrt(length(x)));
    zgrid = reshape(z, sqrt(length(x)), sqrt(length(x)));

    f = figure('DefaultAxesFontSize', 12);
    contourf(xgrid, ygrid, zgrid,15, ...
        'LineStyle', 'none'); 
    %caxis([0, 100]);
    colorbar();
    colormap(jet);
    %set(gca, 'Xscale', 'log');
    %set(gca, 'Yscale', 'log');
    hold on;
    contour(xgrid, ygrid, zgrid,[34, 34], ...
        'color','k', 'LineWidth',2);
    %hold on;
    %scatter(x, y, 20, 'k', 'filled')

    xlabel('% TE')
    ylabel('Duration , hour')%10^-7 papp


    %xlim([0 1])
    %ylim([0 24])
    exportgraphics(f, 'mode1TeDuration.png', 'resolution', 300);
end

function response = getResponseVar(sd)
    %tlast = d.Interval * (d.RepeatCount - 1);
    
    response = zeros(1, length(sd));
    for ii = 1:length(sd)
        cgrp = selectbyname(sd{ii}, 'P');
        response(ii) = 1000*cgrp.Data(end);
    end
end

function plotProfiles(sd, fname)
arguments
    sd;
    fname=false;
end
f = figure('DefaultAxesFontSize', 12);
tiledlayout('flow')
for ii = 1:length(sd{1}.DataNames)
    nexttile;
    for jj = 1:length(sd)
        if strcmpi(sd{jj}.DataNames{ii}, 'CGRP')
            scale = 1000;
            units = 'pg/mL';
        else
            scale = 1;
            units = sd{jj}.DataInfo{ii}.Units;
        end
        plot(sd{jj}.Time, sd{jj}.Data(:, ii)*scale, 'LineWidth',2); hold on;
    end
    xlabel(sprintf('Time, %s', sd{jj}.TimeUnits));
    ylabel(sprintf('%s, %s',sd{jj}.DataInfo{ii}.Name, units));
    xlim([0, 480]);

end

if fname
    exporgraphics(f, fname, 'resolution', 300);
end
end


