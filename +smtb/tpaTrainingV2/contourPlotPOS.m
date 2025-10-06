function f = contourPlotPOS(res_med, xlimit, ylimit)
% CONTOURPLOTPOS allows to generate a contour plot of probability of
% success with Peff_sol on xaxis and IC50_CLint on yaxis
res_med.peff_sol = res_med.Peff.*res_med.solubility;%absorption
res_med.ic50_clint = res_med.RC50.*res_med.CLint_u;%

med = res_med.MED_mg;

xvar = log10(res_med.ic50_clint);
yvar = log10(res_med.peff_sol);

n_pts = 20; % number of grid points
n_ctr = 15; % number of contours

xpoints = linspace(min(xvar),max(xvar),n_pts);
ypoints = linspace(min(yvar),max(yvar),n_pts);

xx = nan(1,n_pts^2);
yy = nan(1,n_pts^2);

dose = nan(1,n_pts*n_pts);
dose_prob = nan(1,n_pts*n_pts);
log_dose_median = nan(1,n_pts*n_pts);
log_dose_mean = nan(1,n_pts*n_pts);
log_dose_90prctile = nan(1,n_pts*n_pts);

er = log10(3);
count = 0;
for i = 1:length(xpoints)
    for j = 1:length(ypoints)
        count = count + 1;
        xx(count) = xpoints(i);
        yy(count) = ypoints(j);
        xidx = xvar > (xpoints(i)-er) & xvar < (xpoints(i)+er);
        yidx = yvar > (ypoints(j)-er) & yvar < (ypoints(j)+er);
        idx = xidx & yidx;
        if any(idx)
            dose(count) = mean(med(idx));
            dose_prob(count) = sum(med(idx) <= 500) / sum(idx);
            log_dose_median(count) = median(log10(med(idx)));
            log_dose_mean(count) = mean(log10(med(idx)));
            log_dose_90prctile(count) = prctile(log10(med(idx)),90);
        end
    end
end

[xgrid,ygrid] = meshgrid(10.^xpoints,10.^ypoints,n_pts);
z = griddata(10.^xx,10.^yy,dose_prob,xgrid,ygrid,'cubic');


f = figure('DefaultAxesFontSize', 16);
[~,h] = contourf(xgrid,ygrid,z,n_ctr);
set(h,'LineColor','none')
clim([0 1]);
xlim(xlimit);
ylim(ylimit);
colorbar(); colormap('jet');

xlabel('CLint (mL/min/g) x IC50 (ng/ml)')
ylabel('Solubility (mg/mL) x Peff (x10^-4cm/s)')%10^-7 papp
h.Parent.YScale = 'log';
h.Parent.XScale = 'log';
end