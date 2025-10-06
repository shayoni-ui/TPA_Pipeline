m = pubtpa.models.pbpkpd.createTMDD;
cs = m.getconfigset();
cs.StopTime = 7*10; % days
cs.TimeUnits = 'day';

d1 = sbiodose('d1','repeat');
%d1.Amount = dosenm;
d1.AmountUnits = 'nanomole';
d1.TargetName = 'p.D';
d1.Interval = 7;
d1.TimeUnits = 'day';
d1.RepeatCount = 0;
dose = logspace(-1, 2, 30);

for jj = 1:length(dose)
mw=100e3; % Dalton
dosemg = dose(jj); % mg
dosenm = 1e9*dosemg*1e-3/mw; % nanomole
d1.Amount = dosenm;
% p1 = sbioselect(m, 'Type', 'Parameter', 'Name', 'Kpt');
% p1.Value = 0.186;
% p2 = sbioselect(m, 'Type', 'Parameter', 'Name', 'Ktp');
% p2.Value = 0.184;
% p3 = sbioselect(m, 'Type', 'Parameter', 'Name', 'KelD');
% p3.Value = 0.085;
p4 = sbioselect(m, 'Type', 'Parameter', 'Name', 'KD');
p4.Value = 0.01;
% 
% p = sbioselect(m, 'Type', 'Parameter', 'Name', 't12T');
% p.Value = 0.1;
sd(jj) = sbiosimulate(m, d1);
end




f = figure('DefaultAxesFontSize', 16);
tiledlayout('flow'); 
ylabels = ["Dt, nM", "Dp, nM", 'T, nM', 'DpT,nM'];
color = ['k', 'r', 'b', 'g', 'm'];
for jj = 1:length(dose)
for ii = 1:4
subplot(2,2,ii);
plot(sd(jj).Time/7, sd(jj).Data(:, ii), 'color', 'k', ...
    'LineWidth',2);
hold on;
xticks(0:1:10);
xlabel('Time, week'); ylabel(ylabels(ii));
if ii == 2
set(gca, 'Yscale', 'log');
end
end
end





