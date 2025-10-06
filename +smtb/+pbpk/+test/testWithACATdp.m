
% First build DrugProps object
dp = smtb.DrugProps('GSK3925204A');

%% 
% Build the PBPK base model
m = smtb.pbpk.buildPBPKmodelTemplate;
me = export(m);

%%
% Now parameterize the variants and add to non exported simbiology model
[inputs,v] = smtb.pbpk.getPBPKVariant(dp, 'species', ...
    {'human'});

%%
x0 = smtb.simbio.initialValuesFromVariant(me, v(1));

%%
d = me.getdose('PO'); % note this is not a simbiology dose object
d.Amount = 100;
d.Interval = 24; % recall the 

me.SimulationOptions.StopTime = 24; % set simulation stop time at 24 hours
sd = simulate(me,x0, d); % use simulate, not sbiosimulate

%%
cp = selectbyname(sd,'plasma.drug');

f = figure ('DefaultAxesFontSize', 14);
plot(cp.Time,cp.Data,'LineWidth',2, 'DisplayName', 'Updated code');
hold on;
if exist("cpDp.mat",'file')
    load("cpDp.mat");
    plot(cp.Time, cp.Data, 'LineWidth',2, 'LineStyle','--', ...
    'color', 'c', 'DisplayName','Original');
end

xlabel('Time, hours');
ylabel(sprintf('%s, %s', cp.DataInfo{1}.Name,cp.DataInfo{1}.Units));
legend();



