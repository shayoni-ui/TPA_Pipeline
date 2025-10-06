
% First build DrugProps object
dp = DrugProps('GSK3925204A');

% Build the PBPK base model
m = BuildPBPKStructure;
me = export(m);

% Now parameterize the variants
[m,inputs,v] = parameterizePBPK(m,dp);
% The output variable inputs are the inputs with sources that were used to
% parameterize the model from DrugProps
% The output variable v are the 1x4 cell arrays of variant content - these
% are needed to run with an exported model

%% This section is an example running with a "normal" (not exported) 
% SimBiology model

% You need to specify the species and state with which you would like to
% run. This is accomplished via variants

v_rat = getvariant(m,'Rat Fasted');
v_human = getvariant(m,'Human Fasted');

% Now get the dose (in this case run PO). The dose is a simbiology dose
% object
d_po = getdose(m,'PO');

% Define the configuration set
cs = getconfigset(m);
cs.StopTime = 24;

% Run rat with 1 mpk dose (0.25 mg)
d_po.Amount = 0.25;
sd_rat = sbiosimulate(m,cs,v_rat,d_po);

cp = selectbyname(sd_rat,'plasma.drug');
figure
plot(cp.Time,cp.Data)
title('Rat PK')

% Run human with 100 mg dose
d_po.Amount = 100;
d_po.Interval = 24; %NOTE - the
% function calculate_time_series_simulation_results requires
% a dose interval to calculate steady-state concentration
sd_human = sbiosimulate(m,cs,v_human,d_po);

cp = selectbyname(sd_human,'plasma.drug');
figure
plot(cp.Time,cp.Data)
title('Human PK')

% Now calculate time series endpoints (e.g., PK endpoints)
%res = calculate_time_series_simulation_results(sd_human,m,d_po);

%% This section is an example running with an exported SimBiology model

% To run with an exported model, you need to select the proper variant from
% variable "v" - A SIMBIOLOGY VARIANT WILL NOT WORK WITH THIS

% run with human as an example - get human variant content from variable v
v_human_ex = v(strcmpi(v(:,1),'Human Fasted'),:);

% applyVariant takes the 1x4 cell arrays and overwrites the original model
% values - note this is only passing the 1x4 cell arrays
me = applyVariant(me,v_human_ex{2});
d_po_ex = me.getdose('PO'); % note this is not a simbiology dose object
d_po_ex.Amount = 100;
d_po_ex.Interval = 24; % recall the 
% calculate_time_series_simulation_results requires a dose interval
me.SimulationOptions.StopTime = 24; % set simulation stop time at 24 hours

sd_human_ex = simulate(me,d_po_ex); % use simulate, not sbiosimulate
cp = selectbyname(sd_human_ex,'plasma.drug');
figure
plot(cp.Time,cp.Data)
title('Human PK - exported model')
save('cpDp.mat', 'cp');

%res_ex = calculate_time_series_simulation_results(sd_human_ex,me,d_po_ex);



