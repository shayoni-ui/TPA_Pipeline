%% Perform mode2 analysis. See below runMode2Analysis function for details


%%
N = 100;
% 1. Indirect response model fast turnover rate
fup =0.01; B2P = 1;
pdConstPars.Imax = 1; pdConstPars.RC50 = 1; pdConstPars.ksyn_P = 1;
pdConstPars.kdeg_P = 1;
runMode2Analysis('indirectresponse', fup, B2P, "P", "indirectFast01.mat", ...
    pdConstPars, N);

%% 
fup =0.05; B2P = 1;
pdConstPars.Imax = 1; pdConstPars.RC50 = 1; pdConstPars.ksyn_P = 1;
pdConstPars.kdeg_P = 1;
runMode2Analysis('indirectresponse', fup, B2P, "P", "indirectFast05.mat", ...
    pdConstPars, N);

fup =0.1; B2P = 1;
pdConstPars.Imax = 1; pdConstPars.RC50 = 1; pdConstPars.ksyn_P = 1;
pdConstPars.kdeg_P = 1;
runMode2Analysis('indirectresponse', fup, B2P, "P", "indirectFast1.mat", ...
    pdConstPars, N);


fup =0.2; B2P = 1;
pdConstPars.Imax = 1; pdConstPars.RC50 = 1; pdConstPars.ksyn_P = 1;
pdConstPars.kdeg_P = 1;
runMode2Analysis('indirectresponse', fup, B2P, "P", "indirectFast2.mat", ...
    pdConstPars, N);

%%
% 2. Indirect response model slow turnover rate
fup =0.01; B2P = 1;
pdConstPars.Imax = 1; pdConstPars.RC50 = 1; pdConstPars.ksyn_P = 0.01;
pdConstPars.kdeg_P = 0.01;
runMode2Analysis('indirectresponse', fup, B2P, "P", "indirectSlow01.mat", ...
    pdConstPars, N); 

fup =0.05; B2P = 1;
pdConstPars.Imax = 1; pdConstPars.RC50 = 1; pdConstPars.ksyn_P = 0.01;
pdConstPars.kdeg_P = 0.01;
runMode2Analysis('indirectresponse', fup, B2P, "P", "indirectSlow05.mat", ...
    pdConstPars, N); 


fup =0.1; B2P = 1;
pdConstPars.Imax = 1; pdConstPars.RC50 = 1; pdConstPars.ksyn_P = 0.01;
pdConstPars.kdeg_P = 0.01;
runMode2Analysis('indirectresponse', fup, B2P, "P", "indirectSlow1.mat", ...
    pdConstPars, N); 


fup =0.2; B2P = 1;
pdConstPars.Imax = 1; pdConstPars.RC50 = 1; pdConstPars.ksyn_P = 0.01;
pdConstPars.kdeg_P = 0.01;
runMode2Analysis('indirectresponse', fup, B2P, "P", "indirectSlow2.mat", ...
    pdConstPars, N); 
%%
% 3. Precursor activation model
clear pdConstPars;
fup = 0.1; B2P = 1;
pdConstPars.Imax = 1; pdConstPars.RC50 = 1;
runMode2Analysis('precursoractivation', fup, B2P, "A", "precursor.mat", ...
    pdConstPars, N);
%%
% 4. Transduction model
fup =0.1; B2P = 1;
pdConstPars.Imax = 1; pdConstPars.RC50 = 1;
runMode2Analysis('transduction', fup, B2P, "M1", "transduction.mat", ...
    pdConstPars, N);
%%
% 5. Tolerance model
fup =0.1; B2P = 1;
pdConstPars.Imax = 1; pdConstPars.RC50 = 1;
runMode2Analysis('tolerance', fup, B2P, "P1", "tolerance.mat", ...
    pdConstPars, N);
%%
% 6. Tumor growth model
fup =0.1; B2P = 1;
pdConstPars.Imax = 0.003; pdConstPars.RC50 = 1;
runMode2Analysis('tumor', fup, B2P, "TGI", "tumor.mat", ...
    pdConstPars, N);



%% Helper functions
function results = runMode2Analysis(pd, fup, B2P, ...
    responsename, fname, pdConstPars, N)
% Get the integrated PBPK-PD exported and non-exported model
[m, mexp] = pubtpa.models.pbpkpd.createMergedPBPKPD(pd, ...
    true);


% Sample the PK and PD parameters that we want to vary and generate
% variant (virtual compounds) to perform the simulation
disp('Starting parameterization to get dervied parameters...')
tic
% Varying only CLint and logP.. additional PK parameters to be varied
% can be added here
sampleSetPK = smtb.gsa.getSamples('methodtype', 'lhs', 'N', N, ...
    'paramDist', ...
        {'CLint', 'LogUniform', {'Lower', 0.1, 'Upper', 300}, {},...
         'logP', 'uniform', {'Lower', 0.1, 'Upper', 5},         {},...
        }...
    );
pkParamNames = sampleSetPK.ParamNames;
pkParamValues = sampleSetPK.Samples{1};

% Other pk parameters set to constant value
pkConstPars.sol_mgmL = 10;
pkConstPars.Peff0 = 10;
pkConstPars.fup = fup;
pkConstPars.B2P = B2P;
pkConstPars.mw = 500;
pkConstPars.Deff = 5e-7;
pkConstPars.SF = 20;
pkConstPars.species = {'human'};
pkConstPars.fub = 'fup/B2P';
pkConstPars.fup_type = 'experimental'; % adjusted vs experimental
pkConstPars.ASF_method = 'Theoretical SAV';
pkConstPars.para_abs = 'off';
pkConstPars.RenalFiltration = 'off';
pkConstPars.FaSSGF = 0;
pkConstPars.FaSSIF = 0;
pkConstPars.FeSSIF = 0;
% Default adjusted however significant impact on the range of estimated
% Weighted avg Kp... experimental provides a wider range


%pdParamNames = {'Imax', 'RC50'};
%pdParamValues = ones(N, 1).*[Imax, XC50];
pdParamNames = {}; pdParamValues = [];
% pdConstPars.Imax = 1;
% pdConstPars.RC50 = XC50;

% generate VCs allows to generate virtual compounds given a name of the PK
% parameters that are varied (pkParamNames) and its values (pkParamValues)
% a structure with the name value of PK constant parameters (pkConstPars)
% and similar for the PD, pdParamNames, pdParamValues and pdConstPars
% The generateVCs function is part of the SMTB package
% written by Jaimit and not our SMTBToolbox used for the tutorial
[variantArr, qoiArr] = smtb.tpa.generateVCs(pkParamNames,...
    pkParamValues, pkConstPars, {}, [], pdConstPars);

toc
disp('Finished parameterization')
% Perform the simulation using the model and the generated virtual
% compound variants
disp('Starting simulation of virtual compounds')
tic


d = mexp.getdose('PO'); % Oral dosing
d.Amount = 100; d.Interval = 24; d.RepeatCount = 20;
mexp.SimulationOptions.StopTime = d.Interval * d.RepeatCount;


% The performVirtualExploration function is part of the SMTB package
% written by Jaimit and not our SMTBToolbox used for the tutorial
sd=smtb.tpa.performVirtualExploration(mexp, d,'variantArr',variantArr, ...
    'selectNames',["plasma.drug", "TargetTissue.drugU", ...
    responsename], "shrinkData",1);
toc
disp('Finished simulation of virtual compounds');

% Lot of flexiblity around saving 
disp('Saving results')
tic
results.model = m; results.mexp = mexp; results.qoiArr = qoiArr; 
results.sd = sd; results.N = N; results.variantArr = variantArr;
results.dose = d; results.file = which('mode2Analysis');
results.sampleSetPK = sampleSetPK;
results.pdParam.names = pdParamNames;
results.pdParam.values = pdParamValues;
results.pdParam.const = pdConstPars;
%clearvars -except results
save(strcat('./matFiles/', fname), 'results');
toc
disp('Finished saving results')
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
%         plot(sd{jj}.Time, sd{jj}.Data(:, ii)*scale, 'LineWidth',2);
%         hold on;
%     end
%     xlabel(sprintf('Time, %s', sd{jj}.TimeUnits));
%     ylabel(sprintf('%s, %s',sd{jj}.DataInfo{ii}.Name, units));
%     %xlim([0, 480]);
% 
% end
% 
% if fname
%     exporgraphics(f, fname, 'resolution', 300);
% end
% end
% 
% 

    

