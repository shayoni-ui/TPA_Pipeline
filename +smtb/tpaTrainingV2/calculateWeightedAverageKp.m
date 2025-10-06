function kpAvg= calculateWeightedAverageKp(kp, tis)
% calcuateWeightedAverageKp is a function to calculate the weighted
% average kp to estimate fut
% Parameters
%   kp: strucuture containting estimated kp values for different tissue
%   tis: tissue contains informatio nof tissue volumes and density saved in
%   tis.mat
% Outputs:
%   kp_avg: Average kp value for the compound
arguments
    kp; % kp structure for the virtual compound
    tis; % tis contains information of tissue volumes and density
end

    tis = [[tis(1:13).volume]' [tis(1:13).density]'];
    wgts = tis(:,1).*tis(:,2)./(sum(tis(:,1).*tis(:,2)));

    kpAvg = sum(kp'.*wgts);
    %kp_fut_avg.fut_avg = sum([kp.Kp.fut]'.*wgts);
end