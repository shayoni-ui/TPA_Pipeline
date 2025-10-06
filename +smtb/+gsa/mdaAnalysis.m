function S = mdaAnalysis(X, Y, regType, rensembleArgs)
arguments
    X; %  n x k ... n input samples with k parameters
    Y; % n x 1 ... nsamples with 1 output variable of interest
    regType = "rensemble"; % available options are
    rensembleArgs =[]; % structure containting argument-value pairs for fitrensemble
end
% mdaAnalysis perform mean decrease accuracy analysis to estimate the sensitivity of
% the parameters X to the response variable Y using nonlinear ensemble
% regression models
% Parameters:
% 	 X:  n x k ... n input samples with k parameters
% 	 Y: n x 1 ... nsamples with 1 output variable of interest
% 	 regType: Default "rensemble", available options are ..
% 	 rensembleArgs: strucuture to pass any name-value args supported by
% 	 fitrtree function
% Outputs
% 	 S: Structure with sensitivity indices
% Date modified: 27-Sep-2022
% File modified by SMTB user: jbp17697
% File created by SMTB user
% Examples:
%   s = mdaAnalysis(X, Y);

regType = string(regType);

X = normalize(X);
Y = normalize(Y);


switch regType % TODO: Add different regression
    % types in addition to random forest
    case "rensemble"
        if isempty(rensembleArgs)
            mdl = fitrtree(X, Y);
        else
            % convert structure wity P fields to 1 x 2P cell array for
            % arguments
            rensembleArgsCell = namedargs2cell(resensembleArgs);

            mdl = fitrensemble(X, Y, rensembleArgsCell{:});
        end

        Ynew = mdl.predict(X);
        S.acc = corr(Y, Ynew)^2;
        fprintf("The original R2 score is: %d \n", S.acc);
        S.ST = estimateMDA(mdl, X, Ynew, S);
        S.ST = S.ST / sum(S.ST);
end

end


function mda = estimateMDA(mdl, X, Y, S)
mda = zeros(1, size(X, 2));
for ii = 1:size(X, 2)
    Xshuffled = permuteCols(X, ii);
    Ypred = mdl.predict(Xshuffled);
    mda(ii) = S.acc - corr(Y, Ypred)^2;
end
end


function X = permuteCols(X, idx)
X(:, idx) = X(randperm(length(X)), idx);
end

