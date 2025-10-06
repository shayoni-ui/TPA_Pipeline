function [AUC, AUMC, MRT]  = trapzlog(t,c,method) 
arguments 
	 t; 
	 c; 
	 method {mustBeMember(method, {'linuplogdown', 'linlog'})} = ...
         'linuplogdown'; 
end 
% trapzlog provides rapezoidal numerical integration using the logarithmic
% trapezoidal rule after Cmax. TRAPZLOG(t,c) computes the integral of c
% with respect to t using the linuplogdown trapezoidal method.  
% t and c must be vectors of the same length, 
% or t must be a column vector and c an array whose first 
% non-singleton dimension is length(X).  
% TRAPZLOG operates along this dimension. 
% trapzlog calculates the area under the curve between
% two timepoints as AUC_interval = delta_t * ( c2 - c1 ) / ln (c2/c1). For
% intervals prior to tmax, or where c1 or c2 are zero, or the data is 
% increasing, the linear trapezoidal rule is used.
% Parameters: 
% 	 t: time
% 	 c: conc
% 	 method: 
% Outputs 
% 	 auc: 
% 	 aumc: 
% Date modified: 09-Oct-2022 
% File modified by SMTB user: jbp17697 
% File created by SMTB user 
% Examples: 
% 
%   See also TRAPZ, SUM, CUMSUM, CUMTRAPZ, INTEGRAL.



if sum(isnan(c))
    t = t(~isnan(c));
    c = c(~isnan(c));
    warning('Ignored NaN values for AUC calculations')
end

if length(t) ~= length(c), error(['The length of the x and y' ...
        ' vectors must be the same']); end

[t,idx] = sort(t);
c = c(idx);

z_int = zeros(length(t) - 1,1);
z_mc  = zeros(length(t) - 1,1);
for i = 2:length(t)
    if c(i) == 0 || c(i-1) == 0
        [z_int(i),z_mc(i)] = calcInterval(c(i-1:i),t(i-1:i),'lin');
    else
        switch lower(method)
            case 'linlog'
                [z_int(i),z_mc(i)] = calcInterval(c(i:i-1),t(i:i-1),'log');
            case 'linuplogdown'
                % Use linear if t is prior to Cmax, the data is 
                % increasing, or the datapoints are the same
                if c(i) >= c(i-1)
                    [z_int(i),z_mc(i)] = calcInterval(c(i-1:i), ...
                        t(i-1:i),'lin');
                else
                    [z_int(i),z_mc(i)] = calcInterval(c(i-1:i), ...
                        t(i-1:i),'log');
                end
            otherwise, error('Unrecognized method: %s',method)
        end
    end
end
AUC = sum(z_int);
AUMC = sum(z_mc);
MRT = AUMC ./ AUC;

%NCA = struct('AUC',AUC,'AUMC',AUMC,'MRT',MRT);

end

function [z,z_mc] = calcInterval(y,x,method)

dx = x(2) - x(1);
switch method
    case 'lin'
        z = dx * (y(2) + y(1))/2;
        z_mc = dx * (x(2) * y(2) + x(1) * y(1))/2;
    case 'log'
        z = dx * (y(2) - y(1))/log(y(2)/y(1));
        z_mc = dx * (x(2) * y(2) - x(1) * y(1))/log(y(2)/y(1)) -...
            dx^2 * (y(2) - y(1))/log(y(2)/y(1))^2;
    otherwise, error('Invalid interval method')
end

end




