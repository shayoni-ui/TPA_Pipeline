function variant = updateVariantContent(variant, name, value)
arguments
    variant cell;
end
arguments(Repeating)
    name; % name of the parameter to be updated
    value;% value of the parameter to be updated
end
% updateVariantContent allows to easily update the variant content
% cell array
%
% Parameters:
%   variant cell; cell array of variant
%   repeating:
%   name; name of parameter to be updated
%   value; value of the paramter to be updated
% Outputs:
%   variant cell; updated variant cell array
%
% File modified by SMTB user:jbp17697
% File modified date: Oct-11-2022
% File created by SMTB user: Jaimit Parikh
% 
% Examples:
%   newVariant = updateVariantContent(variant, 'fup', 0.1);
%   newVariant = updateVariantContent(variant, 'fup, 0.1, 'B2P', 1.5);

for ii = 1:length(name)
    idx = strcmpi(cellfun(@(x)x{2}, ...
        variant, 'UniformOutput',false), name{ii});
    variant{idx}{4} = value{ii};
end
    
end