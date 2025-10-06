function val = getVariantContent(v,name) 
arguments 
	 v SimBiology.Variant;  
	 name char; 
end 
% getVariantContent allows to get the value of the parameter in 
% simbiology variant object
% Parameters: 
% 	 v: SimBiology.Variant object
% 	 name: name of the parameter you want the value of 
% Outputs 
% 	 val: Value of the parameter
% Date modified: 19-Oct-2022 
% File modified by SMTB user: jbp17697 
% File created by SMTB user 
% Examples: 
%   val = getVariantContent(v, 'fup');
% 
% See also: SimBiology.Variant

idx = strcmpi(cellfun(@(x)x{2}, v.Content, ...
    'UniformOutput', false), name);

values = cellfun(@(x)x{4}, v.Content, 'UniformOutput',false);
val = values{idx};

end 
