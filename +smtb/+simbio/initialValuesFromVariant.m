function x0 = initialValuesFromVariant(me, v) 
arguments 
	 me SimBiology.export.Model; 
	 v SimBiology.Variant; 
end 
% initialValuesFromVariant allows to get the initial values array based on
% the input simbiology variant and exported model 
% Parameters: 
% 	 me: exported simbiology model
% 	 v: simbiology variant
% Outputs 
% 	 icval: initial value array based on the exported model and variant
%       for performing simulation
% Date modified: 18-Oct-2022 
% File modified by SMTB user: jbp17697 
% File created by SMTB user 
% Examples: 
%   icValue = initialValuesFromVariant(me, v);
%
% See also SimBiology.Variant, SimBiology.export.Model

icName = string({me.ValueInfo.QualifiedName})';
x0 = [me.ValueInfo.InitialValue]';
vName = cellfun(@(x)x{2}, v(1).Content, 'UniformOutput', false);
vValue = cell2mat(cellfun(@(x)x{4}, v(1).Content, 'UniformOutput',false));
idx = cell2mat(arrayfun(@(x)find(ismember(icName, x)), vName, ...
    'UniformOutput', false));
x0(idx) = vValue;

end 
