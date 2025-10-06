function  vTb = variantContentTable(v) 
arguments 
	 v {checkValidationOfV(v)}; %
end 
% variantContentTable converts either a simbiology variant obj or cell 
% array of variant information into a matlab table
% Parameters: 
% 	 v: cell array of variant information or simbiology variant object
% Outputs:
%    vTb: variant table
% Date modified: 04-Oct-2022 
% File modified by SMTB user: jbp17697 
% File created by SMTB user 
% Examples: 
% 

if isa(v, "SimBiology.Variant")
    v = v.Content;
end
% convert to table the variant content
vTb = cell2table(reshape([v{:}], 4, [])');
vTb.Properties.VariableNames = ["Type", "Name", "Property", "Value"];



end 


function checkValidationOfV(v)
   if isa(v, "SimBiology.Variant")
   elseif iscell(v)
       if length([v{:}]) ~= length(v) * 4
            error(['v must be a Simbiology Variant object or a '...
           'cell array with each element as 1x4 cell']);
       end
       
   else
       error(['v must be a Simbiology Variant object or a '...
           'cell array with each element as 1x4 cell'])
   end
end