
function addparameters(m1, paramList, paramProps)
arguments
    m1; % model
    paramList; % Cellarray of compartments to be added
    paramProps = {'Name', 'Value', 'Units', 'Constant'};
    % Cellarray of compartment properties to be added
end

for comp = paramList
    temp = [paramProps; comp{:}];
    m1.addparameter(temp{2:end});
end

end

