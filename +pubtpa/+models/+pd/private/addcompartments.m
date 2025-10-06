
function addcompartments(m1, compList, cmpProps)
arguments
    m1; % model
    compList; % Cellarray of compartments to be added
    cmpProps = {'Name', 'Value', 'Units', 'Constant'};
    % Cellarray of compartment properties to be added
end

for comp = compList
    temp = [cmpProps; comp{:}];
    m1.addcompartment(temp{2:end});
end

end
