
function addspecies(m1, spList, spProps)
arguments
    m1; % model
    spList; % Cellarray of compartments to be added
    spProps = {'CompName', 'Name', 'Value', 'Units'};
    % Cellarray of compartment properties to be added
end

for sp = spList
    cmp = sbioselect(m1, 'Type', 'Compartment', 'Name', sp{1});
    temp = [spProps(2:end); sp{1}(2:end)];
    cmp.addspecies(temp{2:end});
end

end
