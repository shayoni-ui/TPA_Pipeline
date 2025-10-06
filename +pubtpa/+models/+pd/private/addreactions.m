
function addreactions(m1, reactionList, reactionProps)
arguments
    m1; % model
    reactionList; % Cellarray of compartments to be added
    reactionProps = {'Reaction', 'ReactionRate', 'Name'};
    % Cellarray of compartment properties to be added
end

for comp = reactionList
    temp = [reactionProps; comp{:}];
    m1.addreaction(temp{2:end});
end

end
