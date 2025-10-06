function m1m2 = mrgMod(m1, m2) 
arguments 
	 m1; 
	 m2; 
end 
% mrgMod merges two simbiolgy models
% Parameters: 
% 	 m1: first simbiology model
% 	 m2: second simbiology model
% Outputs 
% 	 m1m2: merged simbiology model
% Date modified: Oct 8th 2020 
% File modified by SMTB user: jbp17697 
% File created by SMTB user: Jaimit Parikh
% Examples: 
% 

[m1ReactionsList, m1ReactionRatesList] = ...
    correctReactionsIfSinglecompartment(m1);
[m2ReactionsList, m2ReactionRatesList] = ...
    correctReactionsIfSinglecompartment(m2);

m1m2 = initializeMergedModel(m1, m2);

mergeCompartments(m1m2, m1, m2);
mergeParameters(m1m2, m1, m2);
mergeSpecies(m1m2, m1, m2);
mergeRules(m1m2, m1, m2);
mergeReactions(m1m2, m1, m2, ... 
    m1ReactionsList, m1ReactionRatesList, ...
    m2ReactionsList, m2ReactionRatesList);
copydoses(m1m2, m1, m2);

end


%% Helper function to correct for the reaction names given a single 
% compartment 
function [reactionsList, reactionRatesList] = ...
    correctReactionsIfSinglecompartment(m)
reactionsList = []; reactionRatesList = [];
if numel(m.Compartments) > 1
    return
end
reactionsList = string({m.Reactions.Reaction});
reactionRatesList = string({m.Reactions.ReactionRate});
rulesList = string({m.Rules.Rule});
speciesList = string({m.Species.Name});

for ii = speciesList
    reactionsList = regexprep(reactionsList, sprintf('(?<!\\w)%s(?!\\w)', ...
        ii), string(m.Compartments.Name)+ "." + ii);
    reactionRatesList = regexprep(reactionRatesList, ...
        sprintf('(?<!\\w)%s(?!\\w)', ii), ...
        string(m.Compartments.Name)+ "." + ii);
    rulesList = regexprep(rulesList, ...
        sprintf('(?<!\\w)%s(?!\\w)', ii), ...
        string(m.Compartments.Name)+ "." + ii);
end

% for ii = 1:length(m.Reactions)
%     m.Reactions(ii).Reaction = reactionsList(ii);
% end
% 
% m.Reactions

end

%%
function model = initializeMergedModel(m1, m2)

% can be improved to check if the tags, notes or name exists and construct
% names accordingly
newName = string(get(m1).Name) + " + "+ string(get(m2).Name);
newNotes = string(get(m1).Notes) + " + "+ string(get(m2).Notes);
newTag  = string(get(m1).Tag) + " + "+ string(get(m2).Tag);
model = sbiomodel(newName, 'Notes', newNotes, 'Tag', newTag);
end

%%
function mergeCompartments(m1m2, m1, m2)

% Add all compartments of m1 to m1m2
for ii = 1:numel(m1.Compartment)
    copyCompartments(m1m2, m1, ii)
end
% Add unique compartments of m2 not present in m1 to m1m2

[uniqueCompartments, uniqueCompartmentsIdx] = ... 
    setdiff(string(string({m2.Compartment.Name})), ...
        string(string({m1.Compartment.Name})));
    
for i = 1:numel(uniqueCompartments)
    ii = uniqueCompartmentsIdx(i);
    copyCompartments(m1m2, m2, ii)
end

end


function copyCompartments(m1m2, m, ii)

   addcompartment(m1m2, m.Compartment(ii).Name, ... 
       'Value', m.Compartment(ii).Value, ... 
       'Units', m.Compartment(ii).Units, ...
       'Constant', m.Compartment(ii).Constant);
end

%%
function mergeParameters(m1m2, m1, m2)

% Add all parameters of m1 to m1m2
for ii = 1:numel(m1.Parameters)
    copyParameters(m1m2, m1, ii)
end

% Add unique parameters of m2 not present in m1 to m1m2
% 
[uniqueParameters, uniqueParametersIdx] = ... 
    setdiff(string(string({m2.Parameters.Name})), ...
        string(string({m1.Parameters.Name})));
    
for i = 1:numel(uniqueParameters)
    ii = uniqueParametersIdx(i);
    copyParameters(m1m2, m2, ii)
end

end


function copyParameters(m1m2, m, ii)
    addparameter(m1m2, m.Parameters(ii).Name, ... 
       'Value', m.Parameters(ii).Value, ... 
       'Units', m.Parameters(ii).Units, ...
       'Constant', m.Parameters(ii).Constant, ...
       'ValueUnits', m.Parameters(ii).ValueUnits);

end

%%


function mergeSpecies(m1m2, m1, m2)

compartmentNames = string(string({m1m2.Compartments.Name}));

% Add all species of m1 to m1m2
for ii = 1:numel(m1.Species)
   idxCompartment = find(strcmp(m1.Species(ii).Parent.Name, ...
       compartmentNames));
   compObj = m1m2.Compartments(idxCompartment); 
   copySpecies(compObj, m1, ii)
end

% Add unique parameters of m2 not present in m1 to m1m2
% 
compNameM2 = string(arrayfun(@(x)(x.name), [m2.Species.Parent],...
    'UniformOutput', false));
compNameM1 = string(arrayfun(@(x)(x.name), [m1.Species.Parent],...
    'UniformOutput', false));

[uniqueSpecies, uniqueSpeciesIdx] = ... 
    setdiff(compNameM2 + string({m2.Species.Name}), ...
        compNameM1 + string({m1.Species.Name}));
    
for i = 1:numel(uniqueSpecies)
    ii = uniqueSpeciesIdx(i);
    idxCompartment = find(strcmp(m2.Species(ii).Parent.Name, ...
       compartmentNames));
    compObj = m1m2.Compartments(idxCompartment); 
  
    copySpecies(compObj, m2, ii)
end
    
function copySpecies(compObj, m, ii)
    addspecies(compObj, m.Species(ii).Name, ... 
       'Value', m.Species(ii).Value, ... 
       'Units', m.Species(ii).Units, ...
       'Constant', m.Species(ii).Constant, ...
       'InitialAmount', m.Species(ii).InitialAmount, ...
       'InitialAmountUnits', m.Species(ii).InitialAmountUnits, ...
       'ConstantAmount', m.Species(ii).ConstantAmount, ...
       'Notes', m.Species(ii).Notes, ...
       'Tag', m.Species(ii).Tag);
end


end


%%
function mergeReactions(m1m2, m1, m2, m1R, m1Rt, m2R, m2Rt)

% Add all compartments of m1 to m1m2
for ii = 1:numel(m1.Reactions)
    copyReactions(m1m2, m1, ii, m1R, m1Rt)
end
% Add unique compartments of m2 not present in m1 to m1m2

for ii = 1:numel(m2.Reactions)
    copyReactions(m1m2, m2, ii, m2R, m2Rt)
end
    
end


function copyReactions(m1m2, m, ii, mR, mRt)
    if numel(m.Compartments) > 1
        addreaction(m1m2, m.Reactions(ii).Reaction, ... 
            'ReactionRate', m.Reactions(ii).ReactionRate, ... 
            'Notes', m.Reactions(ii).Notes, ...
            'Name', m.Reactions(ii).Name);
    else
        addreaction(m1m2, mR(ii), ...
            'ReactionRate', mRt(ii), ... 
            'Notes', m.Reactions(ii).Notes, ...
            'Name', m.Reactions(ii).Name);
        
    end
end

%%

function leftSide = getLeftPartRule(rule)
temp = strsplit(rule, '=');
leftSide = temp(1).strip();
end

function mergeRules(m1m2, m1, m2)

% Add all compartments of m1 to m1m2
for ii = 1:numel(m1.Rules)
    copyRules(m1m2, m1, ii)
end
% Add unique compartments of m2 not present in m1 to m1m2

m2RulesLeft = arrayfun(@getLeftPartRule, ...
    string(string({m2.Rules.Rule})));
m1RulesLeft = arrayfun(@getLeftPartRule, ...
    string(string({m1.Rules.Rule})));

if isempty(m2RulesLeft)
        return; 
elseif isempty(m1RulesLeft)
    for ii = 1:numel(m2.Rules)
        copyRules(m1m2, m2, ii)
    end
else
    [uniqueRules, uniqueRulesIdx] = ... 
    setdiff(m2RulesLeft, m1RulesLeft);
    
    for i = 1:numel(uniqueRules)
        ii = uniqueRulesIdx(i);
        copyRules(m1m2, m2, ii)
    end
end
    
end


function copyRules(m1m2, m, ii)
    addrule(m1m2,  m.Rules(ii).Rule, ... 
       'RuleType', m.Rules(ii).RuleType, ... 
       'Notes', m.Rules(ii).Notes, ...
       'Name', m.Rules(ii).Name);
end

function copydoses(m1m2, m1, m2)
    dm1 = m1.getdose;
    dm2 = m2.getdose;

    for ii = 1:length(dm1)
        m1m2.adddose(dm1(ii));
    end

    doseNames = string({dm1.Name});

    for ii = 1:length(dm2)
        if ~any(contains(doseNames, dm2(ii).Name))
        m1m2.adddose(dm2(ii));
        else
            warning (fprintf(['Dose with same name already' ...
                ' added previously %s'], dm2(ii).Name))
        end
    end


end
