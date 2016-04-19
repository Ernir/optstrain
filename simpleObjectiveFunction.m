function model = simpleObjectiveFunction(model)
% Before OptStrain is considered, make customized additions to the model.
% In this case, it means enabling humulene production and solely optimize
% for that.

% Identify relevant reactions
oxygen = strcmp('oxygen exchange', model.rxnNames);
glucose = strcmp('glucose exchange', model.rxnNames);
maltose = strcmp('maltose exchange', model.rxnNames);

% To simplify comparisions, use sugars similar to those available in wort
model = changeRxnBounds(model, model.rxns(oxygen), -1000, 'l');
model = changeRxnBounds(model, model.rxns(glucose), -10, 'l');
model = changeRxnBounds(model, model.rxns(maltose), -20, 'l');

% Add the possibility of producing humulene
humuleneRxn = 'r_9999';
model = addReaction( model, ...
    {'r_9998','farnesyl-diphosphate diphosphate-lyase'}, ...
    ['s_0190 -> s_0633 + ' humuleneRxn] ...
    );
model = addReaction(model, ...
    {humuleneRxn,'humulene exchange'}, ...
    [humuleneRxn ' -> '] ...
    );

model = changeObjective(model, humuleneRxn);

end