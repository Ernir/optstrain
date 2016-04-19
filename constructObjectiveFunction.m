function model = constructObjectiveFunction(model)
% Before OptStrain is considered, make customized additions to the model.
% In this case, it means enabling humulene production

% Identify relevant reactions
growth = strcmp('growth', model.rxnNames);
oxygen = strcmp('oxygen exchange', model.rxnNames);
glucose = strcmp('glucose exchange', model.rxnNames);
maltose = strcmp('maltose exchange', model.rxnNames);
ethanol = strcmp('ethanol exchange', model.rxnNames);
isoAcet = strcmp('isoamyl acetate exchange', model.rxnNames);
glycerol = strcmp('glycerol exchange', model.rxnNames);
urea = strcmp('urea exchange', model.rxnNames);

% Create conditions similar to those in fermenting wort
% Results in 300 available carbon atoms
model = changeRxnBounds(model, model.rxns(oxygen), 0, 'l');
model = changeRxnBounds(model, model.rxns(glucose), -10, 'l');
model = changeRxnBounds(model, model.rxns(maltose), -20, 'l');

% Add the possibility of producing humulene
humuleneRxnId = 'r_9999';
model = addReaction( model, ...
    {'r_9998','farnesyl-diphosphate diphosphate-lyase'}, ...
    ['s_0190 -> s_0633 + ' humuleneRxnId] ...
    );
model = addReaction(model, ...
    {humuleneRxnId,'humulene exchange'}, ...
    [humuleneRxnId ' -> '] ...
    );
humulene = strcmp('humulene exchange', model.rxnNames);

% Define some objectives for the default strain
% Metabolites with more carbon atoms receive a higher priority to
% compensate
model.c = zeros(size(model.c));
model.c(growth) = 0.1; % Some growth is OK
model.c(ethanol) = 1/30*2; % C2H6O
model.c(isoAcet) = 1/30*7; % C7H14O2
model.c(glycerol) = 1/30*3; % C3H8O3
model.c(urea) = -1; % Urea is highly undesirable in the product.
model.c(humulene) = 1/30*15; % C15H24

end