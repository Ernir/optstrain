% Script for analyzing the effects of humulene production on the normal
% brewing process.
% NOT RELATED to OPTSTRAIN ITSELF.

load('data/yeast7.mat')
clf; close all;

% Identify relevant reactions
growthIdx = find(strcmp('growth', model.rxnNames));
oxygenIdx = find(strcmp('oxygen exchange', model.rxnNames));
glucoseIdx = find(strcmp('glucose exchange', model.rxnNames));
maltoseIdx = find(strcmp('maltose exchange', model.rxnNames));
ethanolIdx = find(strcmp('ethanol exchange', model.rxnNames));
isoAcetIdx = find(strcmp('isoamyl acetate exchange', model.rxnNames));
glycerolIdx = find(strcmp('glycerol exchange', model.rxnNames));
ureaIdx = find(strcmp('urea exchange', model.rxnNames));

% Create conditions similar to those in fermenting wort
% Results in 300 available carbon atoms
model = changeRxnBounds(model, model.rxns(oxygenIdx), 0, 'l');
model = changeRxnBounds(model, model.rxns(glucoseIdx), -10, 'l');
model = changeRxnBounds(model, model.rxns(maltoseIdx), -20, 'l');

% Add the possibility of producing humulene
objectiveRxnId = 'r_9999';
model = addReaction( model, ...
    {'r_9998','farnesyl-diphosphate diphosphate-lyase'}, ...
    ['s_0190 -> s_0633 + ' objectiveRxnId] ...
    );
model = addReaction(model, ...
    {objectiveRxnId,'humulene exchange'}, ...
    [objectiveRxnId ' -> '] ...
    );
humuleneIdx = find(strcmp('humulene exchange', model.rxnNames));

% Define some objectives for the default strain
% Metabolites with more carbon atoms receive a higher priority to
% compensate
model.c = zeros(size(model.c));
model.c(growthIdx) = 1; % Some growth is OK
model.c(ethanolIdx) = 1/30*2; % C2H6O
model.c(isoAcetIdx) = 1/30*7; % C7H14O2
model.c(glycerolIdx) = 1/30*3; % C3H8O3
model.c(ureaIdx) = -1; % Urea is highly undesirable in the product.

% Final drawing

% Compare the fluxes before and after humulene inclusion
interestingRxns = model.rxns([growthIdx, ethanolIdx, isoAcetIdx, ...
    glycerolIdx, humuleneIdx]);
fluxMatrix = [];
[~, maxFlux] = fluxVariability(model, 80, 'max', interestingRxns);
fluxMatrix = [fluxMatrix maxFlux];
model.c(strcmp(objectiveRxnId, model.rxns)) = 1/30*15; % C15H24
[~, maxFlux] = fluxVariability(model, 80, 'max', interestingRxns);
fluxMatrix = [fluxMatrix maxFlux];
bar(fluxMatrix);
ax = gca;
ax.XTickLabel = {'vöxtur', 'etanól', 'isoamýl-acetat', ...
    'glýseról', 'humulene'};
ax.XTickLabelRotation = 45;
ylabel('Hámarksflæði')

% Perform robustness analysis 
figure;
robustnessAnalysis(model, model.rxns(humuleneIdx));
ylabel('Gildi markfalls')
xlabel('Humuleneframleiðsla')