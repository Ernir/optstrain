% Script for analyzing the effects of humulene production on the normal
% brewing process.
% NOT RELATED to OPTSTRAIN ITSELF.

load('data/yeast7.mat')
clf; close all;

% Add humulene to our computational model and define objective function
model = addTargetMetabolite(model);
humulene = strcmp('humulene exchange', model.rxnNames);

% Compare the fluxes with and without humulene
interestingRxns = model.c > 0;
[~, maxFlux] = fluxVariability(model, 80, 'max', model.rxns(interestingRxns));
fluxMatrix = maxFlux;
model.c(humulene) = 0;
[~, maxFlux] = fluxVariability(model, 80, 'max', model.rxns(interestingRxns));
fluxMatrix = [maxFlux fluxMatrix];
bar(fluxMatrix);
ax = gca;
ax.XTickLabel = {'etanól', 'glýseról', 'isoamýl-acetat', 'vöxtur', 'humulene'};
ax.XTickLabelRotation = 45;
ylabel('Hámarksflæði')

% Perform robustness analysis 
figure;
robustnessAnalysis(model, model.rxns(humulene));
ylabel('Gildi markfalls')
xlabel('Humuleneframleiðsla')