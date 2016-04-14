function [model, sol] = calcBaseYield(model)
oxygen = model.rxns(strcmp('oxygen exchange', model.rxnNames));
model = changeRxnBounds(model, oxygen, 0, 'l');
sol = optimizeCbModel(model);
end