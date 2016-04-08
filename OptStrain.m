
cvx_solver gurobi
model = struct(); % Declaration for explicit variable scope
load('data/yeast7.mat'); % Initialize base model
objectiveRxnId = 'r_9999';
keggToYeastPath = 'data/keggToYeast.json';
unidbPath = 'data/unidb.xlsx';

% Optstrain step 1 - add new reactions
model = augmentModel(model, objectiveRxnId, keggToYeastPath, unidbPath);

% Optstrain step 2 - find the base yield for the augmented model
[model, sol] = calcBaseYield(model, objectiveRxnId);

% Optstrain step 3 - remove reactions, bounded by not going below a defined
% ratio of the base yield.
S=model.S;
lb=model.lb;
ub=model.ub;
c=model.c;
[nmets,nrxns]=size(S);
minObjectiveRatio = 0.8;
minObjective = minObjectiveRatio*sol.f;

% Find the minimal set of reactions supporting the above objective
% Let y(i)=1 if reaction i is included, zero otherwise.
cvx_begin
variable v(nrxns)
variable y(nrxns) binary

minimize sum(y)
subject to
S*v == 0
c'*v >= minObjective
for i=1:nrxns
    % Force flux to zero whenever y(i)=0
    lb(i)*y(i) <= v(i) <= ub(i)*y(i)
end
cvx_end

% If the optimization was successful then the variable
% cvx_status contains 'solved'
if ~strcmpi(cvx_status,'solved')
    fprintf('Warning! Optimization problem was not solved sucessfully\n')
    disp(cvx_status)
else
    fprintf('Minimum number of reactions: %1.4f\n', cvx_optval)
    idx=find(abs(y) > 1e-5);
    disp(model.rxns(idx))
end