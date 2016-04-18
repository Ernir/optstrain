cvx_solver gurobi
verbose = true;
load('data/yeast7.mat'); % Initialize base model
keggToYeastPath = 'data/keggToYeast.json';
unidbPath = 'data/unidb.xlsx';

% Make our custom (problem-specfic) additions to the model
model = addTargetMetabolite(model);

% Optstrain step 1 - add a universal database of new reactions
verbosePrint('Starting OptStrain step 1', verbose);
[model] = augmentModel(model, keggToYeastPath, unidbPath);

% Optstrain step 2 - find the base objective function maximum value for 
% the augmented model
verbosePrint('Starting OptStrain step 2', verbose);
[model, sol] = calcBaseYield(model);
disp(['Objective function value is ' num2str(sol.f)])

% Optstrain step 3 - remove reactions, bounded by not going below a defined
% ratio of the pre-calculated maximum value
% This part is based heavily on examples by Steinn Guðmundsson
verbosePrint('Starting OptStrain step 3', verbose);
S=model.S;
lb=model.lb;
lb(isinf(lb))=-1000; % Removing infs, which are problematic for gurobi
ub=model.ub;
ub(isinf(ub))=1000;
c=model.c;
native=model.native;
[nmets,nrxns]=size(S);

minObjectiveRatio = 0.8;
minObjective = minObjectiveRatio*sol.f;

% Find the minimal set of reactions supporting the above objective
% Let y(i)=1 if reaction i is included, zero otherwise.
cvx_begin
variable v(nrxns)
variable y(nrxns) binary

minimize sum(y(~native)) % Minimize, excluding native reactions
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

% Optstrain step 4 - attempt to find growth coupling reaction deletions

% Exclude those reactions we didn't want:
excludedRxns = setdiff(idx, 1:length(model));
for rxn = excludedRxns
    model = changeRxnBounds(model, rxn, 0, 'b');
end
% Use GDLS to find reaction deletion
% Only consider humulene production here, general coupling is infeasible
verbosePrint('Starting OptStrain step 4', verbose);
model.c = double(strcmp('humulene exchange', model.rxnNames));
[gdlsSolution, bilevelMILPProblem, gdlsSolutionStructs] = ...
    GDLS(model, model.rxns(model.c > 0));