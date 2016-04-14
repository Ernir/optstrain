function model = augmentModel(model, keggToYeastPath, unidbPath)
% Step 1 of the OptStrain algorithm involves adding a curated
% database of reactions to the model organism.
keggToYeast = readJson(keggToYeastPath);
[~,~,RAW] = xlsread(unidbPath);
uniDB = RAW(2:end,:);
% Add reactions from the "universal" database to the model
reactionIdColumn = 2;
reactionColumn = 3;
reactionNameColumn = 4;
newReactionNum = 3000; % IDs not present in DB need to be created
newCompoundNum = 20000;
skip = false;
newRxnIds = {};
for i = 1:length(uniDB)
    
    % Find a the new reaction's ID
    try
        newRxnId = keggToYaeast.(uniDB(i, reactionIdColumn));
    catch
        % Create one if it isn't
        newRxnId = ['r_' num2str(newReactionNum)];
        newReactionNum = newReactionNum + 1;
    end
    newRxnIds = [newRxnIds {newRxnId}];
    
    % If the reaction isn't already in the model, add it
    if isempty(find(strcmp(model.rxns,newRxnId), 1))
        
        % Create a new reaction from KEGG IDs
        keggReaction = uniDB(i, reactionColumn);
        parts = strsplit(keggReaction{1});
        n = length(parts);
        newParts = cell(1,n);
        for j = 1:n
            part = parts(j);
            part = part{1}; % Content referencing
            if strcmp(part(1),'C') % Metabolite IDs start with C
                try
                    modelId = keggToYeast.(part);
                catch
                    modelId = ['s_' num2str(newCompoundNum)];
                    newCompoundNum = newCompoundNum + 1;
                end
                newParts{j} = modelId;
            elseif strcmp(part(1),'G') % Some are weird, skip those
                newParts{j} = part;
                skip = true;
            else
                newParts{j} = part;
            end
        end
        newRxn = strjoin(newParts);
        % Add the newly created reaction
        newRxnName = uniDB(i, reactionNameColumn);
        if ~skip
            warning('off','all') % This step is noisy
            addReaction(model, {newRxnId, newRxnName}, newRxn);
            warning('on','all')
        else
            skip = false;
        end
    end
    
end
% Mark new reactions as non-native
model.native = ~ismember(model.rxns, newRxnIds);
end