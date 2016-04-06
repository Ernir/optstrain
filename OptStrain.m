function OptStrain
    model = step1();

    function model = step1()
        % Step 1 of the OptStrain algorithm involves adding a curated
        % database of reactions to the model organism.
        
        load('data/yeast7.mat'); % Creates base model variable named model
        yeastToKegg = readJson('data/yeastToKegg.json');
        keggToYeast = readJson('data/keggToYeast.json');
        [~,~,RAW] = xlsread('data/unidb.xlsx');
        uniDB = RAW(2:end,:);
        
        % TODO fill in missing KEGG entries
        % TODO add the uniDB entries
    end
end