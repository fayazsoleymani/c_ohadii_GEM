function model= modifyGPRsBlast(model, blastpResultsPath)
    
    % description: replacing the model genes that are not for the organism
    % with the ones from the organism
    % input
    %   model
    %   blastpResultsPath   file containeing the blast hits
    % output
    %   model               output model

    blastResults= readtable(blastpResultsPath, "ReadVariableNames",false);
    blastResultsDict= dictionary();
    for i= 1:length(blastResults.Var1)
        blastResultsDict(blastResults.Var1{i})= blastResults.Var2{i};
    end

    blastResultsKeys= keys(blastResultsDict);
    for i= 1: length(model.grRules)
        grRule= model.grRules{i};
        newGenes= [];
        changed= false;
        for j= 1: length(blastResultsKeys)
            gene= blastResultsKeys{j};
            if contains(grRule, gene)
                if ~strcmp(blastResultsDict(gene), "") && ~isempty(find(ismember(model.genes, blastResultsDict(gene))))
                    newGenes= [newGenes; blastResultsDict(gene)];
                end
                changed= true;
            end
        end
        if changed
            newGenes= unique(newGenes);
            if isempty(newGenes)
                model= changeGeneAssoc(model,model.rxns{i},'',true);
            else
                model= changeGeneAssoc(model,model.rxns{i}, char(strjoin(newGenes, ' or ')), true);
            end
        end
    end

end