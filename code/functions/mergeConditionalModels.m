function newModels= mergeConditionalModels(models)
    % description: merging all the models corresponding to different growth
    % conditions to one model
    % input     all input models
    % output    output model

    nModels= length(models);

    unionRxnsDict= dictionary(string([]), double([]));
    for i= 1:nModels
        tempRxns= models{i}.rxns;
        for r= 1:length(tempRxns)
            if ~isKey(unionRxnsDict, tempRxns{r})
                unionRxnsDict(tempRxns{r})= i;
            end
        end
    end
    fprintf("mergedModel nRxns: %d\n", length(keys(unionRxnsDict)));
    
    
    unionMetsDict= dictionary(string([]), double([]));
    for i= 1:nModels
        tempMets= models{i}.mets;
        for m= 1:length(tempMets)
            if ~isKey(unionMetsDict, tempMets{m})
                unionMetsDict(tempMets{m})= i;
            end
        end
    end
    fprintf("mergedModel nMets: %d\n", length(keys(unionMetsDict)));

    unionGenesDict= dictionary(string([]), double([]));
    for i= 1:nModels
        tempGenes= models{i}.genes;
        for g= 1:length(tempGenes)
            if ~isKey(unionGenesDict, tempGenes{g})
                unionGenesDict(tempGenes{g})= i;
            end
        end
    end
    fprintf("mergedModel nGenes: %d\n", length(keys(unionGenesDict)));

    unionRxns= keys(unionRxnsDict);

    newModels= {};
    for i=1:nModels
        newModel= models{i};
        for r= 1:length(unionRxns)
            tempRxn= unionRxns{r};
            if ~ismember(tempRxn, newModel.rxns)
                originalModelIndex= unionRxnsDict(tempRxn);
                originalModel= models{originalModelIndex};
                originalRxnIndex= find(strcmp(tempRxn, originalModel.rxns));
                [newModel, ~] = addReaction(newModel, tempRxn,...
                    'reactionName', originalModel.rxnNames{originalRxnIndex},...
                    'metaboliteList', originalModel.mets(find(originalModel.S(:, originalRxnIndex))),...
                    'stoichCoeffList', full(originalModel.S(find(originalModel.S(:, originalRxnIndex)), originalRxnIndex)),...
                    'reversible', originalModel.rev(originalRxnIndex),...
                    'lowerBound', 0,...
                    'upperBound', 0,...
                    'geneRule', originalModel.grRules{originalRxnIndex});
                if isfield(originalModel, 'rxnECNumbers')
                    newModel.rxnECNumbers(end)= originalModel.rxnECNumbers(originalRxnIndex);
                end
                if isfield(originalModel, 'rxnReferences')
                    newModel.rxnECNumbers(end)= originalModel.rxnReferences(originalRxnIndex);
                end
                if isfield(originalModel, 'rxnKEGGID')
                    newModel.rxnKEGGID(end)= originalModel.rxnKEGGID(originalRxnIndex);
                end
                if isfield(originalModel, 'rxnRheaID')
                    newModel.rxnRheaID(end)= originalModel.rxnRheaID(originalRxnIndex);
                end
                if isfield(originalModel, 'subSystems')
                    newModel.subSystems(end)= originalModel.subSystems(originalRxnIndex);
                end
            end
        end
        newModel.id= strrep(models{i}.id, ...
            num2str(length(models{i}.genes)), ...
            num2str(length(newModel.genes))); 
        newModel= clearComps(newModel, false, false, false);
        newModels= [newModels; newModel];
    end

    

    
    
    
    


end