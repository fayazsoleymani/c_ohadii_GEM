function [model, optimizableMets, gapfillableMets, nonOptimizebleMets, results, errors]= ...
    GapFillBiomassMets(model, models, minFlux, gapfillResultsOutputDir)

    % desctiption: trying to gapfill the model with reference models
    % inputs
    %   model                   input model
    %   models                  reference models
    %   minFlux                 minimum flux for the demand constraint
    %   gapfillREsultsOutputDir directory to save the results
    %
    % outputs
    %   model                   gapfilled model
    %   optimizableMets         biomass precursors that can be produced
    %   gapfillableMets         biomass precursors that are gapfilled
    %   nonOptimizebleMets      biomass precursors that can't be produced
    %   results                 complete results
    %   errors                  errors during gapfilling

    if ~exist(gapfillResultsOutputDir, "dir")
        mkdir(gapfillResultsOutputDir);
    end
    
    nOrgGenes= length(model.genes);
    nOrgRxns= length(model.rxns);
    
    if strcmp(models, 'kegg_metacyc')
        keggModel=getModelFromKEGG([],false,false,false,false);
        keggModel.grRules= cell(keggModel.rxns);
        metacycModel= getModelFromMetaCyc([], false, false, false);
        refModel=combineMetaCycKEGGModels(metacycModel,keggModel);
        refModel=rmfield(refModel,'genes');
        refModel=rmfield(refModel,'rxnGeneMat');
        % this was commented
        refModel=rmfield(refModel, 'pwys');
        balanceStructure=getElementalBalance(refModel);
        refModel=removeReactions(refModel,balanceStructure.balanceStatus~=1,true,true);
    else
        refModel= models;
        if ~isfield(model, 'combinedMetBiggIDs')
            model.combinedMetBiggIDs= cell(length(model.mets), 1);
        end
        for metIndex= 1: length(model.mets)
            if metIndex <= length(model.metBiGGID)
                if ~isempty(model.metBiGGID{metIndex})
                    model.combinedMetBiggIDs(metIndex)= model.metBiGGID(metIndex);
                else
                    model.combinedMetBiggIDs(metIndex)= model.mets(metIndex);
                end
            else
                model.combinedMetBiggIDs(metIndex)= model.mets(metIndex);
            end
        end
        originalMets= model.mets;
        model.mets= model.combinedMetBiggIDs;
    end

    biomassIndices= findBiomassIndices(model);
    biomassIndex= biomassIndices(1);
    biomassMetIndices= find(model.S(:, biomassIndex));
    optimizableMets= {};
    errors= {};
    gapfillableMets= {};
    nonOptimizebleMets= {};
    results= dictionary;

    for metIndex= 1:length(biomassMetIndices)
        biomassMetIndex= biomassMetIndices(metIndex);
        biomassMetName= model.mets(biomassMetIndex);
        model = addDemandReaction(model, biomassMetName);
        model.lb(end)= minFlux;
        model.c(end)= 1;
        solution = optimizeCbModel(model);
        if solution.stat == 1
            optimizableMets= [optimizableMets; biomassMetName{1}];
            results(biomassMetName{1})= solution;
        else
            % gapfile with error handling
            try
                [~, ~, addedRxns, newModel, exitFlag]= ...
                fillGaps(model, refModel, false , true, false);
                solution = optimizeCbModel(newModel);
                tempStruct= struct();
                if ~isfield(tempStruct, 'metIndex')
                    tempStruct.('metIndex')= biomassMetIndex;
                end
                if ~isfield(tempStruct, 'metName')
                    tempStruct.('metName')= biomassMetName{1};
                end
                if ~isfield(tempStruct, 'addedRxns')
                    tempStruct.('addedRxns')= addedRxns;
                end
                if ~isfield(tempStruct, 'newModel')
                    tempStruct.('newModel')= newModel;
                end
                if ~isfield(tempStruct, 'exitFlag')
                    tempStruct.("exitFlag")= exitFlag;
                end
                if ~isfield(tempStruct, 'solution')
                    tempStruct.('solution')= solution;
                end
                if solution.stat == 1
                    gapfillableMets= [gapfillableMets; biomassMetName{1}];
                    results(biomassMetName{1})= tempStruct;
                    model= newModel;
                    resultsOutputFile= fullfile(gapfillResultsOutputDir, strcat(biomassMetName{1}, '.txt'));
                    fileID= fopen(resultsOutputFile, 'w');
                    fprintf(fileID, '%s\n', addedRxns{:});
                    fclose(fileID);
                else
                    nonOptimizebleMets= [nonOptimizebleMets; biomassMetName{1}];
                end
            catch err
                errors= [errors; err.message];
                nonOptimizebleMets= [nonOptimizebleMets; biomassMetName{1}];
            end   
        end
        model= removeRxns(model, strcat('DM_', biomassMetName{1}));
    end

    %%% deleteing genes which are added through gapfilling second from other
    %%% models
    model.genes(nOrgGenes+1:end)= [];
    model.rxnGeneMat(:, nOrgGenes+1:end)= [];
%     model.geneFrom(nOrgRxns+1:end)= [];
    model.grRules(nOrgRxns+1:end)=repmat({''}, length(model.rxns)-nOrgRxns, 1);
    
    resultsOutputFile= fullfile(gapfillResultsOutputDir, 'nonOptimizebleMets.txt');
    fileID= fopen(resultsOutputFile, 'w');
    fprintf(fileID, '%s\n', nonOptimizebleMets{:});
    fclose(fileID);

end

