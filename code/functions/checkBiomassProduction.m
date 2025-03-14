function [solution, optimizableMets, nonOptimizebleMets, demandRxnsFluxes]= ...
    checkBiomassProduction(model, modelType, checkAllBiomassMets, verbose)
    % This function checks if the model can produce the biomass
    % inputs
    %       model                   input model
    %       modelType               pan or one
    %       checkAllBiomassMets     check which biomass precursors
    %                               can/can't be produced 
    %       verbose                 print
    % outputs
    %       solution                LP solution
    %       optimizableMets         producable biomass precursors
    %       nonOptimizebleMets      not producable biomass precursors
    %       demandRxnsFluxes        flux of the demand reaction of biomass
    %                               precursor                   

    if nargin < 2
        modelType= 'one';
    else
        modelType= modelType;
    end
    if nargin <3
        checkAllBiomassMets = false;
    else
        checkAllBiomassMets= checkAllBiomassMets;
    end

    if nargin <4
        verbose= 0;
    end

    if modelType == 'pan'
        biomassIndex= find(model.c);
        biomassMetIndices= zeros(length(biomassIndex), 1);
        biomassMetNames= cell(length(biomassIndex), 1);
        for i=1:length(biomassIndex)
            rxnIndex= biomassIndex(i);
            biomassMetIndices(i)= find(model.S(:, rxnIndex));
            biomassMetNames{i}= model.mets(find(model.S(:, rxnIndex)));
        end

    else
        biomassIndex= findBiomassIndices(model);
        if length(biomassIndex)>1
            biomassIndex= biomassIndex(1);
        end
        biomassMetIndices= find(model.S(:, biomassIndex));
        biomassMetNames= model.mets(biomassMetIndices);
    end

    optimizableMets= {};
    nonOptimizebleMets= {};
    demandRxnsFluxes= zeros(length(biomassMetIndices), 1);
    if checkAllBiomassMets
        model.c= zeros(size(model.rxns));
        for metIndex= 1:length(biomassMetIndices)
            biomassMetIndex= biomassMetIndices(metIndex);
            biomassMetName= model.mets(biomassMetIndex);
            model = addDemandReaction(model, biomassMetName);
            model.lb(end)= 1e-3;
            model.c(end)= 1;
            solution = optimizeCbModel(model);
            if solution.stat == 1
                optimizableMets= [optimizableMets; biomassMetName{1}];
            else
                nonOptimizebleMets= [nonOptimizebleMets; biomassMetName{1}];
            end
            demandRxnsFluxes(metIndex)= solution.x(end);   
            model= removeRxns(model, strcat('DM_', biomassMetName{1}));
        end
    end


    model.c(biomassIndex)= 1;
    solution= optimizeCbModel(model);
    if verbose
        if solution.stat== 1 && solution.f > 0
            fprintf("Model can produce biomass with the value equal to %5.4f\n", solution.f)
        else
            fprintf("Model can't produce Biomass\nStatus: %d\tBiomass Value: %f\n", solution.stat, solution.f)
        end
        disp(repmat('-', 1, 55));
    end
end