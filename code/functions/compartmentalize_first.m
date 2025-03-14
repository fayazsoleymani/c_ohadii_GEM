function model= compartmentalize_first(model, GSS, threshold, modelType, ...
    predefinedRxnComps)
    % descripttion: This function uses GSS struct which has been created
    % with constructGSS function to operate the first step of the
    % compartmentalization. This process is based on deciding compartment
    % for each reaction based on the majority voting of the localization
    % predictions of its genes. Compartments which have more than threshold
    % of the votes as the threshold are considered, and as a result there 
    % are some rxns which are in more than one compartment and named with
    % prefix of "dRXN_". At the end for the metabolites that are present in
    % more than one compartment, all of the possible transport rxns are
    % added to the model which have prefix of tRXN.

    % inputs:
    %   model:              uncompartmentalized model 
    %   GSS:                GSS file created with construct GSS function
    %   threshold           the threshold of the votes for adding reaction to
    %                       the compartment (default 0.25)
    %   modelType           one/pan (default 'one')
    %   predefinedRxnComps  file contained the predefined compartments for
    %                       rxns
    %   
    % output:
    %   model               compartmentalized model
    
    if nargin <3
        threshold= 0.25;
    end
    if nargin < 4
        modelType= 'one';
    end
    if nargin < 5
        predefinedRxnComps.rxns= {''};
        predefinedRxnComps.comp= {''};
    end

    while length(model.rxns) ~= length(model.rxnECNumbers)
        model.rxnECNumbers(end+1)= {''};
    end

    % creating fields comps and compNames for the model based on GSS file
    model.compNames= GSS.compartments';
    model.comps= cell(size(model.compNames));
    for compIndex= 1:numel(model.comps)
        model.comps{compIndex}= model.compNames{compIndex}(1:2);
    end

    % displaying the number of mets and rxns in the original model
    fprintf("info of the orignal model:\nnRXN  %d\tnMets:  %d\n", ...
    length(model.rxns), length(model.mets));

    % finding biomass index
    biomassIndex= findBiomassIndices(model);
    % creating some arrays to keep track of the process
    originalMets= model.mets;
    rxnComps= cell(size(model.rxns, 1), 1);
    metHasMoved= zeros(size(model.mets, 1), 1);

    nDupRxns= 0;
    
    % looping through the rxns
    for rxnIndex= 1:numel(model.rxns)
        % the biomass reaction is remained in the default compartment which
        % is the first compartment in the GSS.
        if any(biomassIndex(:)==rxnIndex)
            continue
        end
        
        % grabing the genes of the rxn from rxnGeneMat field of the model
        % if exists, otherwise from the grRules field.
        if isfield(model, 'rxnGeneMat')
            genes= model.genes(find(model.rxnGeneMat(rxnIndex, :)));
        else
            genes= cellfun(@strtrim, split(model.grRules(rxnIndex), 'or'), ...
                UniformOutput=false);
        end
        % creating the dictionary to count the compartment votes
        % initialized with zero.
        compsDict= dictionary(1:numel(GSS.compartments), zeros(size(GSS.compartments)));
        comp2genes= dictionary(double([]), cell([]));
        if ~strcmp(model.grRules{rxnIndex}, '') && ...
                ~strcmp(model.grRules{rxnIndex}, 'AT2G05990') && ...
                ~strcmp(model.grRules{rxnIndex}, 'AT4G30580') && ...
                ~strcmp(model.grRules{rxnIndex}, 'CHLREDRAFT_32523')
            for i= 1:numel(genes)
                gene= genes{i};
                geneIndexGSS= find(strcmp(GSS.genes, gene));
                comp= GSS.geneComp(geneIndexGSS);
                compsDict(comp)= compsDict(comp) + 1;

                if ~isKey(comp2genes, comp)
                    comp2genes(comp) = {geneIndexGSS};
                else
                    temp= comp2genes(comp);
                    temp= [temp{1}, geneIndexGSS];
                    comp2genes(comp) = {temp};
                end
            end
            
            % sorting the predicted values
            probValues = values(compsDict);
            [sortedValues, sortedIndices] = sort(probValues, 'descend');
            % taking the most probable compartment and for others if the
            % votes are more than the threshold.
            rxnComps{rxnIndex}= [sortedIndices(1)];
            for compIndex=2:numel(sortedIndices)
                if sortedValues(compIndex)/length(genes) >= threshold
                    rxnComps{rxnIndex}= [rxnComps{rxnIndex}; sortedIndices(compIndex)];
                end
            end
        else
            % if the reaction is non enzymatic, its compartment would be
            % the default compartment which is the first compartment in GSS
            rxnComps{rxnIndex}= 1;
        end

        % checking weather it the reactions should be fixed in a specific
        % compartment
        if ismember(model.rxns{rxnIndex}, predefinedRxnComps.rxn)
            
            predefinedRxn= model.rxns{rxnIndex};
            memberIndex= find(strcmp(predefinedRxn, predefinedRxnComps.rxn));
            predefinedComp=predefinedRxnComps.comp(memberIndex);
            predefinedCompIndex= find(strcmp(predefinedComp,model.comps));

            if ismember(predefinedCompIndex, rxnComps{rxnIndex})
                index= find(rxnComps{rxnIndex} == predefinedCompIndex);
                rxnComps{rxnIndex} = [];
            end

            if startsWith(predefinedRxn, 'EXC_') || startsWith(predefinedRxn, 'Ex_')
                rxnComps{rxnIndex}= predefinedCompIndex;
            else
                rxnComps{rxnIndex}= [predefinedCompIndex; rxnComps{rxnIndex}];
            end     
        end
    
        % moving the reaction to its first predicted compartment
        compToMove= rxnComps{rxnIndex}(1);
        [model, metHasMoved, originalMets]= moveRxnToComp(model, rxnIndex, compToMove, metHasMoved, originalMets);
    
        % for other compartments add duplicate reactions
        movedGPRs= {};
        for compIndex= 2:size(rxnComps{rxnIndex}, 1)
            compToMove= rxnComps{rxnIndex}(compIndex);
    
            % creating new reaction and add it to model
            newRxnID= strcat('dRXN', '_', model.rxns(rxnIndex), '_', upper(model.compNames(compToMove)));
            newRxnName= strcat(model.rxnNames(rxnIndex), '_', model.compNames(compToMove));
            
            if ~isempty(genes)
                gssIndices= comp2genes(compToMove);
                gssIndices= gssIndices{1};
                genesOfComp= cell(length(gssIndices), 1);
                for i= 1:length(gssIndices)
                    genesOfComp(i)= GSS.genes(gssIndices(i));
                end
                newRxnGeneRule= strjoin(genesOfComp, ' or ');
                movedGPRs= [movedGPRs; genesOfComp];
            else
                newRxnGeneRule= '';
            end
            [model, ~] = addReaction(model, newRxnID{1},...
                'reactionName', newRxnName{1},...
                'metaboliteList', model.mets(find(model.S(:, rxnIndex))),...
                'stoichCoeffList', full(model.S(find(model.S(:, rxnIndex)), rxnIndex)),...
                'reversible', model.rev(rxnIndex),...
                'lowerBound', model.lb(rxnIndex),...
                'upperBound', model.ub(rxnIndex),...
                'geneRule', newRxnGeneRule);
            if isfield(model, 'rxnECNumbers')
                model.rxnECNumbers(end)= model.rxnECNumbers(rxnIndex);
            end
            if isfield(model, 'rxnReferences')
                model.rxnECNumbers(end)= model.rxnReferences(rxnIndex);
            end
            if isfield(model, 'rxnKEGGID')
                model.rxnKEGGID(end)= model.rxnKEGGID(rxnIndex);
            end
            if isfield(model, 'rxnRheaID')
                model.rxnRheaID(end)= model.rxnRheaID(rxnIndex);
            end
            if isfield(model, 'subSystems')
                model.subSystems(end)= model.subSystems(rxnIndex);
            end
            if model.lb(end) < 0 && model.rev(end) == 0
                model.rev(end)= 1;
            end
            
            [model, metHasMoved, originalMets]= ...
                moveRxnToComp(model, size(model.rxns, 1), compToMove, metHasMoved, originalMets);
            nDupRxns= nDupRxns+1;
        end
        if size(rxnComps{rxnIndex}, 1) > 1
            model.rxns{rxnIndex}= strcat('dRXN_', model.rxns{rxnIndex}, '_', ...
                upper(model.compNames{rxnComps{rxnIndex}(1)}));
            remainedGPRs= setdiff(genes, movedGPRs);
            
            if modelType == 'pan'
                model.grRules{rxnIndex} = strjoin(remainedGPRs, ' or ');
            else
                if ~isempty(remainedGPRs)
                    model= changeGrRules(model, model.rxns(rxnIndex), strjoin(remainedGPRs, ' or '));
                end
            end
        end
    end
    fprintf("Total number of added duplicate rxns:  %d\n", nDupRxns);
        
    nTransRxns= 0;
    % adding transport rxns
    for metIndex= 1:length(model.mets)
        met= originalMets(metIndex);
        % finding the same metabolite in other compartments
        metsOtherComps= setdiff(find(strcmp(originalMets, met)), metIndex);
        if ~isempty(metsOtherComps)
            % loopint through the duplicated metabolites and adding a
            % forward rxn to them.
            for i= 1:length(metsOtherComps)
                metOtherCompIndex= metsOtherComps(i);
                firstMetName= model.mets(metIndex);
                metOtherCompName= model.mets(metOtherCompIndex);
                transportRxnID= strcat('tRXN', '_', met, '_', ...
                    model.comps(model.metComps(metIndex)), '_', ...
                    model.comps(model.metComps(metOtherCompIndex)));
                transportRxnMetaboliteList= {firstMetName{1}; metOtherCompName{1}};
                [model, ~] = addReaction(model, transportRxnID{1}, ...
                    'metaboliteList', transportRxnMetaboliteList, ...
                    'stoichCoeffList', [-1;1], ...
                    'reversible', false, 'lowerBound', 0, 'upperBound', 1000);
                if length(model.rxns) ~=length(model.rev)
                    model.rev(end+1) = 0;
%                     model.eccodes(end+1)= {''};
                end
                nTransRxns= nTransRxns+1;
            end
        end
    end
    
    
    fprintf("Total number of added transport rxns:  %d\n", nTransRxns);
    disp(repmat('-', 1, 40));

    % short brief of the compartmentalized model
    % number of rxns in each compartment 
    rxnCompsDict= dictionary(1:length(model.comps), zeros(1, length(model.comps)));
    for i=1:length(rxnComps)
        for j= 1:length(rxnComps{i})
            comp= rxnComps{i}(j);
            rxnCompsDict(comp)= rxnCompsDict(comp)+ 1;
        end
    end

    disp("RXN compartments:");
    for compIndex= 1:length(keys(rxnCompsDict))
        fprintf("%s --> %d\n", model.compNames{compIndex}, rxnCompsDict(compIndex))
    end
    disp(repmat('-', 1, 40));

    % number of metabolites in each compartment.
    metCompDict= dictionary(1:length(model.comps), zeros(1, length(model.comps)));
    for metIndex= 1:length(model.metComps)
        comp= model.metComps(metIndex);
        metCompDict(comp)= metCompDict(comp)+1;
    end

 
    disp("Met compartments:");
    for compIndex= 1:length(keys(metCompDict))
        fprintf("%s --> %d\n", model.compNames{compIndex}, metCompDict(compIndex))
    end
    disp(repmat('-', 1, 40));

    %%% nRxns, nMets
    fprintf("info of the compartmentalized model:\nnRxns:  %d\tnMets:  %d\n",...
        length(model.rxns), length(model.mets));

    disp(repmat('=', 1, 40));
end