function info= gemAnalyse(model, operateEssenGenes, operateEssenRxns)
    % description: This is a function for analying the genome-scale
    % metabolic model
    % input
    %   model               input model
    %   operateEssenGenes   binary argument to operate the gene
    %                       essentiallity analyses or not
    %   operateEssenRxns    binary argument to operate the rxn essentiality 
    %                       analysis or not 
    % output
    %   info                results of the analysis, containing fields:
    %                       nGenes, nMets, nUniqueMets, deadendMets,
    %                       metComps, nRxns, nEnzymaticRxns, nRevRxns,
    %                       modelIrrev, blocked rxns, activeGenes,
    %                       rxnComps, nTransRxns, enzymaticTransportRxns,
    %                       exchangeRxns, media, essentialGenes,
    %                       essentialRxns
    
    if nargin < 2
        operateEssenGenes = false;
    end
    if nargin < 3
        operateEssenRxns = false;
    end

    info= struct();

    biomassIndex=findBiomassIndices(model);
    solution= checkBiomassProduction(model, 'one', false, true);

    if ~isfield(info, 'nGenes')
        info.nGenes= length(model.genes);
        fprintf("Number of Genes:  %d\n", info.nGenes)
    end
    disp(repmat('-', 1, 40));
    
    if ~isfield(info, 'nMets')
        info.nMets= length(model.mets);
        fprintf("Number of Mets:  %d\n", info.nMets);
    end

    if ~isfield(info, 'nUniqueMets')
        clearedCompsModel= clearComps(model, true, true, false);
        uniqueMets= unique(clearedCompsModel.mets);
        info.nUniqueMets= length(uniqueMets);
        fprintf("Number of Unique Mets:  %d\n", info.nUniqueMets);
    end

    % dead-end metabolites
    if ~isfield(info, 'deadendMets')
        info.deadendMets= zeros(length(model.mets), 1);
    end
    for metIndex= 1:length(model.mets)
        sRow= model.S(metIndex, :);
        if all(sRow >= 0) || all(sRow <= 0)
            info.deadendMets(metIndex) = 1;
        end
    end
    
    fprintf("Number of dead-end Metabolites:  %d\n", length(find(info.deadendMets)));

    % met compartments
    if ~isfield(info, 'metComps')
        info.metComps= histcounts(model.metComps);
        disp("Metabolites Compartments:")
        for i= 1:length(info.metComps)
            fprintf("\t%s:  %d\n", model.compNames{i}, info.metComps(i));
        end
    end
    disp(repmat('-', 1, 40));


    if ~isfield(info, 'nRxns')
        info.nRxns= length(model.rxns);
        fprintf("Number of Rxns:  %d\n", info.nRxns);
    end

    % rxns with gpr
    if ~isfield(info, 'nEnzymaticRxns')
        if isfield(model, 'grRules')
            enzymaticRxns= ~cellfun(@isempty, model.grRules);
        elseif isfield(model, 'rules')
            enzymaticRxns= ~cellfun(@isempty, model.rules);
        end
        info.enzymaticRxns= model.rxns(enzymaticRxns);
        info.nEnzymaticRxns= sum(enzymaticRxns);
        fprintf("Number of enzymatic Rxns:  %d\t\tOthers:  %d\n", ...
            info.nEnzymaticRxns, length(model.rxns)- info.nEnzymaticRxns);
    end
    
    if ~isfield(info, 'nRevRxns')
        if isfield(model, 'rev')
            info.nRevRxns= length(find(model.rev==1));
        else
            tempRev= zeros(length(model.rxns), 1);
            for i=1:length(model.rxns)
                if model.lb(i)<0
                    tempRev(i)= 1;
                end
            end
            info.nRevRxns= length(find(tempRev==1));
        end
        fprintf("Number of rev reactions:  %d\t\tNumber of Irrev reactions:  %d\n", ...
            info.nRevRxns, length(model.rxns)- info.nRevRxns);
    end

    modelIrrev= convertToIrreversible(model);
    if ~isfield(model, 'modelIrrev')
        info.modelIrrev= modelIrrev;
    end




    % blocked rxns
    if ~isfield(info, 'noFlux')
        [minFlux, maxFlux] = fluxVariability(modelIrrev);
        info.noFlux= find(minFlux==0 & maxFlux == 0);
%         fprintf("Number of active Rxns:  %d\n", ...
%             length(modelIrrev.rxns)- length(info.noFlux))



        activeRxnsRaw= setdiff(modelIrrev.rxns, modelIrrev.rxns(info.noFlux));
%         activeEnzymaticRxns= intersect(info.enzymaticRxns, info.activeRxns);
%         fprintf("Number of active enzymatic Rxns:  %d\n", ...
%             length(activeEnzymaticRxns));
        activeRxns= {};
        for i=1:length(activeRxnsRaw)
            forRxn= activeRxnsRaw(i);
            if endsWith(forRxn, '_f')
                forRxnIndex= find(strcmp(modelIrrev.rxns, forRxn));
                maxFluxForward= maxFlux(forRxnIndex);
                backRxn= strrep(forRxn, '_f', '_b');
                backRxnIndex= find(strcmp(modelIrrev.rxns, backRxn));
                j= find(strcmp(activeRxns, backRxn));
                maxFluxBackward= maxFlux(backRxnIndex);
                if maxFluxForward ~= maxFluxBackward
                    activeRxns= [activeRxns;forRxn];
                end
            elseif endsWith(forRxn, '_b')
                continue
            else
                activeRxns= [activeRxns;forRxn];
            end
        end
        info.activeRxns= activeRxns;
        fprintf("Number of active Rxns:  %d\n", ...
            length(activeRxns));

    end
    
    if ~isfield(info, 'activeGenes')
        temp= [];
        for i=1:length(info.activeRxns)
            rxnIndex= strcmp(model.rxns, info.activeRxns(i));
            genes= find(model.rxnGeneMat(rxnIndex, :));
            for j=1:length(genes)
                temp= [temp;genes(j)];
            end
        end
        info.activeGenes= model.genes(unique(temp));
    end

    % rxn compartmets
    if ~isfield(info, 'rxnComps')
        info.rxnComps= dictionary(1:length(model.compNames), zeros(1, length(model.compNames)));
        info.rxnCompartments= [];
        for rxnIndex= 1:length(model.rxns)
            if rxnIndex == biomassIndex
                continue
            end
            rxnMetIndices= find(model.S(:, rxnIndex));
            metComps= zeros(length(rxnMetIndices), 1);
            for i=1:length(rxnMetIndices)
                metComps(i)= model.metComps(rxnMetIndices(i));
            end
            metComps= unique(metComps);
            if length(metComps) == 1
                info.rxnComps(metComps)= info.rxnComps(metComps) + 1;
                info.rxnCompartments= [info.rxnCompartments;metComps];
            elseif length(metComps)>=1
                info.rxnCompartments= [info.rxnCompartments ; 0];
                if ~isfield(info, 'nTransRxns')
                    info.nTransRxns= 1;
                else
                    info.nTransRxns= info.nTransRxns + 1;
                end
                if ~isempty(model.grRules{rxnIndex})
                    if ~isfield(info, 'enzymaticTransportRxns')
                        info.enzymaticTransportRxns= [rxnIndex];
                    else
                        info.enzymaticTransportRxns= [info.enzymaticTransportRxns; rxnIndex];
                    end
                end
            end
        end
        disp("Reaction Compartments");
        for i=1:length(keys(info.rxnComps))
            compName= model.compNames{i};
            compCount= info.rxnComps(i);
            fprintf("\t%s-->  %d\n", compName, compCount);
        end
        fprintf("Number of Transport Rxns:  %d\n", info.nTransRxns);
    end
    
    % exchange  and demand rxns
    if ~isfield(info, 'exchangeRxns')
        [info.exchangeRxns, info.exchangeRxnIndices]= getExchangeRxns(model); %[exchangeRxns,exchangeRxnsIndexes]  
    end

    % checkMedium
    if ~isfield(info, 'media')
        info.exchangeRxns= getMedia(model);
    end


    disp(repmat('-', 1, 40));


    % essential genes
    if operateEssenGenes
        if ~isfield(info, 'essentialGenes')
            essentialGenes= {};
            for i= 1:numel(model.genes)
                gene= model.genes{i};
                inRxns= find(model.rxnGeneMat(:, i));
                for j= 1:numel(inRxns)
                    rxnIndex= inRxns(j);
                    orgLb= model.lb(rxnIndex);
                    orgUb= model.ub(rxnIndex);
                    model.lb(rxnIndex)= 0;
                    model.ub(rxnIndex)= 0;
                    solution= checkBiomassProduction(model);
                    model.lb(rxnIndex)= orgLb;
                    model.ub(rxnIndex)= orgUb;
                    if isnan(solution.f) || solution.f < 10^-4
                        essentialGenes= [essentialGenes; gene];
                        break;
                    end
                end
            end
            info.essentialGenes= unique(essentialGenes);
            fprintf("Number of Essential Genes:  %d\n", length(essentialGenes));
        end
    end

    
    % essential reactions
    if operateEssenRxns
        if ~isfield(info, 'essentialRxns')
            essentialRxns= {};
            for i= 1:numel(model.rxns)
                rxn= model.rxns{i};
                orgLb= model.lb(i);
                orgUb= model.ub(i);
                model.lb(i)= 0;
                model.ub(i)= 0;
                solution= checkBiomassProduction(model);
                model.lb(i)= orgLb;
                model.ub(i)= orgUb;
                if isnan(solution.f) || solution.f < 10^-4
                    essentialRxns= [essentialRxns; rxn];
                end
            end
            info.essentialRxns= essentialRxns;
            fprintf("Number of Essential Rxns:  %d\n", length(essentialRxns));
        end
    end
   
    %%% biomass production
    info.solution= checkBiomassProduction(model, 'one', false);
    disp(repmat('=', 1, 50));


end