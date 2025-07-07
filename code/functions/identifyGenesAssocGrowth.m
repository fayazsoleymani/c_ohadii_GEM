function [correlations, corresPvals, targetGenes, modelIrrev]= identifyGenesAssocGrowth(model)
    % description: function for identifying the genes which are associated
    % with high growth rates based on pFBA flux distribution with 100 steps
    % input
    %   model           model of interest
    % output
    %   correlations    pearson correlation between of the flux of all
    %                   reactions and growth 
    %   corresPvals     corresponding p-values
    %   targetGenes     genes participating in rxns whose correlation is
    %                   greater than 0.8
    %   modelIrrev      converted model with all irreversible rxns  

    modelIrrev= convertToIrreversible(model);

    biomassIndex= findBiomassIndices(model);
    maxGrowth= checkBiomassProduction(modelIrrev).f;
    step= maxGrowth/100;
    growthRates=step:step:maxGrowth; 
    fluxes= zeros(length(growthRates), length(modelIrrev.rxns));
    
    for i= 1:length(growthRates)
        growth= growthRates(i);
        model.lb(biomassIndex)= growth;
        model.ub(biomassIndex)= growth;
        [~, ~, ~, MinimizedFlux]= pFBA(model, 'geneoption', 0, 'skipclass', 1);
        fluxes(i, :)= max(MinimizedFlux.x(1:end-1), 0);
    end

%     logTransformedGrowthRates= log(growthRates');
    
    correlations= zeros(length(modelIrrev.rxns), 1);
    pvals= zeros(length(modelIrrev.rxns), 1);
    for i= 1:length(modelIrrev.rxns)
        if all(fluxes(:, i))
%             logTransFluxes= log(fluxes(:, i));
            [rho, pval]= corr(fluxes(:, i), growthRates');
            if ~isnan(rho)
                correlations(i)= abs(rho);
                pvals(i)= pval;
            else
                correlations(i)= 0;
                pvals(i)= 0;
            end
        else
            correlations(i)= 0;
            pvals(i)= 0;
        end
    end
    
   [correlationsSorted, I]= sort(correlations, 1, 'descend');
   corresPvals= pvals(I);


   if ~isfield(model, 'rxnGeneMat')
       modelIrrev= standardizeGPRs(modelIrrev);
   end

   targetRxnsIndices= [];

   index=1;
   while correlationsSorted(index)>0.80
       targetRxnsIndices= [targetRxnsIndices;I(index)];
       index= index + 1;
   end

   targetGeneIndices= [];
   for i=1:length(targetRxnsIndices)
       rxnIndex= targetRxnsIndices(i);
       genes= find(modelIrrev.rxnGeneMat(rxnIndex, :));
       for j= 1:length(genes)
           gene= genes(j);
           targetGeneIndices= [targetGeneIndices; gene];
       end
   end
   targetGeneIndices=unique(targetGeneIndices);
   targetGenes= modelIrrev.genes(targetGeneIndices);

end