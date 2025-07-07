function topSubsystems= subsystemsSum(inputs, model, topN, mode)
    % description: get sorted dictionary of the subsystems of rxns or genes
    % of interest

    % inputs
    %   inputs          array of genes or rxns
    %   model           input GEM          
    %   topN            number of the top subsystems
    %   mode            'genes' or 'rxns'

    % output
    %   topSubsystems   dictionary whose keys are the subsystems and values
    %                   are the number of repetitions

    if strcmp(mode, 'rxns')
        rxns= inputs;
    elseif strcmp(mode,'genes')
        rxns= [];
        genes= inputs;
        for i= 1:length(genes)
            geneIndex= find(strcmp(model.genes, genes{i}));
            catalRxns= find(model.rxnGeneMat(:, geneIndex));
            for j= 1:length(catalRxns)
                rxns= [rxns;model.rxns(catalRxns(j))];
            end
        end
        rxns= unique(rxns);
%         printRxnFormula(model, 'rxnAbbrList', rxns, 'metNameFlag', true);
        fprintf("Number of unique rxns: %d\n", length(rxns))
    end

    allSubsystems= dictionary('', 0);
    for i= 1:length(rxns)
        rxn= rxns{i};
        rxnIndex= find(strcmp(model.rxns, rxn));
        subsystems= model.subSystems{rxnIndex};
        for j= 1:size(subsystems, 2)
            subsystem= subsystems{j};
            if ~isKey(allSubsystems, subsystem)
                allSubsystems(subsystem) = 1;
            else
                allSubsystems(subsystem)= allSubsystems(subsystem) + 1;
            end
        end
    end
    allSubsystems('') = [];
    
    k= keys(allSubsystems);
    v= values(allSubsystems);

    [vSorted, indices]= sort(v, 'descend');
    kSorted= k(indices);

    topSubsystems= dictionary(kSorted(1:topN), vSorted(1:topN));


end