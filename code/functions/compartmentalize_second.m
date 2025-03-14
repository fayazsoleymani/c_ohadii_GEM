function model= compartmentalize_second(model, modelType, lambda, ...
    countLowerLimit, predefinedRxnComps)
    % description: this function solves the compartmentalization 
    % optimization problem and deletes the redundant reactions
    % inputs
    %   model               input model
    %   modelType           'one'/'pan'
    %   lambda              weight for the replicate rxns in the objective
    %                       function (default= 0.1)
    %   countLowerLimit     minimum number of replicates for the reactions
    %   predefinedRxnComps  predefined compartments for reactions if exist
    % output
    %   model               compartmentalized model    
    
    if nargin <3
        lambda= 0.1;
    end
    if nargin < 4
        countLowerLimit= 1;
    end
    

    defaultSolution= checkBiomassProduction(model, modelType, false);
    defualtBiomassValue= defaultSolution.f;
    initialFluxes= defaultSolution.x;
    biomassIndices= findBiomassIndices(model);
    if modelType == 'one'
        biomassIndex= biomassIndices(1);
        model.lb(biomassIndex)= defualtBiomassValue;
    end


    % defining some useful counts 
    nRxns= length(model.rxns);
    nMets= length(model.mets);
    nDupRxns= length(find(cellfun(@(x) contains(x, 'dRXN'), model.rxns)));
    nTransRxns= length(find(cellfun(@(x) contains(x, 'tRXN'), model.rxns)));
    nOrgRxns= length(model.rxns)- nDupRxns- nTransRxns;
    

    % finding instances of duplicate rxns and save them in a dictionary
    dupInstancesDict= dictionary(string([]), {});
    dupRxnsIndicesModel= [];
    dupRxnsIndicesModel2D= dictionary();
    D2ModelIndex= dictionary();
    
    counter=1;
    for i= 1:size(model.rxns, 1)
        if startsWith(model.rxns(i), 'dRXN')
            dupRxnsIndicesModel= [dupRxnsIndicesModel, i];
            dupRxnsIndicesModel2D(i)= counter;
            D2ModelIndex(counter)= i;
            splited= strsplit(model.rxns{i}, '_');
            mainRxn= strjoin(splited(2:end-1), '_');
            if ~isKey(dupInstancesDict, mainRxn)
                dupInstancesDict(mainRxn)= {i};
            else
                temp= dupInstancesDict(mainRxn);
                temp= mat2cell([cell2mat(temp), i], 1);
                dupInstancesDict(mainRxn)= temp;
            end
            counter= counter+1;
        end
    end
    nGp= size(keys(dupInstancesDict), 1);

    
    % defining variables for minimum and maximum flux for duplicate and transport rxns.


    % objective function
    f= [zeros(nRxns, 1); lambda*ones(nDupRxns, 1); ones(nTransRxns, 1)];
    
    % equality constraints: Aeq@x == beq
    % - Aeq= [S 0 0]    shape: (nMets , (nRxns + nDups + nTrans))
    % - x= [v; d; t]      shape: (nRxns+ nDups + nTrans, 1) 
    % - beq= [0]        shape: (nMets, 1)
    Aeq= [model.S, zeros(nMets, nDupRxns+nTransRxns)];
    beq= model.b;

    % creating array which defines a constraint for number of instances for
    % a rxn that can be present in the final model.     shape: (nGp, nDup)
    dupRxns= keys(dupInstancesDict);
    AGPMin= zeros(nDupRxns, nOrgRxns+nDupRxns);
    AGPMax= zeros(nDupRxns, nDupRxns+nDupRxns);
    AGPCount= zeros(nGp, nDupRxns);

    vDupMin= zeros(nDupRxns);
    vDupMax= zeros(nDupRxns);
    vTransMin= model.lb(nOrgRxns+nDupRxns+1:end);
    vTransMax= model.ub(nOrgRxns+nDupRxns+1:end);

    for groupIndex= 1:nGp
        rxn= dupRxns{groupIndex};
        indicesModel= cell2mat(dupInstancesDict(rxn));
        for i= 1:length(indicesModel)
            rxnIndex= indicesModel(i);
            indexD= dupRxnsIndicesModel2D(rxnIndex);
            AGPCount(groupIndex, indexD)= -1;
            AGPMin(indexD, rxnIndex)= -1;
            AGPMax(indexD, rxnIndex)= 1;
            vDupMin(indexD, indexD)= model.lb(rxnIndex);
            vDupMax(indexD, indexD)= -model.ub(rxnIndex);
        end
    end
    
    % inequality constraint: A@x <= b
    % A= |AGPMin  0  v_i^min  0      |    
    %    |AGPMax  0 -v_i^max  0      |
    %    |0   0  -I  0        v_i^min|
    %    |0   0   I  0       -v_i^max|
    %    |0   0   0  AGP      0      |
    %               shape: (2nDup + 2nTrans + nGP, nOrg+nDup+nTran+nDup+nTrans)
    % x= [v; d; t]  shape: (nRxns+ nDups + nTrans, 1) 
    % b= [0; -1]    shape: (2nDup + 2nTrans + nGp, 1)

    A= [AGPMin, zeros(nDupRxns, nTransRxns), vDupMin, zeros(nDupRxns, nTransRxns);
        AGPMax, zeros(nDupRxns, nTransRxns), vDupMax, zeros(nDupRxns, nTransRxns);
        zeros(nTransRxns, nOrgRxns+nDupRxns), -eye(nTransRxns), zeros(nTransRxns, nDupRxns), diag(vTransMin);
        zeros(nTransRxns, nOrgRxns+nDupRxns), eye(nTransRxns), zeros(nTransRxns, nDupRxns), -diag(vTransMax);
        zeros(nGp, nRxns), AGPCount, zeros(nGp, nTransRxns)];
    b= [zeros(2*(nDupRxns+ nTransRxns), 1); -countLowerLimit*ones(nGp, 1)];

    intcon= nRxns+1:nRxns+nDupRxns+nTransRxns;
    
    % defining lower bounds and upper bounds
    lb= [model.lb; zeros(nDupRxns+nTransRxns, 1)];
    
    % I have added the condition for Marius model
    if length(find(model.c)) == 1
        lb(biomassIndex)= defualtBiomassValue;
    end
    
    predefinedRxns= predefinedRxnComps.rxn;
    predefinedComps= predefinedRxnComps.comp;

    for i=1:numel(predefinedRxns)
        predefinedRxn= predefinedRxns{i};
        if find(strcmp(predefinedRxn, keys(dupInstancesDict)))
            comp= predefinedComps{i};
            compIndexModel= find(strcmp(comp, model.comps));
            compSuffix= upper(model.compNames(compIndexModel));
            predefinedRxnModelId= strcat('dRXN_', predefinedRxn, '_', compSuffix);
            predefinedRxnModelIndex= find(strcmp(predefinedRxnModelId, model.rxns));
            predefinedRxnsIndexD= dupRxnsIndicesModel2D(predefinedRxnModelIndex);
            lb(nRxns+predefinedRxnsIndexD)= 1;
        end
    end
    
        
    ub= [model.ub; ones(nDupRxns+nTransRxns, 1)];
    
    x0= [initialFluxes;ones(nDupRxns+nTransRxns, 1)];
    
%     statusOfx0= checkFeasibility(x0, A, b, Aeq, beq, 1e-5);
    
    params.FeasibilityTol= 1e-9;
    params.intFeasTol= 1e-9;
    
    %%% first solving the problem as a LP with Gurobi
    problemLP.obj= f;
    problemLP.A= [sparse(A); sparse(Aeq)];
    problemLP.vtype= repmat('C', nRxns, 1);
    problemLP.vtype(intcon)= 'C';
    problemLP.sense= [repmat('<', size(A, 1), 1); repmat('=', size(Aeq, 1), 1)];
    problemLP.rhs= full([b(:); beq(:)]);
    problemLP.modelsense= 'min';
    problemLP.start= x0;
    problemLP.lb= lb;
    problemLP.ub= ub;
    resultLP= gurobi(problemLP, params);
    
    % statusOfLPsolution= checkFeasibility(resultLP.x, A, b, Aeq, beq, 1e-5);
    
    
    %%% solving the MILP with Gurobi
    problemMILP.obj= f;
    problemMILP.A= [sparse(A); sparse(Aeq)];
    problemMILP.vtype= repmat('C', nRxns, 1);
    problemMILP.sense= [repmat('<', size(A, 1), 1); repmat('=', size(Aeq, 1), 1)];
    problemMILP.rhs= full([b(:); beq(:)]);
    problemMILP.modelsense= 'min';
    problemMILP.lb= lb(1:nRxns);
    problemMILP.ub= ub(1:nRxns);
    
    % from the LP solution, d and t variables that are near 0 and 1, with
    % the threshold distance are fixed to 0 and 1 respectively, and the
    % others are considered as binary values.
    threshold_int= 0;
    roundedToOne= 0; roundedToZero= 0; remainedCon= 0;
    for i=nRxns+1:nRxns+nDupRxns+nTransRxns
        if resultLP.x(i)>= 1- threshold_int
            problemMILP.vtype(i)= 'C';
            problemMILP.lb(i)= 1;
            problemMILP.ub(i)= 1;
            roundedToOne= roundedToOne+1;
        elseif resultLP.x(i) <= threshold_int
            problemMILP.vtype(i)= 'C';
            problemMILP.lb(i)= 0;
            problemMILP.ub(i)= 0;
            roundedToZero= roundedToZero+1;
        else
            problemMILP.vtype(i)= 'B';
            problemMILP.lb(i)= 0;
            problemMILP.ub(i)= 1;
            remainedCon= remainedCon+1;
        end
    end
    
    fprintf("Rounded to one:  %d\tRounded to zero:  %d\tRemained:  %d\n", ...
        roundedToOne, roundedToZero, remainedCon); 
    
    resultMILP= gurobi(problemMILP, params);

    model.lb(biomassIndex)= 0;

    D= resultMILP.x(nRxns+1:nRxns+nDupRxns);
    T= resultMILP.x(nRxns+nDupRxns+1:end);
    
    fprintf("Duplicate:  %d out of %d\tTransports:  %d out of %d\n", ...
        length(find(D)), nDupRxns, length(find(T)), nTransRxns);
    
    ToRemoveRxns= [];
    for i= 1:length(D)
        if D(i) == 0
            ToRemoveRxns= [ToRemoveRxns; D2ModelIndex(i)];
        end
    end
    for i= 1:length(T)
        if T(i) == 0
            ToRemoveRxns= [ToRemoveRxns; nOrgRxns+nDupRxns+i];
        end
    end

    % concat the gpur of the rxns that are going to be removed, to the
    % duplicate versions of the rxn
    for i= 1: length(ToRemoveRxns)
        rxnIndex= ToRemoveRxns(i);
        if ~isempty(model.grRules{rxnIndex})
%             disp(rxnIndex);
            splited= strsplit(model.rxns{rxnIndex}, '_');
            mainRxn= strjoin(splited(2:end-1), '_');
            otherInstances= dupInstancesDict(mainRxn);
            otherInstances= otherInstances{1};
            index= find(otherInstances==rxnIndex);
            otherInstances(index) = [];
            if length(otherInstances) == 1
                model.rxns{otherInstances} = mainRxn;
            end

            dupInstancesDict(mainRxn) = {otherInstances};
            firstRemained= otherInstances(1);
            removedGenes= cellfun(@strtrim, split(model.grRules(rxnIndex), 'or'),...
                'UniformOutput', false);
            mainGenes= cellfun(@strtrim, split(model.grRules(firstRemained), 'or'),...
                'UniformOutput', false);

            if modelType== 'pan'
                model.grRules{firstRemained}= strjoin([mainGenes;removedGenes], ' or ');
                model.grRules{rxnIndex}= ''; 
            else
                model= changeGrRules(model, model.rxns(firstRemained), ...
                    strjoin(unique([mainGenes;removedGenes]), ' or '));
                model= changeGrRules(model, model.rxns(rxnIndex), {''});
            end
        end
    end
    
    
    %%% forth argument: my model: true, marius: False
    model= removeReactions(model, ToRemoveRxns, true, true, true);
    
    model= rmfield(model, "csense");
    for metIndex= 1:length(model.metNames)
        model.metNames{metIndex}= strcat(model.metNames{metIndex}, ...
            ' [', model.compNames{model.metComps(metIndex)}, ']');

    end

    solutionAfter= checkBiomassProduction(model, 'one', false);

end