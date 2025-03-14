function [modelAfter, finalResults]= eliminateSBCs(model)
    % description: function for removing the Stoichiometric Balanced Cycles
    % (SBCS) / Thermodynamically Infeasible Cycles (TICs)
    % input
    %   model       initial model
    % output
    %   modelAfter      model without SBCs
    %   finalResults    complete results of the function

    finalResults= [];

    model= convertToIrreversible(model);

    biomassIndex= find(model.c);
    [solutionBiomass, ~, ~, ~]= checkBiomassProduction(model, 'one', false, 0);

    exchangeRxns= find(findExcRxns(model));

    orgLBs= model.lb(exchangeRxns);
    orgUBs= model.ub(exchangeRxns);
    
    model.lb(exchangeRxns)= 0;
    model.ub(exchangeRxns)= 0;
    
    SCBs= [];
    for i= 1:numel(model.rxns)
        if isempty(find(i==exchangeRxns))
            model.c= zeros(length(model.rxns),1);
            model.c(i)= 1;
            solution= optimizeCbModel(model);
            vMax= solution.f;
            if vMax > 0
                SCBs= [SCBs; i];    
            end
        end
    end

    model.lb(exchangeRxns)= orgLBs;
    model.ub(exchangeRxns)= orgUBs;

    model.c= zeros(length(model.rxns),1);
    model.c(biomassIndex)= 1;

    if isempty(SCBs)
        disp("there is no SCB in the model")
        return
    else
        disp("Model has SCBs")
    end

    % defining the problem
    nRxns= length(model.rxns);
    nMets= length(model.mets);
    nSCBs= length(SCBs);


    % variables
    % x= |v v' delta| shape: (nRxns + nRxns + nSCBs) * 1
    % f= |0/1 0  0|     shape: (nRxns + nRxns + nSCBs) * 1

    fSCBs= zeros(nRxns, 1);
    for i= 1:nSCBs
        index= SCBs(i);
        fSCBs(index, 1)= 1;
    end

    f= [fSCBs; zeros(nRxns + nSCBs, 1)];

    % inequality constraints: Aeq @ x <= beq
    % Aeq= |S 0 0|
    %      |0 S 0|
    %      |0 0 1| shape: (2nMets + 1) * (2nRxns + nSCBS)
    % beq= |0 0 l|   shape: (2nMets + 1) * 1
    Aeq= [model.S, zeros(nMets, nRxns + nSCBs);
          zeros(nMets, nRxns), model.S, zeros(nMets, nSCBs);
          zeros(1, 2 * nRxns), ones(1, nSCBs)];

    ASCB= zeros(nSCBs, nRxns);
    vLBSCB= zeros(nSCBs, nSCBs); 
    vUBSCB= zeros(nSCBs, nSCBs);
    for i= 1:nSCBs
        index= SCBs(i);
        ASCB(i, index)= 1;
        vLBSCB(i, i)= model.lb(index);
        vUBSCB(i, i)= model.ub(index);
    end



    % Aub= |-ASCB 0      vLBSCB|
    %      |ASCB  0      -vUBSCB|
    %      |0     -ASCB  vLBSCB|
    %      |0     ASCB   -vUBSCB| shape: 4nSCBs * (2nRxns + nSCBs)
    % bub= |0 0 0 0|     shape: 4nSCB * 1
    Aub= [-ASCB, zeros(nSCBs, nRxns)  vLBSCB;
          ASCB, zeros(nSCBs, nRxns), -vUBSCB;
          zeros(nSCBs, nRxns), -ASCB, vLBSCB;
          zeros(nSCBs, nRxns), ASCB, -vUBSCB];
    bub= zeros(4 * nSCBs, 1);

    

    lb= [model.lb; model.lb; zeros(nSCBs, 1)];
    lb(exchangeRxns)= 0;

    ub= [model.ub; model.ub; ones(nSCBs, 1)];
    ub(exchangeRxns)= 0;

    ngamIndex= find(strcmp(model.rxns, 'NGAM'));
%     lb(ngamIndex)=1;
%     ub(ngamIndex)=1;

    params.FeasibilityTol= 1e-9;
    params.intFeasTol= 1e-9;
    params.OutputFlag= 0;
    params.Threads= 8;

    l= 0;
    biomassThreshold= 1;

%     x0= [zeros(nRxns, 1); solutionBiomass.x; ones(nSCBs, 1)];

    count= 1;
    increment= 1;

    while true

        beq= [zeros(2 * nMets, 1); nSCBs-l];
        lb(biomassIndex+nRxns)= solutionBiomass.f * biomassThreshold;

        problemMILP.obj= f;
        problemMILP.A= [sparse(Aub); sparse(Aeq)];
        problemMILP.vtype= [repmat('C', 2 * nRxns, 1); repmat('B', nSCBs, 1)];
        problemMILP.sense= [repmat('<', size(Aub, 1), 1); repmat('=', size(Aeq, 1), 1)];
        problemMILP.rhs= full([bub(:); beq(:)]);
        problemMILP.modelsense= 'max';
        problemMILP.lb= lb;
        problemMILP.ub= ub;
    
        resultMILP= gurobi(problemMILP, params);

        if strcmp(resultMILP.status, 'OPTIMAL')
            totalFlux= sum(resultMILP.x(1:nRxns));
%             fprintf("Feasible, knocking out %d rxns, biomass percentage: %0.2f total SCB fluxes: %d\n", ...
                l, biomassThreshold, totalFlux);

            D= resultMILP.x(2*nRxns+1:end);
            if totalFlux <= 1e-3
                finalResults(count, 1)= biomassThreshold;
                finalResults(count, 2)= l;
                finalResults(count, 3)= totalFlux;
                
                ToRemoveRxns= [];
                for i= 1:length(D)
                    if D(i) == 0
                        index= SCBs(i);
                        ToRemoveRxns= [ToRemoveRxns; index];
                    end
                end
                modelAfter= removeRxns(model, model.rxns(ToRemoveRxns));

                for i= 1:length(modelAfter.rxns)
                    if endsWith(modelAfter.rxns(i), '_b')
                        baseName= strrep(modelAfter.rxns(i), '_b', '');
                        fwRxnId= strcat(baseName, '_f');
                        if ~any(find(strcmp(modelAfter.rxns, fwRxnId)))
                            modelAfter.rxns(i)= baseName;
                            modelAfter.S(:, i)= -modelAfter.S(:, i);
                            tempLb= model.lb(i);
                            modelAfter.lb(i)= -modelAfter.ub(i);
                            modelAfter.ub(i)= tempLb;
                        end
                    end
                end

                modelAfter= convertToReversible(modelAfter);
                [solutionAfter, ~, ~, ~]= checkBiomassProduction(modelAfter);
                
                save(fullfile('..', 'data', 'SCBs', 'models', ...
                    strcat(model.description, '_', ...
                    num2str(biomassThreshold), '_', ...
                    num2str(l), '.mat')), ...
                    'modelAfter')

                biomassThreshold = biomassThreshold - 0.1;
                if biomassThreshold < 0.95
                    break
                end
                fprintf("Biomass threshold reduced to %0.2f\n", biomassThreshold);
                l= 0;
                count = count+1;
                ub= [model.ub; model.ub; ones(nSCBs, 1)];
                ub(exchangeRxns)= 0;
                lb(ngamIndex)=1;
                ub(ngamIndex)=1;

                

            else
                l= l + increment;
                lastResult= resultMILP;
    
                knockedOutIndices= find(D == 0);
                ub(knockedOutIndices + 2*nRxns) = 0;
            end

        elseif strcmp(resultMILP.status, 'INFEASIBLE')
            fprintf("Infeasible with knocking out %d reactions with biomass percentage %0.2f\n",...
                l, biomassThreshold)
            totalFlux= sum(lastResult.x(1:nRxns));
            fprintf("Total flux of model without uptake: %f\n", totalFlux);           

            finalResults(count, 1)= biomassThreshold;
            finalResults(count, 2)= l- increment;
            finalResults(count, 3)= totalFlux;
            
            biomassThreshold = biomassThreshold - 0.1;
            if biomassThreshold < 0.95
                break
            end
            fprintf("Model is infeasible, biomass threshold reduced to %0.2f\n", biomassThreshold);
            l= 0;
            count = count+1;
            ub= [model.ub; model.ub; ones(nSCBs, 1)];
            ub(exchangeRxns)= 0;
            lb(ngamIndex)=1;
            ub(ngamIndex)=1;

        end
        
    end
end


