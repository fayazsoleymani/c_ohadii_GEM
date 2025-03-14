function model= modifyExchangeRxns(model, ex2lbDict)
    % first, removeing all the exchange reactions that are present in the
    % model and second adding exchange reaction to the model based on the
    % provided exchange metabolites and corresponding lower bound of their
    % uptake
    % input:
    %       model       input model
    %       ex2lbDict   dictionary containing the uptake metabolites and
    %                   corresponding lower bound
    % output:
    %       model       final model

    [exchangeRxns,~]=getExchangeRxns(model);
    model= removeReactions(model, exchangeRxns);

    ex2lbDictKeys= keys(ex2lbDict);
    for i= 1:length(ex2lbDictKeys)
        met= ex2lbDictKeys(i);
        [model, ~]= addExchangeRxns(model,'both',met);
        model.lb(end)= ex2lbDict(met);
    end
    
end