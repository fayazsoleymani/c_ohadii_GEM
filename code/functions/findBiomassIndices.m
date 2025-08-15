function biomassIndices= findBiomassIndices(model, verbose)
    % finding the biomass reaction index #TODO

    if nargin < 2
        verbose= 0;
    end

    if isfield(model, 'description')
        name= model.description;
    elseif isfield(model, 'modelID')
        name= model.modelID;
    elseif isfield(model, 'id')
        name= model.id;
    end
    if verbose
        fprintf("Ref Model: %s\n", name);
    end
    if any(model.c)
        biomassIndices= find(model.c);
    else
        biomassIndices = [];
        for rxnIndex = 1:length(model.rxns)
            if contains(model.rxns{rxnIndex}, 'biomass', 'IgnoreCase', true)
                biomassIndices = [biomassIndices, rxnIndex];
                if verbose
                    fprintf("%d->%s\n", rxnIndex, model.rxns{rxnIndex})
                end
            end
        end
    end
end