function model= clearComps(model, mergeComps, clearMets, clearMetNames)
    % the description should be updated
    % description: remove all the compartment  suffixes from model.met and
    % model.metNames and move it to the model.metComps field.. This could 
    % be applied to the suffixies like "[*]" and "_*". Also the option to
    % merge all compartments to one compartment.
    % 
    % input
    %   model   
    %   mergeComps  merging all compartments to one compartment named System
    % output
    %   model       model without compartments in the model.mets and model.metNames


    if ~isfield(model, 'metComps')
        model.metComps= [];
    end

    % These two fields should be updated so that filled with regard to
    % the metabolite compartment patterns
    if ~isfield(model, 'comps')
        model.comps= {'s'};
    end
    if ~isfield(model, 'compNames')
        model.compNames= {'System'};
    end
    
    % finding the pattern of suffix in model.mets, same as model.metnames
    if model.mets{1}(end-1) == '_'
        pattern = '(_[a-zA-Z]+)$';
    elseif model.mets{1}(end) == ']'
        pattern = '(\[[a-zA-Z]+\]$)';
    else
        return
    end

    for i=1:size(model.mets, 1)
        metMatches= regexp(model.mets{i}, pattern, 'tokens', 'once');
        comp= regexp(metMatches{1}, '[a-z]+', 'match');

        if clearMets
            model.mets{i}= erase(model.mets{i}, metMatches{end});
        end

        if clearMetNames
            metNameMatches= regexp(model.metNames{i}, pattern, 'tokens', 'once');
            model.metNames{i}= erase(model.metNames{i}, metNameMatches{1});
        end
        compIndex= find(strcmp(comp, model.comps));
        if mergeComps
            model.metComps(i, 1)= 1;
        else
            model.metComps(i, 1)= compIndex;
        end
    end
    if mergeComps
        model.comps= {'s'};
        model.compNames= {'System'};
    end
end